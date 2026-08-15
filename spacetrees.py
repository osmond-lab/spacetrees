from scipy.optimize import minimize
import scipy.sparse as sp
import time
import numpy as np
import math
from tqdm import tqdm

def locate_ancestors(samples, ancestor_times,
                     processed_times, sample_locations, sigma=None,
                     BLUP=False, quiet=False):

    """
    Locate genetic ancestors given sample locations and shared times.

    processed_times: the per-tree data for a single locus, as loaded from
    process-times' *.times.pkl -- a list of M dicts, each with a 'subtrees'
    list of {sample_ids, shared_times, shared_times_logdet, shared_times_inv},
    one per isolated subtree (samples still sharing time since the cutoff T),
    plus that tree's 'branching_times' and 'logpcoal'.

    sigma: raw dispersal (and branching) rate, as written by estimate-dispersal
    (sdx[,sdy,rho],phi -- the branching rate phi is always the last entry).
    Used both to importance-weight each tree by how probable its branching
    times are relative to the standard coalescent, and (outside BLUP) as the
    dispersal covariance in the location likelihood. If None, trees are
    weighted equally (no importance sampling), and BLUP=True is required
    (the full-likelihood path needs an actual dispersal covariance).
    """

    if not quiet: print('\n%%%%%%%%%%%% locating ancestors with spacetrees %%%%%%%%%%%%')

    M = len(processed_times)
    try:
        n, d = sample_locations.shape
    except:
        n = len(sample_locations)
        d = 1
    if not quiet: print('number of trees per locus:',M,'\nnumber of samples:',n,'\nnumber of spatial dimensions:',d)
    if not quiet: print('samples:',samples,'\ntimes:',ancestor_times)

    # importance weights, from each tree's branching times and log coalescent-time probability;
    # with no sigma there's no branching rate to weight trees by, so weight them equally
    if sigma is None:
        log_weights = np.zeros(M)
    else:
        phi = sigma[-1] #branching rate is always the last entry
        sigma = _sds_rho_to_sigma(sigma[:-1]) #dispersal covariance matrix
        lbds = np.array([_log_birth_density(tree['branching_times'], phi, n) for tree in processed_times])
        lpcs = np.array([tree['logpcoal'] for tree in processed_times])
        log_weights = lbds - lpcs

    ancestor_locations = []
    samples_iter = samples if quiet else tqdm(samples)
    for sample in samples_iter:

        # a sample can fall in a different subtree of each tree, so precompute,
        # per tree, the subtree containing this sample and the quantities that
        # don't depend on the ancestor time
        precomp = []
        for tree in processed_times:
            subtree, idx = _get_focal_subtree(sample, tree['subtrees'])
            k = len(subtree['sample_ids'])
            st = subtree['shared_times'].astype(float) #shared times in subtree
            stmr = np.mean(st, axis=1) #average times in each row
            stm = np.mean(stmr) #average times in whole matrix
            stci = subtree['shared_times_inv']
            Tmat = np.identity(k) - [[1/k for _ in range(k)] for _ in range(k)]; Tmat = Tmat[:-1] #mean centering matrix
            locs = sample_locations[subtree['sample_ids']] #locations of samples in subtree
            loc_mean = np.mean(locs, axis=0) #mean location
            stcilc = np.matmul(stci, np.matmul(Tmat, locs)) #a product we will use
            precomp.append((idx, st, stmr, stm, stci, Tmat, loc_mean, stcilc))

        for anc_time in ancestor_times:

            # calculate likelihoods or mles over trees
            fs = []
            mles = []
            bvars = []
            for (idx, st, stmr, stm, stci, Tmat, loc_mean, stcilc), log_weight in zip(precomp, log_weights):

                at = _anc_times(st, anc_time, idx) #shared times between samples and ancestor of sample at time
                atc = np.matmul(Tmat, (at[:-1] - stmr)) #center this
                taac = at[-1] - 2*np.mean(at[:-1]) + stm #center shared times of ancestor with itself
                mle = loc_mean + np.matmul(atc.transpose(), stcilc) #most likely location

                # if getting best linear unbiased predictor we collect the mles (and variance) at each tree
                if BLUP:
                    mles.append(mle)
                    var = taac - np.matmul(np.matmul(atc.transpose(), stci), atc) #variance in loc (you can multiply by sigma later)
                    bvars.append(var)
                # and otherwise we get the full likelihood at each tree
                else:
                    var = (taac - np.matmul(np.matmul(atc.transpose(), stci), atc)) * sigma #variance in loc
                    fs.append(lambda x, mle=mle, var=var: _lognormpdf(x, mle, var)) #append likelihood (mle/var bound at append time, not lookup time)

            # combine information across trees
            if BLUP:
                blup = np.zeros(d)
                tot_weight = 0
                # weighted average of mles
                for mle, log_weight in zip(mles, log_weights):
                     blup += mle * np.exp(log_weight)
                     tot_weight += np.exp(log_weight)
                mle = blup/tot_weight
                # weighted average of variances
                blup_var = 0
                for bvar, log_weight in zip(bvars, log_weights):
                     blup_var += bvar * np.exp(log_weight)
                mle = np.append(mle, blup_var/tot_weight)
            else:
                # find min of negative of log of summed likelihoods (weighted by importance)
                def g(x):
                    return -_logsumexp([f(x) + log_weight for f,log_weight in zip(fs, log_weights)])
                x0 = sample_locations[sample]
                mle = minimize(g, x0=x0).x

            ancestor_locations.append([sample,anc_time] + [float(i) for i in mle])

    return ancestor_locations

def estimate_dispersal(locations, processed_times, method='L-BFGS-B',
                       important=True, quiet=False, BLUP=False):

    """
    Numerically estimate maximum likelihood dispersal rate (and possibly branching rate) given sample locations and shared times.

    processed_times: one entry per locus, each the per-tree data for that locus as
    loaded from process-times' *.times.pkl -- a list of M dicts, each with a
    'subtrees' list of {sample_ids, shared_times, shared_times_logdet,
    shared_times_inv} (one per isolated subtree since the cutoff T), plus
    that tree's 'branching_times' and 'logpcoal'.
    """

    if not quiet: print('\n%%%%%%%%%%%% inferring dispersal with spacetrees %%%%%%%%%%%%')

    L = len(processed_times)
    M = len(processed_times[0])
    try:
        n, d = locations.shape
    except:
        n = len(locations)
        d = 1
    if not quiet: print('number of loci:',L,'\nnumber of trees per locus:',M,'\nnumber of samples:',n,'\nnumber of spatial dimensions:',d,'\n')

    # find decent initial dispersal rate: a pooled average of per-subtree MLEs,
    # weighted by each subtree's degrees of freedom (k-1), since subtree sizes vary
    if not quiet: print('initializing dispersal rate...')
    guess = np.zeros((d,d))
    kw_total = 0
    loci_iter = processed_times if quiet else tqdm(processed_times)
    for trees in loci_iter: #loop over loci
        for tree in trees: #loop over trees
            for subtree in tree['subtrees']: #loop over subtrees
                k = len(subtree['sample_ids'])
                if k > 1: #need more than 1 sample in a subtree to estimate dispersal
                    Tmat = np.identity(k) - [[1/k for _ in range(k)] for _ in range(k)]; Tmat = Tmat[:-1] #mean centering matrix
                    loc_mc = np.matmul(Tmat, locations[subtree['sample_ids']]) #mean centered locations
                    mle = _mle_dispersal_tree(loc_mc, subtree['shared_times_inv'])
                    guess += (k - 1) * mle #weighted mle dispersal rate for subtree
                    kw_total += k - 1 #and weight
    guess = guess/kw_total #pooled mle dispersal rate over all subtrees, trees, and loci
    x0 = _sigma_to_sds_rho(guess) #convert initial dispersal rate to standard deviations and correlation, to feed into numerical search
    if BLUP:
        return x0 #best linear unbiased predictor (returned as sds and corr, like numerical search below)
    x0 = [i/2 for i in x0] #heuristic because the estimate seems to be a consistent overestimate
    if not quiet: print('initial dispersal rate:',x0)

    # initializing branching rate
    if important:
        phi0 = np.mean([np.log(n/(n-len(tree['branching_times'])+1))/tree['branching_times'][-1]
                         for trees in processed_times for tree in trees]) #initial guess at branching rate, from n(t)=n(0)e^(phi*t)
        if not quiet: print('initial branching rate:',phi0)
        scale_phi = x0[0]/phi0 #we will search for the value of phi*scale_phi that maximizes the likelihood (putting phi on same scale as dispersal accelarates search)
        x0.append(phi0*scale_phi)
    else:
        scale_phi = None

    # negative composite log likelihood ratio, as function of x
    f = _sum_mc(locations=locations, processed_times=processed_times, important=important, scale_phi=scale_phi)

    # impose bounds on parameters
    bnds = [(1e-6,None)] #sdx
    if d==2:
        bnds.append((1e-6,None)) #sdy
        bnds.append((-0.99,0.99)) #corr
    if important:
        bnds.append((1e-6,None)) #scaled phi

    # find mle
    if not quiet: print('\nsearching for maximum likelihood parameters...')
    t0 = time.time()
    m = minimize(f, x0=x0, bounds=bnds, method=method) #find MLE
    if not quiet: print(m)
    if not quiet: print('finding the max took', time.time()-t0, 'seconds')

    mle = m.x
    if important:
        mle[-1] = mle[-1]/scale_phi #unscale phi
    if not quiet:
        if important:
            sigma = _sds_rho_to_sigma(mle[:-1]) #convert to covariance matrix
            print('\nmaximum likelihood branching rate:',mle[-1])
        else:
            sigma = _sds_rho_to_sigma(mle)
        print('\nmaximum likelihood dispersal rate:\n',sigma)

    return mle 

def _get_focal_subtree(focal_sample, subtrees):

    """
    get the subtree containing focal_sample, and its index within that subtree
    (subtrees is a tree's list of {sample_ids, shared_times, ...} dicts)
    """

    for subtree in subtrees:
        ids = list(subtree['sample_ids'])
        if focal_sample in ids:
            return subtree, ids.index(focal_sample)

def _anc_times(shared_times, ancestor_time, sample):

    """
    get shared times with ancestor 
    """
    
    taa = shared_times[0,0] - ancestor_time #shared time of ancestor with itself 

    anc_times = [] 
    for t in shared_times[sample]:
        anc_times.append(min(t, taa)) # shared times between ancestor and each sample lineage

    anc_times.append(taa) #add shared time with itself
        
    return np.array(anc_times)

def _lognormpdf(x, mu, S):

    """
    Calculate log probability density of x, when x ~ N(mu,S)
    """

    norm_coeff = np.linalg.slogdet(S)[1] #just care about relative likelihood so drop the constant

    # term in exponential (times -2)
    err = x - mu #difference between mean and data
    if sp.issparse(S):
        numerator = spln.spsolve(S, err).T.dot(err) #use faster sparse methods if possible
    else:
        numerator = np.linalg.solve(S, err).T.dot(err) #just a fancy way of calculating err.T * S^-1  * err

    return -0.5 * (norm_coeff + numerator) #add the two terms together and multiply by -1/2

def _mle_dispersal_tree(locations, shared_times_inverted):

    """
    Maximum likelihood estimate of dispersal rate given locations and (inverted) shared times between lineages in a tree.
    """

    return np.matmul(np.matmul(np.transpose(locations), shared_times_inverted), locations) / len(locations)

def _sigma_to_sds_rho(sigma):

    """
    Convert 1x1 or 2x2 covariance matrix to sds and correlation
    """
    d = len(sigma)
 
    sdx = sigma[0,0]**0.5
    if d==1:
        return [sdx]
    elif d==2:
        sdy = sigma[1,1]**0.5
        rho = sigma[0,1]/(sdx * sdy) #note that small sdx and sdy will raise errors
        return [sdx, sdy, rho]

def _sum_mc(locations, processed_times, important=False, scale_phi=None):

    """
    Negative log composite likelihood of parameters x given the locations and shared times at all loci, trees, and subtrees, as function of x.
    """

    def sumf(x):

        # reformulate parameters
        if important:
            sigma = _sds_rho_to_sigma(x[:-1])
            phi = x[-1]/scale_phi
        else:
            sigma = _sds_rho_to_sigma(x)
            phi = None
        log_det_sigma = np.linalg.slogdet(sigma)[1] #log of determinant
        sigma_inverted = np.linalg.inv(sigma) #inverse

        # calculate negative log composite likelihood ratio
        # by subtracting log likelihood ratio at each locus
        g = 0
        for trees in processed_times: #loop over loci
            g -= _mc(locations=locations, trees=trees, sigma_inverted=sigma_inverted, log_det_sigma=log_det_sigma,
                     important=important, phi=phi)
        return g

    return sumf

def _sds_rho_to_sigma(x):

    """
    Convert sds and correlation to 1x1 or 2x2 covariance matrix
    """
    sdx = x[0]
    if len(x) == 1:
        sigma = np.array([[sdx**2]])
    else:
        sdy = x[1]
        rho = x[2]
        cov = sdx*sdy*rho
        sigma = np.array([[sdx**2, cov], [cov, sdy**2]])

    return sigma

def _mc(locations, trees, sigma_inverted, log_det_sigma, important=False, phi=None):

    """
    Monte Carlo estimate of log of likelihood ratio of the locations given parameters (sigma,phi) vs data given standard coalescent, for a given locus
    """

    # loop over trees at a locus, summing the log likelihood ratio over each tree's subtrees
    LLRs = [_log_likelihoodratio(locations=locations, tree=tree, sigma_inverted=sigma_inverted, log_det_sigma=log_det_sigma,
                                 important=important, phi=phi)
            for tree in trees]

    return _logsumexp(np.array(LLRs)) #sum likelihood ratios over trees then take log

def _logsumexp(a):

    """
    Take the log of a sum of exponentials without losing information.
    """

    a_max = np.max(a) #max element in list a
    tmp = np.exp(a - a_max) #now subtract off the max from each a before taking exponential (ie divide sum of exponentials by exp(a_max))
    s = np.sum(tmp) #and sum those up
    out = np.log(s) #and take log
    out += a_max  #and then add max element back on (ie multiply sum by exp(a_max), ie add log(exp(a_max)) to logged sum)

    return out

def _log_likelihoodratio(locations, tree, sigma_inverted, log_det_sigma, important=False, phi=None):

    """
    Log of likelihood ratio of parameters under branching brownian motion vs standard coalescent, summed over a tree's subtrees.
    """

    # log likelihood of dispersal rate, summed over subtrees
    LLR = 0
    n = 0 #total samples across subtrees (assumes centered-times convention: k-1 per subtree, plus 1)
    ksum = 0 #total degrees of freedom across subtrees
    for subtree in tree['subtrees']:
        k = len(subtree['sample_ids']) - 1 #after mean-centering
        n += k + 1
        if k > 0: #need more than 1 sample in a subtree to contribute a likelihood
            LLR += _location_loglikelihood(locations, subtree, sigma_inverted)
            ksum += k
    d,_ = sigma_inverted.shape
    if ksum > 0: LLR -= ksum/2 * (d*np.log(2*np.pi) + log_det_sigma) #can factor this out over subtrees

    if important and ksum > 0:
        # log probability of branching times given pure birth process with rate phi
        LLR += _log_birth_density(branching_times=tree['branching_times'], phi=phi, n=n)
        # log probability of coalescence times given standard coalescent (precalculated as parameter-independent)
        LLR -= tree['logpcoal']

    return LLR

def _location_loglikelihood(locations, subtree, sigma_inverted):

    """
    Log probability density of locations when locations ~ MVN(0,sigma_inverted*shared_times_inverted), for one subtree.
    """

    k = len(subtree['sample_ids'])
    Tmat = np.identity(k) - [[1/k for _ in range(k)] for _ in range(k)]; Tmat = Tmat[:-1] #mean centering matrix
    locs = np.matmul(Tmat, locations[subtree['sample_ids']]) #mean centered locations for this subtree
    locs_vector = np.transpose(locs).flatten() #make a vector

    # log of coefficient in front of exponential (times -2)
    d,_ = sigma_inverted.shape
    logcoeff = d*subtree['shared_times_logdet'] #just the part that depends on data

    # exponent (times -2)
    exponent = np.matmul(np.matmul(np.transpose(locs_vector), np.kron(sigma_inverted, subtree['shared_times_inv'])), locs_vector)

    return -0.5 * (logcoeff + exponent) #add the two terms together and multiply by -1/2

def _log_birth_density(branching_times, phi, n, condition_on_n=True):

    """
    Log probability of branching times given Yule process with branching rate phi.
    """

    T = branching_times[-1] #storing total time as last entry for convenience
    n0 = n - (len(branching_times) - 1) #initial number of lineages (number of samples minus number of coalescence events)
    
    logp = 0 #initialize log probability
    prevt = 0 #initialize time
    k = n0 #initialize number of lineages
    # probability of each branching time
    for t in branching_times[:-1]: #for each branching time t
        logp += np.log(k * phi) - k * phi *  (t - prevt) #log probability of waiting time t-prevt with k lineages
        prevt = t #update time
        k += 1 #update number of lineages

    # probability of no branching from most recent branching to T
    logp += - k * phi * (T - prevt)

    # condition on having n samples from n0 in time T
    if condition_on_n:
        logp -= np.log(math.comb(k - 1, k - n0) * (1 - np.exp(-phi * T))**(k - n0)) - phi * n0 * T # see page 234 of https://www.pitt.edu/~super7/19011-20001/19531.pdf for two different expressions

    return logp
