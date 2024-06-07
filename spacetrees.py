from scipy.optimize import minimize
import scipy.sparse as sp
import time
import numpy as np
import math
from tqdm import tqdm 

def locate_ancestors(ancestor_samples, ancestor_times, 
                     shared_times_chopped, samples, locations, 
                     log_weights=[0], sigma=1, x0_final=None, BLUP=False, BLUP_var=False):

    """
    Locate genetic ancestors given sample locations and shared times.
    """
    
    d = len(locations[0]) 

    all_ancestor_locations = []
    all_ancestor_variances = []
    # first loop over the samples we want to find ancestors of
#    for sample in tqdm(ancestor_samples):
    for sample in ancestor_samples:

        sts = []
        stmrs = []
        stms = []
        stcis = []
        locs_means = []
        stcilcs = []
        js = []
        # next loop over trees at this locus, getting the quantities we need that are independent of the age of the ancestor
        for stsc, smpls in zip(shared_times_chopped, samples):

            i,j = _get_focal_index(sample, smpls) #subtree and index of sample
            js.append(j)
            st = stsc[i].astype('float') #shared times in subtree (make sure float for faster computation)
            sts.append(st)
            stmr = np.mean(st, axis=1) #average times in each row
            stmrs.append(stmr)
            stm = np.mean(st) #average times in whole matrix
            stms.append(stm)
            n = len(st); 
            # note that if n=1, we get lots of empty matrices below, 
            # but the mle and var are calculated correctly (sample location and sigma*t respectively)

            Tmat = np.identity(n) - [[1/n for _ in range(n)] for _ in range(n)]; Tmat = Tmat[:-1] #mean centering matrix
            stc = np.matmul(Tmat, np.matmul(st, np.transpose(Tmat))) #center shared times matrix
            stci = np.linalg.inv(stc) #invert 
            stcis.append(stci)
            
            locs = locations[smpls[i].astype('int')] #locations of samples in subtree
            locs_mean = np.mean(locs, axis=0) #mean location
            locs_means.append(locs_mean)
            lc = np.matmul(Tmat, locs) #centered locations
            stcilc = np.matmul(stci, lc) #product used below
            stcilcs.append(stcilc)

        ancestor_locations = locations[sample]
        ancestor_variances = [0]
        # now we will loop over the times we want to locate the ancestor at
        for ancestor_time in ancestor_times:

            fs = []
            mles = []
            bvars = []
            # looping over trees at this locus again
            for st,stmr,stm,stci,lm,stcilc,j in zip(sts, stmrs, stms, stcis, locs_means, stcilcs, js):

                n = len(st); Tmat = np.identity(n) - [[1/n for _ in range(n)] for _ in range(n)]; Tmat = Tmat[:-1] #mean centering matrix

                at = _anc_times(st, ancestor_time, j) #shared times between samples and ancestor
                atc = np.matmul(Tmat, (at[:-1] - stmr)) #center this
                taac = at[-1] - 2*np.mean(at[:-1]) + stm #center shared times of ancestor with itself
                mle_loc = lm + np.matmul(atc.transpose(), stcilc) #mean loc
                if BLUP:
                    mles.append(mle_loc)
                    if BLUP_var:
                        var_loc = (taac - np.matmul(np.matmul(atc.transpose(), stci), atc)) * sigma #variance in loc
                        bvars.append(var_loc)
                else:
                    var_loc = (taac - np.matmul(np.matmul(atc.transpose(), stci), atc)) * sigma #variance in loc
                    fs.append(lambda x: _lognormpdf(x, mle_loc, var_loc)) #append likelihood

            if BLUP:
                blup = np.zeros(d) 
                tot_weight = 0
                # weighted average of mles
                for mle, log_weight in zip(mles, log_weights):
                     blup += mle * np.exp(log_weight)
                     tot_weight += np.exp(log_weight)
                blup = blup/tot_weight
                ancestor_locations = np.vstack([ancestor_locations,blup])
                # weighted average of variances
                blup_var = 0
                for bvar, log_weight in zip(bvars, log_weights):
                     blup_var += bvar * np.exp(log_weight)
                blup_var = blup_var/tot_weight
                ancestor_variances = np.vstack([ancestor_variances,blup_var])

            else:
                # find min of negative of log of summed likelihoods (weighted by importance)
                def g(x): 
                    return -_logsumexp([f(x) + log_weight for f,log_weight in zip(fs, log_weights)])
                x0 = locations[sample] 
                if x0_final is not None:
                    x0 = x0 + (x0_final - x0)*ancestor_time/ancestor_times[-1] #make a linear guess
                mle = minimize(g, x0=x0).x
                ancestor_locations = np.vstack([ancestor_locations,mle])

        all_ancestor_locations.append(ancestor_locations)
        if BLUP_var:
            all_ancestor_variances.append(ancestor_variances)

    if BLUP_var:
        return np.array(all_ancestor_locations), np.array(all_ancestor_variances)  
        
    return np.array(all_ancestor_locations)

def mle_dispersal(locations, shared_times_inverted, samples, log_det_shared_times=None, 
                  sigma0=None, bnds=None, method='L-BFGS-B', callbackF=None,
                  important=False, branching_times=None, phi0=None, scale_phi=None, logpcoals=None,
                  quiet=False, BLUP=False):

    """
    Numerically estimate maximum likelihood dispersal rate (and possibly branching rate) given sample locations and shared times.
    """

    L = len(shared_times_inverted)
    M = len(shared_times_inverted[0])
    try: 
        n, d = locations.shape
    except:
        n = len(locations)
        d = 1
    if not quiet: print('number of loci:',L,'\nnumber of trees per locus:',M,'\nnumber of samples:',n,'\nnumber of spatial dimensions:',d)

    # prepare locations
    if not quiet: print('\npreparing locations')
    if sigma0 is not None and not BLUP:
        locsss = []
        for smplss in tqdm(samples): #loci
            locss = []
            for smpls in smplss: #trees
                locs = []
                for smpl in smpls: #subtrees
                    k = len(smpl)
                    loc = np.array([locations[i] for i in smpl]) #locations of that subtree
                    Tmat = np.identity(k) - [[1/k for _ in range(k)] for _ in range(k)]; Tmat = Tmat[0:-1] #mean centering matrix
                    loc_mc = np.matmul(Tmat, loc) #mean centered locations
                    loc_mc_vec = np.transpose(loc_mc).flatten() #make a vector
                    locs.append(loc_mc_vec)
                locss.append(locs)
            locsss.append(locss)
    # and potentially find decent initial dispersal rate
    else:
        if not quiet: print('and finding best linear unbiased predictor (BLUP) of dispersal rate')
        locsss = []
        #kvsum = np.zeros((d,d))
        #ksum = 0
        blup = np.zeros((d,d))
        for stss, smplss in tqdm(zip(shared_times_inverted, samples), total=L): #loci
            locss = []
            #kvisum = np.zeros((d,d))
            #kisum = 0
            for sts, smpls in zip(stss, smplss): #trees
                locs = []
                kvijsum = np.zeros((d,d))
                kijsum = 0
                for st, smpl in zip(sts, smpls): #subtrees
                    k = len(smpl) #number of samples in subtree before mean centering
                    loc = np.array([locations[i] for i in smpl]) #locations of that subtree
                    Tmat = np.identity(k) - [[1/k for _ in range(k)] for _ in range(k)]; Tmat = Tmat[0:-1] #mean centering matrix
                    loc_mc = np.matmul(Tmat, loc) #mean centered locations
                    loc_mc_vec = np.transpose(loc_mc).flatten() #make a vector
                    locs.append(loc_mc_vec)
                    if k>1: #need more than 1 sample in a subtree to estimate dispersal
                        mle = _mle_dispersal_tree(loc_mc, st.astype('float')) #ensure matrix is float for faster computation and broadcasting below
                        kvijsum += (k-1)*mle #add weighted mle dispersal rate for subtree 
                        kijsum += k-1 #and weight
                locss.append(locs)
                #if not BLUP:
                #if kisum == 0: #just use first tree at each locus to initialize sigma
                #kvisum += kvijsum #add weighted mle dispersal rate for tree
                #kisum += kijsum
                #else:
                blup += kvijsum/kijsum #add mle for this tree
            locsss.append(locss)
            #if not BLUP:
            #kvsum += kvisum #add weighted mle dispersal rate for locus
            #ksum += kisum
            
        #if not BLUP:
        #sigma0 = kvsum/ksum #this is the mle dispersal rate over loci and subtrees (using just the first tree at each locus)
        #if not quiet: print('initial dispersal rate:\n',sigma0)

    blup = blup/(L*M) #avg mle over all trees and loci (note that we can avg over all trees and loci simultaneously because same number of trees at every locus)
    if not quiet: print('BLUP dispersal rate:\n',blup)
    x0 = _sigma_to_sds_rho(blup) #convert initial dispersal rate to standard deviations and correlation, to feed into numerical search
    if BLUP:
        return x0 #best linear unbiased predictor (returned as sds and corr, like numerical search below)
    

    # initializing branching rate
    if important:
        if phi0 is None:
            phi0 = np.mean([np.log(n/len(smpls))/bts[-1] for smplss,btss in zip(samples,branching_times) for smpls,bts in zip(smplss,btss)]) #initial guess at branching rate, from n(t)=n(0)e^(phi*t)
            if not quiet: print('initial branching rate:',phi0) 
        if scale_phi is None:
            scale_phi = x0[0]/phi0 #we will search for the value of phi*scale_phi that maximizes the likelihood (putting phi on same scale as dispersal accelarates search) 
        x0.append(phi0*scale_phi)
        
    # negative composite log likelihood ratio, as function of x
    f = _sum_mc(locations=locsss, shared_times_inverted=shared_times_inverted, log_det_shared_times=log_det_shared_times,
                important=important, branching_times=branching_times, scale_phi=scale_phi, logpcoals=logpcoals)

    # impose bounds on parameters
    if bnds is None:
        bnds = [(1e-6,None)] #sdx
        if d==2:
            bnds.append((1e-6,None)) #sdy
            bnds.append((-0.99,0.99)) #corr
        if important:
            bnds.append((1e-6,None)) #scaled phi

    # find mle
    if not quiet: print('\nsearching for maximum likelihood parameters...')
    if callbackF is not None: callbackF(x0)
    t0 = time.time()
    m = minimize(f, x0=x0, bounds=bnds, method=method, callback=callbackF) #find MLE
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

def _get_focal_index(focal_node, listoflists):

    """
    get the subtree and index within that subtree for focal_node (listoflists here is list of samples for each subtree)
    """

    for i,j in enumerate(listoflists):
        if focal_node in j:
            n = i
            for k,l in enumerate(j):
                if focal_node == l:
                    m = k
    return n,m

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

def _sum_mc(locations, shared_times_inverted, log_det_shared_times,
            important=False, branching_times=None, scale_phi=None, logpcoals=None):

    """
    Negative log composite likelihood of parameters x given the locations and shared times at all loci and subtrees, as function of x.
    """

    if not important:
        L = len(log_det_shared_times) #number of loci
        branching_times = [None for _ in range(L)]
        logpcoals = branching_times

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
        for locs, sti, ldst, bts, lpcs in zip(locations, shared_times_inverted, log_det_shared_times, branching_times, logpcoals): #loop over loci
            g -= _mc(locations=locs, shared_times_inverted=sti, log_det_shared_times=ldst,
                     sigma_inverted=sigma_inverted, log_det_sigma=log_det_sigma,
                     important=important, branching_times=bts, phi=phi, logpcoals=lpcs)
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

def _mc(locations, shared_times_inverted, log_det_shared_times, sigma_inverted, log_det_sigma,
        important=False, branching_times=None, phi=None, logpcoals=None):

    """
    Monte Carlo estimate of log of likelihood ratio of the locations given parameters (sigma,phi) vs data given standard coalescent, for a given locus
    """

    LLRs = [] #log likelihood ratios at each tree

    # loop over trees at a locus
    if important:
        for locs, sti, ldst, bts, lpc in zip(locations, shared_times_inverted, log_det_shared_times, branching_times, logpcoals):
            LLRs.append(_log_likelihoodratio(locations=locs, shared_times_inverted=sti, log_det_shared_times=ldst,
                                             sigma_inverted=sigma_inverted, log_det_sigma=log_det_sigma, 
                                             important=important, branching_times=bts, phi=phi, logpcoals=lpc))

    else:
        for locs, sti, ldst in zip(locations, shared_times_inverted, log_det_shared_times):
            LLRs.append(_log_likelihoodratio(locations=locs, shared_times_inverted=sti, log_det_shared_times=ldst,
                                             sigma_inverted=sigma_inverted, log_det_sigma=log_det_sigma,
                                             important=important))
    
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

def _log_likelihoodratio(locations, shared_times_inverted, log_det_shared_times, sigma_inverted, log_det_sigma,
                         important=False, branching_times=None, phi=None, logpcoals=None):

    """ 
    Log of likelihood ratio of parameters under branching brownian motion vs standard coalescent.
    """
  
    # log likelihood of dispersal rate
    LLR = 0
    n = 0
    ksum = 0
    for locs, sts, ldst in zip(locations, shared_times_inverted, log_det_shared_times): #loop over subtrees
        k = len(sts); n += k + 1 #the plus 1 assumes times are mean centered
        if k>0:
            LLR += _location_loglikelihood(locs, sts.astype('float'), ldst, sigma_inverted)
            ksum += k
    d,_ = sigma_inverted.shape
    if ksum>0: LLR -= ksum/2 * (d*np.log(2*np.pi) + log_det_sigma)  #can factor this out over subtrees

    if important and ksum>0:
        # log probability of branching times given pure birth process with rate phi
        LLR += _log_birth_density(branching_times=branching_times, phi=phi, n=n) 
        # log probability of coalescence times given standard coalescent (precalculated as parameter-independent)
        LLR -= logpcoals
    
    return LLR

def _location_loglikelihood(locations, shared_times_inverted, log_det_shared_times, sigma_inverted):

    """
    Log probability density of locations when locations ~ MVN(0,sigma_inverted*shared_times_inverted).
    """
    
    # log of coefficient in front of exponential (times -2)
    d,_ = sigma_inverted.shape
    logcoeff = d*log_det_shared_times #just the part that depends on data

    # exponent (times -2)
    exponent = np.matmul(np.matmul(np.transpose(locations), np.kron(sigma_inverted, shared_times_inverted)), locations)   

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

def _log_coal_density(times, Nes, epochs=None, tCutoff=None):

    """
    log probability of coalescent times under standard neutral/panmictic coalescent
    """

    if epochs is None and len(Nes) == 1:
        epochs = [0, max(times)] #one big epoch
        Nes = [Nes[0], Nes[0]] #repeat the effective population size so same length as epochs 

    logp = 0 #initialize log probability
    prevt = 0 #initialize time
    prevLambda = 0 #initialize coalescent intensity
    k = len(times) + 1 #number of samples
    if tCutoff is not None:
        times = times[times < tCutoff] #ignore old times
    myIntensityMemos = _coal_intensity_memos(epochs, Nes) #intensities up to end of each epoch

    # probability of each coalescence time
    for t in times: #for each coalescence time t
        kchoose2 = k * (k - 1) / 2 #binomial coefficient
        Lambda = _coal_intensity_using_memos(t, epochs, myIntensityMemos, Nes) #coalescent intensity up to time t
        ie = np.digitize(np.array([t]), epochs) #epoch at the time of coalescence
        logpk = np.log(kchoose2 * 1 / (2 * Nes[ie])) - kchoose2 * (Lambda - prevLambda) #log probability (waiting times are time-inhomogeneous exponentially distributed)
        logp += logpk #add log probability
        prevt = t #update time
        prevLambda = Lambda #update intensity
        k -= 1 #update number of lineages

    # now add the probability of lineages not coalescing by tCutoff
    if k > 1 and tCutoff is not None: #if we have more than one lineage remaining
        kchoose2 = k * (k - 1) / 2 #binomial coefficient
        Lambda = _coal_intensity_using_memos(tCutoff, epochs, myIntensityMemos, Nes) #coalescent intensity up to time tCutoff
        logPk = - kchoose2 * (Lambda - prevLambda) #log probability of no coalescence
        logp += logPk #add log probability

    return logp[0] #FIX: extra dimn introduced somewhere

def _coal_intensity_using_memos(t, epochs, intensityMemos, Nes):

    """
    add coal intensity up to time t
    """

    iEpoch = int(np.digitize(np.array([t]), epochs)[0] - 1) #epoch 
    t1 = epochs[iEpoch] #time at which the previous epoch ended
    Lambda = intensityMemos[iEpoch] #intensity up to end of previous epoch
    Lambda += 1 / (2 * Nes[iEpoch]) * (t - t1) #add intensity for time in current epoch
    return Lambda

def _coal_intensity_memos(epochs, Nes):

    """
    coalescence intensity up to the end of each epoch
    """

    Lambda = np.zeros(len(epochs))
    for ie in range(1, len(epochs)):
        t0 = epochs[ie - 1] #start time
        t1 = epochs[ie] #end time
        Lambda[ie] = (t1 - t0) #elapsed time
        Lambda[ie] *= 1 / (2 * Nes[ie - 1]) #multiply by coalescence intensity
        Lambda[ie] += Lambda[ie - 1] #add previous intensity

    return Lambda
