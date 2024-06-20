# pipeline to infer dispersal rates and locate genetic ancestors with spacetrees (Osmond & Coop 2024)

datadir = 'data/' #relative path to data directory
relatedir = 'relate' #path to your version of relate

# we start by assuming you have run Relate's EstimatePopulationSize to get the following files (see https://myersgroup.github.io/relate/modules.html#CoalescenceRate)
anc = datadir + 'test_chr{CHR}.anc' #name of anc files, with wildcard for chromosome (chr)
mut = datadir + 'test_chr{CHR}.mut' #name of mut files
dist = datadir + 'test_chr{CHR}.dist' #name of dist files, only needed if analyzing a subregion of the chromosome, which we are here so that filesizes are small
coal = datadir + 'test.coal' #name of coal file

# you also need the locations of every sample in the same order you gave those samples to relate
locations = datadir + 'test.locations' #if individuals are diploid you need to repeat each location twice

CHRS = [1] #list of chromosomes you have anc/mut files for
m = '1e-8' #estimated mutation rate
dispersal_loci = [1] #which loci to use to infer dispersal

# ---------------- get positions of all loci ------------------------------

loci = anc.replace('.anc','.loci') #filename for list of loci positions, just changing the suffix from 'anc' to 'loci'

rule loci_positions:
  input:
    mut=mut
  output:
    loci=loci
  threads: 1
  resources:
    runtime=15 #will be much shorter, but my server has 15m minumum 
  run:
    from utils import loci_positions
    loci_positions(input.mut, output.loci)

# the output is a space delimited file with the position of the first and last mutation at each locus, with each locus in a separate row

# ---------------- sample trees at a locus ------------------------------

# now we sample trees at a given locus
# more specifically, Relate fixes the topology at each locus and resamples the branch lengths
# see https://myersgroup.github.io/relate/modules.html#ReEstimateBranchLengths

newick = loci.replace('.loci','_{locus}locus_{M}M.newick')

rule sample_trees:
  input:
    loci=loci,
    anc=anc,
    mut=mut,
    dist=dist,
    coal=coal
  output:
    newick 
  params:
    prefix_in = anc.replace('.anc',''), #prefix of anc and mut files (relate searches for anc/mut files with this prefix)
    prefix_out = newick.replace('.newick','') #prefix of outfile (relate adds its own suffix)
  threads: 1
  resources:
    runtime=15
  shell:
    '''
    start=$( awk 'NR=={wildcards.locus} {{print $1}}' {input.loci} ) #position of first snp at locus
    stop=$( awk 'NR=={wildcards.locus} {{print $2}}' {input.loci} ) #position of last snp at locus
    module load gcc/11.3.0 #on my server i had to load the same version i built relate with to run it
    {relatedir}/scripts/SampleBranchLengths/SampleBranchLengths.sh \
                 -i {params.prefix_in} \
                 --dist {input.dist} \
                 --coal {input.coal} \
                 -o {params.prefix_out} \
                 -m {m} \
                 --format n \
                 --num_samples {wildcards.M} \
                 --first_bp $start \
                 --last_bp $stop \
                 --seed 1 
    '''

# ---------------- extract times from trees -----------------------------

# now we will extract the information we need from the trees, the shared times between each pair of lineages and the coalescence times

shared_times = newick.replace('.newick','.stss')
coal_times = newick.replace('.newick','.ctss')

rule extract_times:
  input:
    newick=newick 
  output:
    stss=shared_times,
    ctss=coal_times
  threads: 1
  resources:
    runtime=15
  run:
    # prevent numpy from using more than {threads} threads (useful for parallizing on my server)
    import os
    os.environ["OMP_NUM_THREADS"] = str(threads)
    os.environ["GOTO_NUM_THREADS"] = str(threads)
    os.environ["OPENBLAS_NUM_THREADS"] = str(threads)
    os.environ["MKL_NUM_THREADS"] = str(threads)
    os.environ["VECLIB_MAXIMUM_THREADS"] = str(threads)
    os.environ["NUMEXPR_NUM_THREADS"] = str(threads)

    # import tools
    import numpy as np
    from tsconvert import from_newick
    from utils import get_shared_times
    from tqdm import tqdm

    # open file of trees to read from
    with open(input.newick, mode='r') as f:
      
      # open files of shared and coalescent times to append to
      with open(output.stss, 'a') as stss:
        with open(output.ctss, 'a') as ctss:

          next(f) #skip header
          for line in tqdm(f, total=int(wildcards.M)): #for each tree sampled

            # import tree
            string = line.split()[4] #extract newick string only (Relate adds some info beforehand)
            ts = from_newick(string, min_edge_length=1e-6) #convert to tskit "tree sequence" (only one tree)
            tree = ts.first() #the only tree

            # get shared times
            samples = [int(ts.node(node).metadata['name']) for node in ts.samples()] #get index of each sample in list we gave to relate
            sample_order = np.argsort(samples) #get indices to put in ascending order
            ordered_samples = [ts.samples()[i] for i in sample_order] #order samples as in relate
            sts = get_shared_times(tree, ordered_samples) #get shared times between all pairs of samples, with rows and columns ordered as in relate
            stss.write(",".join([str(round(i)) for i in sts]) + '\n') #append as new line, rounding times to save space

            # get coalescent times
            cts = sorted([tree.time(i) for i in tree.nodes() if not tree.is_sample(i)]) #coalescence times, in ascending order
            ctss.write(",".join([str(round(i)) for i in cts]) + '\n') #append as new line, rounding times to save space

# ---------------- process shared times -----------------------------

# now we process the shared times, potentially cutting off the tree (to ignore distant past) and getting the exact quantities we need for inference

processed_shared_times = shared_times.replace('.stss','_{T}T.{end}')
ends = ['stss_logdet','stss_inv']

rule process_shared_times:
  input:
    stss = shared_times 
  output:
    expand(processed_shared_times, end=ends, allow_missing=True)
  threads: 1 
  resources:
    runtime=15
  run:
    # prevent numpy from using more than {threads} threads (useful for parallizing on my server)
    import os
    os.environ["OMP_NUM_THREADS"] = str(threads)
    os.environ["GOTO_NUM_THREADS"] = str(threads)
    os.environ["OPENBLAS_NUM_THREADS"] = str(threads)
    os.environ["MKL_NUM_THREADS"] = str(threads)
    os.environ["VECLIB_MAXIMUM_THREADS"] = str(threads)
    os.environ["NUMEXPR_NUM_THREADS"] = str(threads)

    # load tools
    import numpy as np
    from utils import chop_shared_times, center_shared_times
    from tqdm import tqdm

    # determine time cutoff
    T = wildcards.T #get time cutoff
    T = None if T=='None' else float(T) #format correctly

    # open file of shared times to read from
    with open(input.stss, 'r') as stss:

      # open files to write to
      with open(output[0], 'a') as stss_logdet:
        with open(output[1], 'a') as stss_inverted:
      
          # loop over trees at this locus 
          for sts in tqdm(stss, total=int(wildcards.M)):
      
            # load shared time matrix in vector form
            sts = np.fromstring(sts, dtype=float, sep=',') #convert from string to numpy array

            # chop
            sts = chop_shared_times(sts, T=T) #chop shared times to ignore history beyond T

            # convert to matrix form
            k = int((np.sqrt(1+8*len(sts))-1)/2) #get number of samples (from sum_i=0^k i = k(k+1)/2)
            sts_mat = np.zeros((k,k))
            sts_mat[np.triu_indices(k, k=0)] = sts #convert to numpy matrix
            sts_mat = sts_mat + sts_mat.T - np.diag(np.diag(sts_mat))      
            sts = sts_mat

            # center
            sts = center_shared_times(sts) 
      
            # determinant
            sts_logdet = np.linalg.slogdet(sts)[1] #magnitude of log determinant (ignore sign)
            stss_logdet.write(str(sts_logdet) + '\n') #append as new line 
      
            # inverse
            sts = np.linalg.inv(sts) #inverse
            sts = sts[np.triu_indices(k-1, k=0)] #convert to list
            stss_inverted.write(",".join([str(i) for i in sts]) + '\n') #append as new line

# ---------------- process coalescence times -----------------------------

# now we process the coalescence times, potentially cutting off the tree (to ignore distant past) and getting the exact quantities we need for inference

processed_coal_times = processed_shared_times 
ends = ['btss','lpcs']

rule process_coal_times:
  input:
    ctss = coal_times,
    coal = coal
  output:
    expand(processed_coal_times, end=ends, allow_missing=True)
  threads: 1 
  resources:
    runtime=15
  run:
    # prevent numpy from using more than {threads} threads (useful for parallizing on my server)
    import os
    os.environ["OMP_NUM_THREADS"] = str(threads)
    os.environ["GOTO_NUM_THREADS"] = str(threads)
    os.environ["OPENBLAS_NUM_THREADS"] = str(threads)
    os.environ["MKL_NUM_THREADS"] = str(threads)
    os.environ["VECLIB_MAXIMUM_THREADS"] = str(threads)
    os.environ["NUMEXPR_NUM_THREADS"] = str(threads)

    # load tools
    import numpy as np
    from tqdm import tqdm
    from utils import log_coal_density

    # determine time cutoff
    T = wildcards.T #get time cutoff
    T = None if T=='None' else float(T) #format correctly

    # effective population size
    epochs = np.genfromtxt(input.coal, skip_header=1, skip_footer=1) #time at which each epoch starts (and the final one ends)
    Nes = 0.5/np.genfromtxt(input.coal, skip_header=2)[2:] #effective population size during each epoch

    #open file of coalescence times to read from
    with open(input.ctss, 'r') as ctss:

      # open files to write to
      with open(output[0], 'a') as btss:
        with open(output[1], 'a') as lpcs:
        
          # loop over trees at this locus 
          for cts in tqdm(ctss, total=int(wildcards.M)):
            cts = np.fromstring(cts, dtype=float, sep=',') #convert from string to numpy array

            # branching times for branching Brownian model
            Tmax = cts[-1] #time to most recent common ancestor
            if T is not None and T < Tmax:
                Tmax = T #farthest time to go back to
            bts = Tmax - np.flip(cts) #branching times, in ascending order
            bts = bts[bts>0] #remove branching times at or before T
            bts = np.append(bts,Tmax) #append total time as last item      
            btss.write(",".join([str(i) for i in bts]) + '\n') #append as new line

            # probability of coalescence times under neutral coalescent
            lpc = log_coal_density(times=cts, Nes=Nes, epochs=epochs, T=Tmax) #log probability density of coalescence times
            lpcs.write(str(lpc) + '\n') #append as new line 

# ----------- composite dispersal rates ------------------------

# and now we bring in our processed times across chromosomes and loci to estimate a dispersal rate

composite_dispersal_rate = processed_shared_times.replace('_chr{CHR}','').replace('_{locus}locus','').replace('{end}','sigma')

rule composite_dispersal_rate:
  input:
    stss_logdet = expand(processed_shared_times, end=['stss_logdet'], CHR=CHRS, locus=dispersal_loci, allow_missing=True),
    stss_inv = expand(processed_shared_times, end=['stss_inv'], CHR=CHRS, locus=dispersal_loci, allow_missing=True),
    btss = expand(processed_coal_times, end=['btss'], CHR=CHRS, locus=dispersal_loci, allow_missing=True),
    lpcs = expand(processed_coal_times, end=['lpcs'], CHR=CHRS, locus=dispersal_loci, allow_missing=True),
    locations = locations
  output:
    sigma = composite_dispersal_rate
  threads: 1
  resources:
    runtime=15
  run:
    # prevent numpy from using more than {threads} threads (useful for parallizing on my server)
    import os
    os.environ["OMP_NUM_THREADS"] = str(threads)
    os.environ["GOTO_NUM_THREADS"] = str(threads)
    os.environ["OPENBLAS_NUM_THREADS"] = str(threads)
    os.environ["MKL_NUM_THREADS"] = str(threads)
    os.environ["VECLIB_MAXIMUM_THREADS"] = str(threads)
    os.environ["NUMEXPR_NUM_THREADS"] = str(threads)

    # load tools
    import numpy as np
    from tqdm import tqdm
    from spacetrees import mle_dispersal

    # load input data
    stss_logdet = [] #log determinants of chopped and centered shared times matrices
    for f in tqdm(input.stss_logdet):
      stss_logdet.append(np.loadtxt(f))
    stss_inv = [] #inverse of chopped and centered shared times matrices, in vector form
    for f in tqdm(input.stss_inv):
      sts_inv = np.loadtxt(f, delimiter=',') #list of vectorized matrices
      k = int((np.sqrt(1+8*len(sts_inv[0]))-1)/2) #get size of matrix (from sum_i=0^k i = k(k+1)/2)
      sts_inv_mat = [] #list of inverses in matrix form
      for st_inv in sts_inv:
        mat = np.zeros((k,k))
        mat[np.triu_indices(k, k=0)] = st_inv #convert to numpy matrix
        mat = mat + mat.T - np.diag(np.diag(mat))      
        sts_inv_mat.append(mat)
      stss_inv.append(sts_inv_mat)
    btss = [] #branching times
    for f in tqdm(input.btss):
      btss.append(np.loadtxt(f, delimiter=','))
    lpcs = [] #log probability of coalescence times
    for f in tqdm(input.lpcs):
      lpcs.append(np.loadtxt(f))
    locations = np.loadtxt(input.locations) #location of each sample

    # estimate dispersal rate
    def callbackF(x):
      '''print updates during numerical search'''
      print('{0: 3.6f}   {1: 3.6f}   {2: 3.6f}   {3: 3.6f}'.format(x[0], x[1], x[2], x[3]))
    sigma = mle_dispersal(locations=locations, shared_times_inverted=stss_inv, shared_times_logdet=stss_logdet,
                          branching_times=btss, logpcoals=lpcs,
                          callbackF=callbackF)
    with open(output.sigma, 'w') as f:
      f.write(','.join([str(i) for i in sigma])) #save variances and covariance

# ----------------------- locate ancestors -----------------------

# finally, we use our processed times and dispersal rate to locate the genetic ancestor at a particular locus for a particular sample and time
# TODO: it might be better to locate internal nodes of a tree

ancestor_locations = processed_shared_times.replace('.{end}','_{s}s_{t}t.locs')

rule locate_ancestors:
  input:
    stss = shared_times,
    stss_inv = processed_shared_times.replace('{end}','stss_inv'),
    btss = processed_coal_times.replace('{end}','btss'),
    lpcs = processed_coal_times.replace('{end}','lpcs'),
    locations = locations,
    sigma = composite_dispersal_rate
  output:
    ancestor_locations
  threads: 1
  resources:
    runtime=15
  run:
    # prevent numpy from using more than {threads} threads (useful for parallizing on my server)
    import os
    os.environ["OMP_NUM_THREADS"] = str(threads)
    os.environ["GOTO_NUM_THREADS"] = str(threads)
    os.environ["OPENBLAS_NUM_THREADS"] = str(threads)
    os.environ["MKL_NUM_THREADS"] = str(threads)
    os.environ["VECLIB_MAXIMUM_THREADS"] = str(threads)
    os.environ["NUMEXPR_NUM_THREADS"] = str(threads)

    # load tools
    import numpy as np
    from tqdm import tqdm
    from spacetrees import locate_ancestors, _log_birth_density, _sds_rho_to_sigma 
    from utils import chop_shared_times

    # load input data
    stss = np.loadtxt(input.stss, delimiter=',') #list of vectorized shared times matrices
    k = int((np.sqrt(1+8*len(stss[0]))-1)/2) #get size of matrix (from sum_i=0^k i = k(k+1)/2)
    T = wildcards.T #get time cutoff
    T = None if T=='None' else float(T) #format correctly
    stss_mat = [] #list of chopped shared times matrices in matrix form
    for sts in stss:
      sts = chop_shared_times(sts, T=T) #chop shared times to ignore history beyond T
      mat = np.zeros((k,k))
      mat[np.triu_indices(k, k=0)] = sts #convert to numpy matrix
      mat = mat + mat.T - np.diag(np.diag(mat))      
      stss_mat.append(mat)
    stss = stss_mat
    stss_inv = np.loadtxt(input.stss_inv, delimiter=',') #list of vectorized chopped centered inverted shared times matrices
    k = k-1 #get size of matrix
    stss_inv_mat = [] #list of chopped shared times matrices in matrix form
    for sts_inv in stss_inv:
      mat = np.zeros((k,k))
      mat[np.triu_indices(k, k=0)] = sts_inv #convert to numpy matrix
      mat = mat + mat.T - np.diag(np.diag(mat))      
      stss_inv_mat.append(mat)
    stss_inv = stss_inv_mat
    btss = np.loadtxt(input.btss, delimiter=',') #branching times
    lpcs = np.loadtxt(input.lpcs) #log probability of coalescence times
    locations = np.loadtxt(input.locations) #location of each sample
    sigma = np.loadtxt(input.sigma, delimiter=',') #mle dispersal rate and branching rate
    phi = sigma[-1] #branching rate
    sigma = _sds_rho_to_sigma(sigma[:-1]) #dispersal as covariance matrix

    # importance weights
    lbds = np.array([_log_birth_density(bts, phi, len(locations)) for bts in btss]) #log probability densities of birth times
    log_weights = lbds - lpcs #log importance weights

    # locate ancestors
    ancestor_locations = locate_ancestors(samples=[int(wildcards.s)], times=[float(wildcards.t)], 
                                          shared_times_chopped=stss, shared_times_chopped_centered_inverted=stss_inv, locations=locations, 
                                          sigma=sigma, log_weights=log_weights)
    print(ancestor_locations) 
