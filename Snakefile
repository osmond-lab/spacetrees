# pipeline to infer dispersal rates and locate genetic ancestors with spacetrees (Osmond & Coop 2024)

datadir = 'data/' #relative path to data directory

# we start by assuming you have run Relate's EstimatePopulationSize to get the following files (see https://myersgroup.github.io/relate/modules.html#CoalescenceRate)
anc = datadir + 'test_chr{CHR}.anc' #name of anc files, with wildcard for chromosome (chr)
mut = datadir + 'test_chr{CHR}.mut' #name of mut files
coal = datadir + 'test.coal' #name of coal file

CHRS = [1] #list of chromosomes you have anc/mut files for
m = '1e-8' #estimated mutation rate

# ---------------- choose loci ------------------------------

# given these files, the first step is to choose which trees (herafter loci) to use for inference
# we want to sample loci sparsely so that the trees at any two loci are roughly independent
# we do this by choosing every nth locus, but you could imagine doing this differently (every nth basepair, 100 evenly spaced loci, etc)

skip = 10 #every 10th locus. in real datasets you will likely want this to be greater. we used 100 in our paper. the tradeoff is indepdence vs info.
loci = anc.replace('.anc','_10skip.loci') #filename for list of sampled loci, just changing the suffix from 'anc' to 'loci'

# we will run this rule once for each chromosome
rule choose_loci:
  input:
    anc=anc, 
    mut=mut, 
    coal=coal
  output:
    loci=loci
  threads: 1
  resources:
    runtime=15 #will be much shorter, but my server has 15m minumum 
  run:
    from utils import choose_loci
    choose_loci(input.anc, input.mut, skip, output.loci)

# the output is a space delimited file with the position of the first and last mutation in each of the chosen loci, with each locus in a separate row

# ---------------- sample trees at each chosen locus ------------------------------

# now we sample trees at each chose locus
# more specifically, Relate fixes the topology at each locus and resamples the branch lengths
# see https://myersgroup.github.io/relate/modules.html#ReEstimateBranchLengths

relatedir = '~/Applications/relate_v1.2.1_MacOSX_Intel' #directory where your version of relate is

M = 10 #number of trees to sample at each locus
newick = loci.replace('.loci','_{locus}locus_{M}M.newick') #will output all M trees at a locus as a single Newick file, each tree being a new line

rule sample_trees:
  input:
    loci=loci,
    anc=anc,
    mut=mut,
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
    # module load gcc/11.3.0 #on my server i had to load the same version i built relate with to run it
    {relatedir}/scripts/SampleBranchLengths/SampleBranchLengths.sh \
                 -i {params.prefix_in} \
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

shared_times = newick.replace('.newick','_sts.npy')
coal_times = newick.replace('.newick','_cts.npy')

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
    # tame numpy if needed (required on my server)
    #import os
    #os.environ["OMP_NUM_THREADS"] = str(threads)
    #os.environ["GOTO_NUM_THREADS"] = str(threads)
    #os.environ["OPENBLAS_NUM_THREADS"] = str(threads)
    #os.environ["MKL_NUM_THREADS"] = str(threads)
    #os.environ["VECLIB_MAXIMUM_THREADS"] = str(threads)
    #os.environ["NUMEXPR_NUM_THREADS"] = str(threads)
    import numpy as np

    from npy_append_array import NpyAppendArray
    from tsconvert import from_newick
    from utils import get_shared_times
    from tqdm import tqdm

    with open(input.newick, mode='r') as f:
      
      # will save times as numpy arrays, but want to save tree by tree to reduce memory requirement, so appending
      with NpyAppendArray(output.stss, delete_if_exists=True) as stss:
        with NpyAppendArray(output.ctss, delete_if_exists=True) as ctss:

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
            stss.append(np.array([sts])) #append to numpy file

            # get coalescent times
            cts = sorted([tree.time(i) for i in tree.nodes() if not tree.is_sample(i)]) #coalescence times, in ascending order
            ctss.append(np.array([cts])) #append to numpy file
