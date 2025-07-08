# pipeline to infer dispersal rates and locate genetic ancestors with spacetrees (Osmond & Coop 2024)

datadir = 'data/' #relative path to data directory
relatedir = 'relate' #path to your version of relate

# we start by assuming you have run Relate's EstimatePopulationSize to get the following files (see https://myersgroup.github.io/relate/modules.html#CoalescenceRate)
#prefix = 'test' #only contemporary samples
prefix = 'SGDP_contemporary_only_chr2' #contemporary and ancient samples
# anc = datadir + 'SGDP_aDNA_new_filtered_chr2.anc' #name of anc files, with wildcard for chromosome (chr)
# mut = datadir + 'SGDP_aDNA_new_filtered_chr2.mut' #name of mut files
anc = datadir + "SGDP_contemporary_only_chr2.anc"
mut = datadir + "SGDP_contemporary_only_chr2.mut"
#dist = datadir + prefix + '_chr{CHR}.dist' #name of dist files, only needed if analyzing a subregion of the chromosome, which we are here so that filesizes are small
coal = datadir + 'SGDP_v1_annot_ne.coal' #name of coal file

# you also need the locations of every sample in the same order you gave those samples to relate
#locations = datadir + prefix + '.locations' #if individuals are diploid you need to repeat each location twice
locations = datadir + "SGDP_contemporary_only.locations"

CHRS = [2] #list of chromosomes you have anc/mut files for
m = '4e-9' #estimated mutation rate
#dispersal_loci = [1] #which loci to use to infer dispersal
#ancestor_loci = dispersal_loci #which loci to locate ancestors at
#ancestor_times = [10,100,1000] #times in the past to locate ancestors at (if t='All')
# from config_contemporary import dispersal_loci
#from config_contemporary import ancestor_loci
#from config_contemporary import ancestor_times

dispersal_loci = [523, 1047, 1571, 2094, 2618, 3142, 3665, 4189, 4713, 5236, 5760, 6284, 6807, 7331, 7855, 8378, 8902, 9426, 9950, 10473, 10997, 11521, 12044, 12568, 13092, 13615, 14139, 14663, 15186, 15710, 16234, 16757, 17281, 17805, 18329, 18852, 19376, 19900, 20423, 20947, 21471, 21994, 22518, 23042, 23565, 24089, 24613, 25136, 25660, 26184, 26708, 27231, 27755, 28279, 28802, 29326, 29850, 30373, 30897, 31421, 31944, 32468, 32992, 33515, 34039, 34563, 35087, 35610, 36134, 36658, 37181, 37705, 38229, 38752, 39276, 39800, 40323, 40847, 41371, 41894, 42418, 42942, 43466, 43989, 44513, 45037, 45560, 46084, 46608, 47131, 47655, 48179, 48702, 49226, 49750, 50273, 50797, 51321]
ancestor_loci = [51, 103, 155, 207, 259, 311, 363, 415, 467, 518, 570, 622, 674, 726, 778, 830, 882, 934, 986, 1037, 1089, 1141, 1193, 1245, 1297, 1349, 1401, 1453, 1505, 1556, 1608, 1660, 1712, 1764, 1816, 1868, 1920, 1972, 2023, 2075, 2127, 2179, 2231, 2283, 2335, 2387, 2439, 2491, 2542, 2594, 2646, 2698, 2750, 2802, 2854, 2906, 2958, 3010, 3061, 3113, 3165, 3217, 3269, 3321, 3373, 3425, 3477, 3528, 3580, 3632, 3684, 3736, 3788, 3840, 3892, 3944, 3996, 4047, 4099, 4151, 4203, 4255, 4307, 4359, 4411, 4463, 4515, 4566, 4618, 4670, 4722, 4774, 4826, 4878, 4930, 4982, 5033, 5085, 5137, 5189, 5241, 5293, 5345, 5397, 5449, 5501, 5552, 5604, 5656, 5708, 5760, 5812, 5864, 5916, 5968, 6020, 6071, 6123, 6175, 6227, 6279, 6331, 6383, 6435, 6487, 6539, 6590, 6642, 6694, 6746, 6798, 6850, 6902, 6954, 7006, 7057, 7109, 7161, 7213, 7265, 7317, 7369, 7421, 7473, 7525, 7576, 7628, 7680, 7732, 7784, 7836, 7888, 7940, 7992, 8044, 8095, 8147, 8199, 8251, 8303, 8355, 8407, 8459, 8511, 8562, 8614, 8666, 8718, 8770, 8822, 8874, 8926, 8978, 9030, 9081, 9133, 9185, 9237, 9289, 9341, 9393, 9445, 9497, 9549, 9600, 9652, 9704, 9756, 9808, 9860, 9912, 9964, 10016, 10067, 10119, 10171, 10223, 10275, 10327, 10379, 10431, 10483, 10535, 10586, 10638, 10690, 10742, 10794, 10846, 10898, 10950, 11002, 11054, 11105, 11157, 11209, 11261, 11313, 11365, 11417, 11469, 11521, 11573, 11624, 11676, 11728, 11780, 11832, 11884, 11936, 11988, 12040, 12091, 12143, 12195, 12247, 12299, 12351, 12403, 12455, 12507, 12559, 12610, 12662, 12714, 12766, 12818, 12870, 12922, 12974, 13026, 13078, 13129, 13181, 13233, 13285, 13337, 13389, 13441, 13493, 13545, 13596, 13648, 13700, 13752, 13804, 13856, 13908, 13960, 14012, 14064, 14115, 14167, 14219, 14271, 14323, 14375, 14427, 14479, 14531, 14583, 14634, 14686, 14738, 14790, 14842, 14894, 14946, 14998, 15050, 15101, 15153, 15205, 15257, 15309, 15361, 15413, 15465, 15517, 15569, 15620, 15672, 15724, 15776, 15828, 15880, 15932, 15984, 16036, 16088, 16139, 16191, 16243, 16295, 16347, 16399, 16451, 16503, 16555, 16607, 16658, 16710, 16762, 16814, 16866, 16918, 16970, 17022, 17074, 17125, 17177, 17229, 17281, 17333, 17385, 17437, 17489, 17541, 17593, 17644, 17696, 17748, 17800, 17852, 17904, 17956, 18008, 18060, 18112, 18163, 18215, 18267, 18319, 18371, 18423, 18475, 18527, 18579, 18630, 18682, 18734, 18786, 18838, 18890, 18942, 18994, 19046, 19098, 19149, 19201, 19253, 19305, 19357, 19409, 19461, 19513, 19565, 19617, 19668, 19720, 19772, 19824, 19876, 19928, 19980, 20032, 20084, 20135, 20187, 20239, 20291, 20343, 20395, 20447, 20499, 20551, 20603, 20654, 20706, 20758, 20810, 20862, 20914, 20966, 21018, 21070, 21122, 21173, 21225, 21277, 21329, 21381, 21433, 21485, 21537, 21589, 21641, 21692, 21744, 21796, 21848, 21900, 21952, 22004, 22056, 22108, 22159, 22211, 22263, 22315, 22367, 22419, 22471, 22523, 22575, 22627, 22678, 22730, 22782, 22834, 22886, 22938, 22990, 23042, 23094, 23146, 23197, 23249, 23301, 23353, 23405, 23457, 23509, 23561, 23613, 23664, 23716, 23768, 23820, 23872, 23924, 23976, 24028, 24080, 24132, 24183, 24235, 24287, 24339, 24391, 24443, 24495, 24547, 24599, 24651, 24702, 24754, 24806, 24858, 24910, 24962, 25014, 25066, 25118, 25169, 25221, 25273, 25325, 25377, 25429, 25481, 25533, 25585, 25637, 25688, 25740, 25792, 25844, 25896, 25948, 26000, 26052, 26104, 26156, 26207, 26259, 26311, 26363, 26415, 26467, 26519, 26571, 26623, 26675, 26726, 26778, 26830, 26882, 26934, 26986, 27038, 27090, 27142, 27193, 27245, 27297, 27349, 27401, 27453, 27505, 27557, 27609, 27661, 27712, 27764, 27816, 27868, 27920, 27972, 28024, 28076, 28128, 28180, 28231, 28283, 28335, 28387, 28439, 28491, 28543, 28595, 28647, 28698, 28750, 28802, 28854, 28906, 28958, 29010, 29062, 29114, 29166, 29217, 29269, 29321, 29373, 29425, 29477, 29529, 29581, 29633, 29685, 29736, 29788, 29840, 29892, 29944, 29996, 30048, 30100, 30152, 30203, 30255, 30307, 30359, 30411, 30463, 30515, 30567, 30619, 30671, 30722, 30774, 30826, 30878, 30930, 30982, 31034, 31086, 31138, 31190, 31241, 31293, 31345, 31397, 31449, 31501, 31553, 31605, 31657, 31709, 31760, 31812, 31864, 31916, 31968, 32020, 32072, 32124, 32176, 32227, 32279, 32331, 32383, 32435, 32487, 32539, 32591, 32643, 32695, 32746, 32798, 32850, 32902, 32954, 33006, 33058, 33110, 33162, 33214, 33265, 33317, 33369, 33421, 33473, 33525, 33577, 33629, 33681, 33732, 33784, 33836, 33888, 33940, 33992, 34044, 34096, 34148, 34200, 34251, 34303, 34355, 34407, 34459, 34511, 34563, 34615, 34667, 34719, 34770, 34822, 34874, 34926, 34978, 35030, 35082, 35134, 35186, 35237, 35289, 35341, 35393, 35445, 35497, 35549, 35601, 35653, 35705, 35756, 35808, 35860, 35912, 35964, 36016, 36068, 36120, 36172, 36224, 36275, 36327, 36379, 36431, 36483, 36535, 36587, 36639, 36691, 36743, 36794, 36846, 36898, 36950, 37002, 37054, 37106, 37158, 37210, 37261, 37313, 37365, 37417, 37469, 37521, 37573, 37625, 37677, 37729, 37780, 37832, 37884, 37936, 37988, 38040, 38092, 38144, 38196, 38248, 38299, 38351, 38403, 38455, 38507, 38559, 38611, 38663, 38715, 38766, 38818, 38870, 38922, 38974, 39026, 39078, 39130, 39182, 39234, 39285, 39337, 39389, 39441, 39493, 39545, 39597, 39649, 39701, 39753, 39804, 39856, 39908, 39960, 40012, 40064, 40116, 40168, 40220, 40271, 40323, 40375, 40427, 40479, 40531, 40583, 40635, 40687, 40739, 40790, 40842, 40894, 40946, 40998, 41050, 41102, 41154, 41206, 41258, 41309, 41361, 41413, 41465, 41517, 41569, 41621, 41673, 41725, 41777, 41828, 41880, 41932, 41984, 42036, 42088, 42140, 42192, 42244, 42295, 42347, 42399, 42451, 42503, 42555, 42607, 42659, 42711, 42763, 42814, 42866, 42918, 42970, 43022, 43074, 43126, 43178, 43230, 43282, 43333, 43385, 43437, 43489, 43541, 43593, 43645, 43697, 43749, 43800, 43852, 43904, 43956, 44008, 44060, 44112, 44164, 44216, 44268, 44319, 44371, 44423, 44475, 44527, 44579, 44631, 44683, 44735, 44787, 44838, 44890, 44942, 44994, 45046, 45098, 45150, 45202, 45254, 45305, 45357, 45409, 45461, 45513, 45565, 45617, 45669, 45721, 45773, 45824, 45876, 45928, 45980, 46032, 46084, 46136, 46188, 46240, 46292, 46343, 46395, 46447, 46499, 46551, 46603, 46655, 46707, 46759, 46811, 46862, 46914, 46966, 47018, 47070, 47122, 47174, 47226, 47278, 47329, 47381, 47433, 47485, 47537, 47589, 47641, 47693, 47745, 47797, 47848, 47900, 47952, 48004, 48056, 48108, 48160, 48212, 48264, 48316, 48367, 48419, 48471, 48523, 48575, 48627, 48679, 48731, 48783, 48834, 48886, 48938, 48990, 49042, 49094, 49146, 49198, 49250, 49302, 49353, 49405, 49457, 49509, 49561, 49613, 49665, 49717, 49769, 49821, 49872, 49924, 49976, 50028, 50080, 50132, 50184, 50236, 50288, 50339, 50391, 50443, 50495, 50547, 50599, 50651, 50703, 50755, 50807, 50858, 50910, 50962, 51014, 51066, 51118, 51170, 51222, 51274, 51326, 51377]
ancestor_times = [4.0, 11.130237608828494, 30.970547307245088, 86.17738760127531, 239.79370012757641, 667.2402148800238, 1856.6355334451114, 5166.198660059532, 14375.254655218505, 40000.00000000001]
Ms = [10] #number of importance samples at each locus
Ts = [None] #time cutoffs
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
    coal=coal
  output:
    newick 
  params:
    prefix_in = anc.replace('.anc',''), #prefix of anc and mut files (relate searches for anc/mut files with this prefix)
    prefix_out = newick.replace('.newick','') #prefix of outfile (relate adds its own suffix)
  threads: 1
  resources:
    runtime=60
  group:
    "sample"
  shell:
    '''
    start=$( awk 'NR=={wildcards.locus} {{print $1}}' {input.loci} ) #position of first snp at locus
    stop=$( awk 'NR=={wildcards.locus} {{print $2}}' {input.loci} ) #position of last snp at locus
    module load gcc/11.3.0 #on my server i had to load the same version i built relate with to run it
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
    # snakemake data/SGDP_aDNA_new_chr2_1locus_10M.newick -c1
    #one instance (data/SGDP_aDNA_new_chr2_1locus_10M) took 77 seconds to run

# ------------------- extract subtrees --------------------------------------

rule extract_subtrees:
  input:
    anc="data/SGDP_aDNA_new_chr2.anc",
    mut="data/SGDP_aDNA_new_chr2.mut",
    poplabels="data/SGDP_contemporary_only.poplabels",
    pops="data/contemporary_populations.txt"  # file with one population per line
  output:
    anc_out="data/SGDP_contemporary_only_chr2.anc",
    mut_out="data/SGDP_contemporary_only_chr2.mut",
    poplabel_out="data/SGDP_contemporary_only_chr2.poplabels"
  resources:
    runtime=480  # minutes
  threads: 2

  shell:
    """
    POPS=$(paste -sd, {input.pops})
    
    relate/bin/RelateExtract --mode SubTreesForSubpopulation \
      --anc {input.anc} \
      --mut {input.mut} \
      --poplabels {input.poplabels} \
      --pop_of_interest $POPS \
      -o data/SGDP_contemporary_only_chr2
    """
#snakemake extract_subtrees --profile slurm --jobs 10

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
    runtime=120
  group:
    "extract"
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
    import sys
    from tsconvert import from_newick
    from utils import get_shared_times
    from tqdm import tqdm
    
    #set recursion limit
    sys.setrecursionlimit(3000)
    # open file of trees to read from
    with open(input.newick, mode='r') as f:
      
      # open files to append to
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
            stss.write(",".join([str(i) for i in sts]) + '\n') #append as new line
 
            # get coalescence times 
            cts = sorted([tree.time(i) for i in tree.nodes() if not tree.is_sample(i)]) #coalescence times, in ascending order
            ctss.write(",".join([str(i) for i in cts]) + '\n') #append as new line

# ---------------- process times -----------------------------

# now we process the times, potentially cutting off the tree (to ignore distant past) and getting the exact quantities we need for inference

processed_times = shared_times.replace('.stss','_{T}T.{end}')
ends = ['stss_logdet','stss_inv','btss','lpcs']

rule process_times:
  input:
    stss = shared_times,
    ctss = coal_times,
    coal = coal
  output:
    expand(processed_times, end=ends, allow_missing=True)
  threads: 1 
  resources:
    runtime=15
  group:
    "process"
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
    from utils import chop_shared_times, center_shared_times, log_coal_density
    from tqdm import tqdm

    # determine time cutoff
    T = wildcards.T #get time cutoff
    T = None if T=='None' else float(T) #format correctly

    # effective population size
    epochs = np.genfromtxt(input.coal, skip_header=1, skip_footer=1) #time at which each epoch starts (and the final one ends)
    Nes = 0.5/np.genfromtxt(input.coal, skip_header=2)[2:] #effective population size during each epoch

    # open file of shared times to read from
    with open(input.stss, 'r') as stss:
      with open(input.ctss, 'r') as ctss:
       
        # open files to write to
        with open(output[0], 'a') as stss_logdet:
          with open(output[1], 'a') as stss_inverted:
            with open(output[2], 'a') as btss:
              with open(output[3], 'a') as lpcs:
            
                # loop over trees at this locus 
                for sts,cts in tqdm(zip(stss,ctss), total=int(wildcards.M)):
            
                  # load shared time matrix in vector form
                  sts = np.fromstring(sts, dtype=float, sep=',') #convert from string to numpy array
    
                  # chop
                  sts = chop_shared_times(sts, T=T) #chop shared times to ignore history beyond T
                  
                  # convert to matrix form
                  #k = int((np.sqrt(1+8*(len(sts)-1))+1)/2) #get number of samples (from len(sts) = k(k+1)/2 - k + 1)
                  k = int((np.sqrt(1+8*len(sts))-1)/2) #get size of matrix (from sum_i=0^k i = k(k+1)/2), allows for non-contemporary samples
                  sts_mat = np.zeros((k,k)) #initialize matrix
                  #sts_mat[np.triu_indices(k, k=1)] = sts[1:] #fill in upper triangle
                  #sts_mat = sts_mat + sts_mat.T + np.diag([sts[0]]*k) #add lower triangle and diagonal
                  sts_mat[np.triu_indices(k, k=0)] = sts #convert to numpy matrix
                  sts_mat = sts_mat + sts_mat.T - np.diag(np.diag(sts_mat)) #fill in all entries
                  sts = sts_mat
                  
                  # sample times
                  x = np.diag(sts)
                  x = np.max(x) - x
                  sample_times = np.sort(x)
    
                  # center
                  sts = center_shared_times(sts) 
            
                  # determinant
                  sts_logdet = np.linalg.slogdet(sts)[1] #magnitude of log determinant (ignore sign)
                  stss_logdet.write(str(sts_logdet) + '\n') #append as new line 
            
                  # inverse
                  sts = np.linalg.inv(sts) #inverse
                  sts = sts[np.triu_indices(k-1, k=0)] #
                  stss_inverted.write(",".join([str(i) for i in sts]) + '\n') #append as new line

                  # branching times
                  cts = np.fromstring(cts, dtype=float, sep=',') 
                  Tmax = cts[-1] #time to most recent common ancestor
                  if T is not None and T < Tmax:
                      Tmax = T #farthest time to go back to
                  bts = Tmax - np.flip(cts) #branching times, in ascending order
                  bts = bts[bts>0] #remove branching times at or before T
                  bts = np.append(bts, Tmax) #append total time as last item      
                  btss.write(",".join([str(i) for i in bts]) + '\n') #append as new line
                 
                  # probability of coalescence times under neutral coalescent --
                  lpc = log_coal_density(coal_times=cts, sample_times=sample_times, Nes=Nes, epochs=epochs, T=Tmax) #log probability density of coalescence times
                  lpcs.write(str(lpc) + '\n') #append as new line 

                  #about 3 minutes per job

# ----------- estimate dispersal ------------------------

# and now we bring in our processed times across chromosomes and loci to estimate a dispersal rate

dispersal_rate = processed_times.replace('_chr{CHR}','').replace('_{locus}locus','').replace('{end}','sigma')

rule dispersal_rate:
  input:
    stss_logdet = expand(processed_times, end=['stss_logdet'], CHR=CHRS, locus=dispersal_loci, allow_missing=True),
    stss_inv = expand(processed_times, end=['stss_inv'], CHR=CHRS, locus=dispersal_loci, allow_missing=True),
    btss = expand(processed_times, end=['btss'], CHR=CHRS, locus=dispersal_loci, allow_missing=True),
    lpcs = expand(processed_times, end=['lpcs'], CHR=CHRS, locus=dispersal_loci, allow_missing=True),
    locations = locations,
    sts = shared_times.replace('{CHR}',str(CHRS[0])).replace('{locus}',str(dispersal_loci[0])) #any chr and locus will do, just getting sampling times
  output:
    sigma = dispersal_rate
  threads: 1
  resources:
    runtime=360
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
    from spacetrees import estimate_dispersal

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
      bts = []
      with open(f, 'r') as fi:
        for line in fi:
          bts.append(np.fromstring(line, dtype=float, sep=','))
      btss.append(bts)
    lpcs = [] #log probability of coalescence times
    for f in tqdm(input.lpcs):
      lpcs.append(np.loadtxt(f))
    locations = np.loadtxt(input.locations) #location of each sample
    # sampling times
    sts = np.loadtxt(input.sts, delimiter=',')[0] #a vectorized shared times matrix to get sample times from
    k = int((np.sqrt(1+8*len(sts)-1)+1)/2) #get size of matrix (from sum_i=0^k i = k(k+1)/2)
    mat = np.zeros((k,k))
    mat[np.triu_indices(k, k=0)] = sts #convert to numpy matrix
    mat = mat + mat.T - np.diag(np.diag(mat))      
    x = np.diag(mat) #shared times with self
    x = np.max(x) - x #sampling times
    sample_times = np.sort(x) #sampling times in asceding order


    
# --- Consistency checks ---
#    num_samples = locations.shape[0]
#   assert k == num_samples, f"Mismatch: sampling matrix size ({k}) ≠ number of samples ({num_samples})"
#   
#    for i, block in enumerate(stss_inv):
#        for j, mat in enumerate(block):
#           assert mat.shape == (num_samples, num_samples), \
#                f"stss_inv[{i}][{j}] has shape {mat.shape}, expected ({num_samples},{num_samples})"
#    
#    for i, logdet in enumerate(stss_logdet):
#        assert logdet.shape[0] == len(stss_inv[i]), \
#            f"logdet[{i}] count ({logdet.shape[0]}) ≠ stss_inv[{i}] count ({len(stss_inv[i])})"
#    
#    assert len(btss) == len(stss_inv) == len(lpcs) == len(stss_logdet), \
#        "Mismatch in the number of loci across input lists"

    # estimate dispersal rate
    def callbackF(x):
      '''print updates during numerical search'''
      print('{0: 3.6f}   {1: 3.6f}   {2: 3.6f}   {3: 3.6f}'.format(x[0], x[1], x[2], x[3]))
    sigma = estimate_dispersal(locations=locations, shared_times_inverted=stss_inv, shared_times_logdet=stss_logdet,
                               branching_times=btss, sample_times=sample_times, logpcoals=lpcs,
                               callbackF=callbackF)
    with open(output.sigma, 'w') as f:
      f.write(','.join([str(i) for i in sigma])) #save

# took 5hr 45 min for 100 dispersal loci

# ---------------- convert files to npy to reduce I/O operations -----------------
# processed_prefix_stss = "data/SGDP_aDNA_new_filtered_chr2_{locus}locus_{M}M"
# rule convert_stss_to_npy:
#   input:
#     stss = processed_prefix_stss + ".stss"
#   output:
#     stss_npy = processed_prefix_stss + "_stss.npy"
#   run:
#     import os

#     os.environ["OMP_NUM_THREADS"] = str(threads)
#     os.environ["GOTO_NUM_THREADS"] = str(threads)
#     os.environ["OPENBLAS_NUM_THREADS"] = str(threads)
#     os.environ["MKL_NUM_THREADS"] = str(threads)
#     os.environ["VECLIB_MAXIMUM_THREADS"] = str(threads)
#     os.environ["NUMEXPR_NUM_THREADS"] = str(threads)

#     import numpy as np
#     arr = np.loadtxt(input.stss, delimiter=',')
#     np.save(output.stss_npy, arr)


# # processed_times_npy = processed_times.replace('.{end}','_{end}.npy')
# # ends = ['stss_inv','btss','lpcs']

# processed_prefix = "data/SGDP_aDNA_new_filtered_chr2_{locus}locus_{M}M_{T}T"

# rule convert_to_npy:
#   input:
#     stss_inv = processed_prefix + ".stss_inv",
#     btss     = processed_prefix + ".btss",
#     lpcs     = processed_prefix + ".lpcs"
#   output:
#     stss_inv_npy = processed_prefix + "_stss_inv.npy",
#     btss_npy     = processed_prefix + "_btss.npy",
#     lpcs_npy     = processed_prefix + "_lpcs.npy"
#   run:
#     import os

#     os.environ["OMP_NUM_THREADS"] = str(threads)
#     os.environ["GOTO_NUM_THREADS"] = str(threads)
#     os.environ["OPENBLAS_NUM_THREADS"] = str(threads)
#     os.environ["MKL_NUM_THREADS"] = str(threads)
#     os.environ["VECLIB_MAXIMUM_THREADS"] = str(threads)
#     os.environ["NUMEXPR_NUM_THREADS"] = str(threads)

#     import numpy as np
#     arr = np.loadtxt(input.stss_inv, delimiter=',')
#     np.save(output.stss_inv_npy, arr)
#     arr = np.loadtxt(input.btss, delimiter=',')
#     np.save(output.btss_npy, arr)
#     arr = np.loadtxt(input.lpcs, delimiter=',')
#     np.save(output.lpcs_npy, arr)

# ----------------------- locate ancestors -----------------------

# finally, we use our processed times and dispersal rate to locate the genetic ancestor at a particular locus for a particular sample and time
# TODO: it might be better to locate internal nodes of a tree

ancestor_locations = processed_times.replace('.{end}','_{s}s_{t}t.locs')
# processed_prefix = "data/SGDP_aDNA_new_filtered_chr2_{locus}locus_10M_{T}T"
# processed_prefix_stss = "data/SGDP_aDNA_new_filtered_chr2_{locus}locus_{M}M"

rule locate_ancestors:
  input:
    stss = shared_times,
    stss_inv = processed_times.replace('{end}','stss_inv'),
    btss = processed_times.replace('{end}','btss'),
    lpcs = processed_times.replace('{end}','lpcs'),
    locations = locations,
    sigma = dispersal_rate
  output:
    ancestor_locations
  threads: 1
  resources:
    runtime=360 #Finished 78% in 2 hours
  group:
    "locate"
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

    T = wildcards.T #get time cutoff
    T = None if T=='None' else float(T) #format correctly

    # load input data
    # shared times
    stss = np.loadtxt(input.stss, delimiter=',') #list of vectorized shared times matrices
    k = int((np.sqrt(1+8*len(stss[0])-1)+1)/2) #get size of matrix (from sum_i=0^k i = k(k+1)/2)
    mat = np.zeros((k,k))
    mat[np.triu_indices(k, k=0)] = stss[0] #convert to numpy matrix
    mat = mat + mat.T - np.diag(np.diag(mat))      
    x = np.diag(mat) #shared times with self
    x = np.max(x) - x #sampling times
    sample_times = np.sort(x) #sampling times in asceding order
    stss_mat = [] #list of chopped shared times matrices in matrix form
    for sts in stss:
      sts = chop_shared_times(sts, T=T) #chop shared times to ignore history beyond T
      mat = np.zeros((k,k))
      #mat[np.triu_indices(k, k=1)] = sts[1:] #convert to numpy matrix
      #mat = mat + mat.T + np.diag([sts[0]]*k)      
      mat[np.triu_indices(k, k=0)] = sts #convert to numpy matrix
      mat = mat + mat.T - np.diag(np.diag(mat))      
      stss_mat.append(mat)
    stss = stss_mat
    # shared times chopped centered inverted
    stss_inv = np.loadtxt(input.stss_inv, delimiter=',') #list of vectorized chopped centered inverted shared times matrices
    k = k-1 #get size of matrix
    stss_inv_mat = [] #list of chopped shared times matrices in matrix form
    for sts_inv in stss_inv:
      mat = np.zeros((k,k))
      mat[np.triu_indices(k, k=0)] = sts_inv #convert to numpy matrix
      mat = mat + mat.T - np.diag(np.diag(mat))      
      stss_inv_mat.append(mat)
    stss_inv = stss_inv_mat
    # branching times
    btss = []
    with open(input.btss, 'r') as f:
      for line in f:
        bts = np.fromstring(line, dtype=float, sep=',') #coalescence times in ascending order
        btss.append(bts)
    # coal probs
    lpcs = np.loadtxt(input.lpcs) #log probability of coalescence times
    #locations 
    locations = np.loadtxt(input.locations) #location of each sample
    # dispersal and branching rate
    sigma = np.loadtxt(input.sigma, delimiter=',') #mle dispersal rate and branching rate
    phi = sigma[-1] #branching rate
    sigma = _sds_rho_to_sigma(sigma[:-1]) #dispersal as covariance matrix

    # calculate importance weights
    lbds = np.array([_log_birth_density(bts, sample_times, phi) for bts in btss]) #log probability densities of birth times
    log_weights = lbds - lpcs #log importance weights

    # locate ancestors
    s = wildcards.s
    if s == 'All': #an option to locate the ancestors of all samples
      samples = range(k+1)   
    else:
      samples = [int(s)]
    t = wildcards.t
    if t == 'All': #an option to locate at pretermined list of times 
      times = ancestor_times
    else: 
      times = [float(t)]
    ancestor_locations = locate_ancestors(samples=samples, times=times, 
                                          shared_times_chopped=stss, shared_times_chopped_centered_inverted=stss_inv, locations=locations, 
                                          sigma=sigma, log_weights=log_weights)
    with open(output[0], 'a') as f:
      for anc_loc in ancestor_locations:
        f.write(','.join([str(int(anc_loc[0]))] + [str(i) for i in anc_loc[1:]]) + '\n') #save

        # Timing out at 2 hours. Need to increase

# ---------------- dummy rule to run everything you need -----------------

rule all:
  input:
    expand(ancestor_locations, CHR=[CHRS], locus=ancestor_loci, M=Ms, T=Ts, s=['All'], t=['All']) 
