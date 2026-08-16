import numpy as np

def loci_positions(mut, outfile):

  with open(outfile, 'w') as fout:
    with open(mut, "r") as fin:
      ix = 1
      next(fin) #skip header
      for i,line in enumerate(fin):
        if i==0:
          start = int(line.split(';')[1]) #position of first snp
        elif int(line.split(';')[4]) == ix: #when we reach next locus 
          end = pos #position of last snp at previous locus
          fout.write(str(start) + ' ' + str(end) + '\n') #position of first and last snp at previous locus
          start = int(line.split(';')[1]) #position of first snp at this locus
          ix = ix + 1 #next locust
        pos = int(line.split(';')[1]) #position of this snp
    fout.write(str(start) + ' ' + str(pos) + '\n') #position of first and last snp at last locus

def get_shared_times(tree, samples):

  """
  Shared times between every pair of samples (including each sample with itself),
  packed as the upper triangle (with diagonal) of the symmetric matrix, row-major.
  A sample's diagonal entry is TMRCA - its own age, so non-contemporary (ancient)
  samples get their own diagonal value instead of one shared TMRCA for everyone.
  """

  TMRCA = tree.time(tree.root) #tmrca of tree
  k = len(samples) #number of samples

  sts = []
  for i in range(k):
    for j in range(i,k): #upper triangular part of symmetric matrix, including diagonal (samples may not be contemporary)
      st = TMRCA - tree.tmrca(samples[i],samples[j]) #shared time of pair, ordered to align with locations
      sts.append(st)

  return sts

def split_shared_times(shared_times, T=None):

  """
  Chop shared times at cutoff T and split samples into isolated subtrees:
  groups of samples whose lineages still share time with each other since T
  (history beyond T disconnects lineages that coalesced before the cutoff,
  so they no longer belong to the same subtree). Returns a list of
  (sample_ids, shared_times_submatrix) pairs, one per subtree.

  shared_times is the upper triangle (with diagonal) of the full pairwise
  matrix, as returned by get_shared_times -- diagonal entries can differ
  per sample (TMRCA - that sample's own age), so this doesn't assume
  contemporary sampling.
  """

  k = int((np.sqrt(1 + 8 * len(shared_times)) - 1) / 2) #number of samples (len(shared_times) = k(k+1)/2)
  mat = np.zeros((k, k))
  mat[np.triu_indices(k, k=0)] = shared_times
  mat = mat + mat.T - np.diag(np.diag(mat)) #full symmetric shared-times matrix

  TMRCA = np.max(mat) #true root time, assuming at least one contemporary sample (its diagonal entry is TMRCA - 0)
  if T is None or T > TMRCA: #dont cut if dont ask or cut time older than MRCA
    return [(np.arange(k), mat)]

  chopped = T - (TMRCA - mat) #shared times since T (negative between samples in different subtrees)

  samples = np.arange(k)
  subtrees = []
  # a sample older than T (negative chopped diagonal) hadn't been sampled yet within the cutoff
  # window, so it can't share time with anyone -- exclude it up front rather than let it try (and
  # fail, since chopped[i] >= 0 would be false even for i itself) to join a subtree
  taken = np.diag(chopped) < 0
  while not taken.all(): #while some samples not yet assigned to a subtree
    i = np.argmax(~taken) #next sample not yet assigned
    withi = chopped[i] >= 0 #samples that still share time with i since T
    taken |= withi
    subtrees.append((samples[withi], chopped[np.ix_(withi, withi)]))

  return subtrees

def center_shared_times(shared_times):
 
  n = len(shared_times) #number of samples in subtree
  Tmat = np.identity(n) - [[1/n for _ in range(n)] for _ in range(n)]; Tmat = Tmat[0:-1]; #matrix for mean centering
  stc = np.matmul(Tmat, np.matmul(shared_times, np.transpose(Tmat))) #center shared times in subtree
 
  return stc

def log_coal_density(coal_times, sample_times, Nes, epochs=None, T=None):

    """
    log probability of coalescent times under standard neutral/panmictic coalescent,
    with samples entering the tree at sample_times (0 for contemporary samples) instead
    of assuming everyone is sampled at time 0.
    """

    if epochs is None and len(Nes) == 1:
        epochs = [0, max(coal_times)] #one big epoch
        Nes = [Nes[0], Nes[0]] #repeat the effective population size so same length as epochs

    sample_times = np.sort(sample_times) #callers don't need to remember to pre-sort

    logp = 0 #initialize log probability
    prevLambda = 0 #initialize coalescent intensity
    if T is not None:
        coal_times = coal_times[coal_times < T] #ignore old coalescence times
        sample_times = sample_times[sample_times < T] #ignore old sampling times
    k = int(np.sum(sample_times == 0)) #number of contemporary samples
    i = k #index of next sample time to enter
    myIntensityMemos = _coal_intensity_memos(epochs, Nes) #intensities up to end of each epoch

    # probability of each coalescence time
    for ct in coal_times: #for each coalescence time
        while i < len(sample_times) and sample_times[i] < ct: #if next sampling happens before next coalescence
            kchoose2 = k * (k - 1) / 2 #binomial coefficient
            Lambda = _coal_intensity_using_memos(sample_times[i], epochs, myIntensityMemos, Nes) #coalescent intensity up to sampling time
            logpk = - kchoose2 * (Lambda - prevLambda) #log probability of no coalescence
            logp += logpk #add log probability
            k += 1 #add the new sample lineage
            i += 1 #move to next sample time
            prevLambda = Lambda #update intensity
        kchoose2 = k * (k - 1) / 2 #binomial coefficient
        Lambda = _coal_intensity_using_memos(ct, epochs, myIntensityMemos, Nes) #coalescent intensity up to coalescence time
        ie = np.digitize(np.array([ct]), epochs) #epoch at the time of coalescence
        logpk = np.log(kchoose2 * 1 / (2 * Nes[ie])) - kchoose2 * (Lambda - prevLambda) #log probability (waiting times are time-inhomogeneous exponentially distributed)
        logp += logpk #add log probability
        prevLambda = Lambda #update intensity
        k -= 1 #update number of lineages

    # deal with any remaining sampling events after the last coalescence
    while i < len(sample_times):
        kchoose2 = k * (k - 1) / 2 #binomial coefficient
        Lambda = _coal_intensity_using_memos(sample_times[i], epochs, myIntensityMemos, Nes) #coalescent intensity up to sampling time
        logpk = - kchoose2 * (Lambda - prevLambda) #log probability of no coalescence
        logp += logpk #add log probability
        k += 1 #add the new sample lineage
        i += 1 #move to next sample time
        prevLambda = Lambda #update intensity

    # now add the probability of any remaining lineages not coalescing by T
    if k > 1 and T is not None: #if we have more than one lineage remaining
        kchoose2 = k * (k - 1) / 2 #binomial coefficient
        Lambda = _coal_intensity_using_memos(T, epochs, myIntensityMemos, Nes) #coalescent intensity up to time T
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

