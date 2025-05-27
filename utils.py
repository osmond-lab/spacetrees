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

  TMRCA = tree.time(tree.root) #tmrca of tree
  k = len(samples) #number of samples

  sts = []
  for i in range(k):
    for j in range(i,k): #just building upper triangular part of symmetric matrix, including diagonal in case of non-contemporary samples
      st = TMRCA - tree.tmrca(samples[i],samples[j]) #shared time of pair, ordered to align with locations
      sts.append(st)

  return sts

def chop_shared_times(shared_times, T=None):

  TMRCA = np.max(shared_times) #tmrca

  if T is None or T > TMRCA: #dont cut if dont ask or cut time older than MRCA
    pass
  else:
    shared_times = shared_times - (TMRCA - T) #calculate shared times since T

  return shared_times

def center_shared_times(shared_times):
 
  n = len(shared_times) #number of samples in subtree
  Tmat = np.identity(n) - [[1/n for _ in range(n)] for _ in range(n)]; Tmat = Tmat[0:-1]; #matrix for mean centering
  stc = np.matmul(Tmat, np.matmul(shared_times, np.transpose(Tmat))) #center shared times in subtree
 
  return stc

def log_coal_density(times, Nes, epochs=None, T=None):

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
    if T is not None:
        times = times[times < T] #ignore old times
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

    # now add the probability of lineages not coalescing by T 
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

