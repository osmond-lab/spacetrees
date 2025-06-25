import numpy as np

total_loci = 51846
num_to_use = 100
dispersal_loci = list(np.linspace(0, total_loci - 1, num=num_to_use, dtype=int))
ancestor_loci = [int(i * 51845 / 999) for i in range(1000)]
ancestor_times = np.logspace(np.log10(4),np.log10(40000),10)
ancestor_loci.pop(0) #Don't want loci 0
dispersal_loci.pop(0)
#print(dispersal_loci[0]," " ,ancestor_loci[0])
