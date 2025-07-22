import numpy as np

# Original total number of loci
total_loci = 51846
max_allowed_locus = 51398  # Filter threshold

# Dispersal loci: 100 evenly spaced (excluding index 0)
num_dispersal = 100
dispersal_loci = list(np.linspace(0, total_loci - 1, num=num_dispersal, dtype=int))
dispersal_loci = [locus for locus in dispersal_loci if locus != 0 and locus < max_allowed_locus]

# Ancestor loci: 1000 evenly spaced (excluding index 0)
num_ancestor = 1000
ancestor_loci = list(np.linspace(0, total_loci - 1, num=num_ancestor, dtype=int))
ancestor_loci = [locus for locus in ancestor_loci if locus != 0 and locus < max_allowed_locus]

# Ancestor times: 10 log-spaced values between 4 and 40000
ancestor_times = list(np.logspace(np.log10(4), np.log10(40000), 10))

print(ancestor_times)