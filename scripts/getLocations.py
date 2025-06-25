import numpy as np

lonlats = []
with open("/scratch/m/mmosmond/raghavs/spacetrees/data/SGDP_metadata.279public.21signedLetter.samples.txt", "rb") as f:
    next(f)
    for line in f:
        l = line.strip().split(b'\t')
        if l != [b'']:
            lonlat = []
            lonlat.append(l[4])
            lonlat.append(l[11])
            lonlat.append(l[12])
            lonlats.append(lonlat)

# get accession names in metadata
accessions = [i[0] for i in lonlats]

IDS = []
with open("/scratch/m/mmosmond/raghavs/spacetrees/data/SGDP_sub_aDNA.poplabels", "rb") as f:
    next(f)
    for line in f:
        IDS.append(line.strip().split(b' ')[0])

order = []
for ID in IDS:
    ix = accessions.index(ID)
    order.append(ix)

locations = np.array([list(map(float, [i[2], i[1]])) for i in lonlats])[order]
np.save("/scratch/m/mmosmond/raghavs/spacetrees/data/projected_lpcation.npy", locations)

