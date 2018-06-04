import pysam
import random
import sys

bam = pysam.AlignmentFile(sys.argv[1], 'rb')
bam2 = pysam.AlignmentFile(sys.argv[2], 'rb')
outbam = pysam.AlignmentFile("last_with_Hqual_scores.sorted.bam", "wb", template=bam2)

all_query_qualities = []
seq_lengths = []
done = False
q = None
for line in bam:
    all_query_qualities.extend(line.query_qualities)
    if len(all_query_qualities) > 1000000:
        break

for line in bam2:
    qual = []
    for position in range(len(line.seq)):
        qual.append(random.choice(all_query_qualities))
    line.query_qualities = qual
    outbam.write(line)

bam.close()
bam2.close()
outbam.close()

