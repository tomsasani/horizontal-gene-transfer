from __future__ import print_function
import pysam
import argparse
from collections import defaultdict, namedtuple, Counter
import glob

p = argparse.ArgumentParser()
p.add_argument("--bam")
p.add_argument("-mq_thresh", default=0, type=int)
p.add_argument("-outname", default='out')
args = p.parse_args()

samfile = pysam.AlignmentFile(args.bam, 'rb')

nfrags = defaultdict(list)

fh = open(args.outname + '.denominator.pickle', 'w')
import pickle

for i, read in enumerate(samfile.fetch()):
    if i % 1000000 == 0 and i != 0: print ("done with {} million reads".format(i / 1000000))
    nfrags[read.query_name].append(read.mapping_quality)

pickle.dump(nfrags, fh)
