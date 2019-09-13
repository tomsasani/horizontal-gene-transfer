import cPickle
import pickle
from quicksect import IntervalTree
from collections import defaultdict
import argparse
import sys

p = argparse.ArgumentParser()
p.add_argument('--pickle')
p.add_argument('--tes')
args = p.parse_args()


def read_exclude(path, chrom=None):
    if path is None:
        return None
    tree = defaultdict(IntervalTree)
    # in case they sent a region.
    if chrom is not None:
        chrom = chrom.split(":")[0]

    added = 0
    for i, line in enumerate((gzip.open if path.endswith(".gz") else open)(path)):
        toks = line.rstrip().split("\t")
        # skip header if necessary.
        if i == 0:
            try:
                int(toks[1])
            except ValueError:
                continue
        if chrom is not None and toks[0] != chrom:
            continue
        added += 1
        # NOTE: made the 'other' value equal to whatever annotation I put as the
        # fourth field in the exclude BED file, such as SINE/LINE identity
        try: 
            tree[toks[0]].add(int(toks[1]), int(toks[2]), other=toks[3])
        except:
            tree[toks[0]].add(int(toks[1]), int(toks[2]))

    if added == 0:
        sys.stderr.write("didnt add any intervals to exclude for chrom:%s\n" % chrom)
    return tree

te_tree = read_exclude(args.tes)

with open(args.pickle, 'r') as fh:
    read_dict = cPickle.load(fh)
    #read_dict = pickle.load(fh)
    n_hq = 0
    for r in read_dict:
        segments = [s for s in read_dict[r] if 'Un' not in s[0] and 'M' not in s[0]]
        if len(segments) == 1: continue
        hq_te = False
        hq_vv = False
        for seg in segments:
            chrom, mq, start = str(seg[0]), int(seg[-1]), int(seg[1])
            if mq < 20: continue
            if 'K3L' not in chrom and len(te_tree[chrom].search(start, start)) > 0:
                hq_te = True
            if 'K3L' in chrom:
                hq_vv = True
        if hq_te and hq_vv: 
            for seg in segments:
                info = [r]
                info.extend(map(str, seg))
                print '\t'.join(info)
            n_hq += 1

    print >> sys.stderr, "total of {} fragments supporting R/V chimeras".format(n_hq)

