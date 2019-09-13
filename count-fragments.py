import pickle
import argparse

p = argparse.ArgumentParser()
p.add_argument('--pickle')
p.add_argument('-mq_thresh', type=int, default=0)
args = p.parse_args()

t = 0

with open(args.pickle, 'r') as fh:
    nfrags = pickle.load(fh)
    for i,frag in enumerate(nfrags):
        if i % 100000 == 0: print ('done with {} fragments'.format(i))
        if all([q >= args.mq_thresh for q in nfrags[frag]]): t += 1

print ('total of {} fragments with all alignments having MQ >= {}'.format(t, args.mq_thresh))

