from __future__ import print_function
import pysam
import argparse
from collections import defaultdict, namedtuple, Counter
import glob

def is_hq_mapped(ct, min_match=0.1):
    """
    >>> is_hq_mapped([(0, 30), (1, 5), (2, 2), (1, 3), (4, 111)])
    True
    >>> is_hq_mapped([(0, 10), (1, 5),  (5, 135)])
    False
    >>> is_hq_mapped([(5, 10), (0, 140)])
    True
    """
    MATCH = 0
    total = sum([c[1] for c in ct])
    
    match_frac = sum([c[1] if c[0] == MATCH else 0 for c in ct]) / float(total)
    
    return match_frac > min_match

def seq_entropy(seq):
    """
    >>> seq_entropy('AAAAAAAAAT')
    0.9
    >>> seq_entropy('CTCCAGCCTC')
    0.6
    """
    base_count = Counter(seq)
    most_common = base_count.most_common()
    total_bases = sum([x[-1] for x in most_common])
    most_freq_frac = most_common[0][-1] / float(total_bases)

    return most_freq_frac

def main():
    
    p = argparse.ArgumentParser()
    p.add_argument("--bam")
    p.add_argument("-outname", default="read_out", help="filename for output BAM")
    p.add_argument("-mq_thresh", default=0, type=int)
    args = p.parse_args()

    samfile = pysam.AlignmentFile(args.bam, 'rb')
    strand_dict = {True:'-', False:'+'}

    for i, read in enumerate(samfile.fetch('chrK3L30', 38200, 38350)):
        if i % 1000000 == 0 and i != 0: print ('done with {} million reads'.format(i / 1000000))

        # first, we'll apply a sane filter on mapping quality.
        # only alignments with MQ >= threshold are included
        if read.mapping_quality < args.mq_thresh: continue

        # for split reads, we can get all of the info
        # we need from the SA tag of the primary alignment.
        # so, we can skip all of the secondary alignments in the BAM
        if read.is_secondary or read.is_proper_pair: continue
        
        # grab the SA tag if one is available (i.e., it's a split read).
        # if it's not a split read, make sure it's not a proper pair.
        # "improper" pairs will include fragments in which the first alignment is to rabbit
        # and the second is to VV (or vice versa), and we want those
        sa_tag = None
        try: 
            sa_tag = read.get_tag('SA')
        except KeyError: pass
       
        contig = read.reference_name

        rb_chimeric = None
        vv_chimeric = None

        print (read.query_name)

        # if the primary alignment has an SA tag, we can quickly determine
        # whether we're interested in the read (i.e., if it has alignments to
        # both rabbit and VV).
        if sa_tag is not None:
            chrs_in_tag = (tag.split(',')[0] for tag in sa_tag.split(';'))
            chrs_in_tag = [c for c in chrs_in_tag if c != '']

            rb_chimeric = ('K3L' not in contig and any(['K3L' in c for c in chrs_in_tag])) 
            vv_chimeric = ('K3L' in contig and any(['K3L' not in c for c in chrs_in_tag])) 

        # if there's no SA tag (i.e., we're looking at only an "improper pair"), check that
        # the mate is aligned to a different organism's contig
        else:
            rb_chimeric = 'K3L' not in contig and 'K3L' in read.next_reference_name
            vv_chimeric = 'K3L' in contig and 'K3L' not in read.next_reference_name

        if not (rb_chimeric or vv_chimeric): continue

        # add info about this alignment to a dictionary keyed on the query (fragment) name
        start, end = read.reference_start, read.reference_end
        seg = (contig, start, strand_dict[read.is_reverse], read.cigarstring, int(read.mapping_quality))
        read_dict[read.query_name].append(seg)

        # add the same information to the dictionary, without having to access the alignments 
        # explicitly, just by using the SA tag
        if sa_tag is not None:
            for s in sa_tag.split(';'):
                if s == '': continue
                read_dict[read.query_name].append(tuple(s.split(',')[:-1]))

    # pickle and dump the output so we can use it quickly later
    fh = open(args.outname + '.pickle', 'w')
    import pickle
    pickle.dump(read_dict, fh)

    print ("dumped to pickle")

if __name__ == "__main__":
    import doctest
    read_dict = defaultdict(list)
    doctest.testmod()
    main()
