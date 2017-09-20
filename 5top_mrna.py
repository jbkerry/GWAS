#!/usr/bin/env python

from __future__ import print_function

import re
import argparse

from Bio.Seq import Seq

def get_top_genes(f_name):
    p = re.compile('\s+')
    p_5TOP = re.compile('^[TC]{4,14}')
    #p_5TOP_like = re.compile('^[AG]{1,3}')
    p_5TOP_like = re.compile('^[AGTC]{0,3}[AG]')
    p_5TOP_2 = re.compile('^[TC]{5,14}')
    top_list = []
    top_like_list = []
    with open(f_name) as f:
        for x in f:
            name, seq = p.split(x.strip())
            seq = seq.upper()
            if name.endswith('-'):
                min_seq = Seq(seq)
                seq = str(min_seq.complement())
            m = p_5TOP_like.match(seq)
            if p_5TOP.match(seq):
                if name[:-2] not in top_list:
                    top_list.append(name[:-2])
                    print('{}\t5\'TOP mRNA\t{}'.format(name[:-2], seq))
            elif m:
                pur_end = m.end()
                if p_5TOP_2.match(seq[pur_end:]):
                    if name[:-2] not in top_like_list:
                        top_like_list.append(name[:-2])
                        print('{}\t5\'TOP-like mRNA\t{}'.format(name[:-2], seq))
                        
if __name__ == '__main__':
    
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-f',
        '--file',
        type = str,
        help = 'Path to tab file.',
        required = True,
    )
    args = parser.parse_args()
    get_top_genes(args.file)
    