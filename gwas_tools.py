'''This module contains functions for parsing GWAS data text files'''

import os
import re
import subprocess

import pandas as pd
import numpy as np

def get_snp_coor(ext, snp_file, bed_dir):
    df_dict = {'JH': '', 'CD': '', 'ST': ''}
    for key in df_dict:
        df = pd.read_table(key+ext, header=1)
        df_dict[key] = df
    p = re.compile('_(?P<initials>[CJS][DHT])_*_d(?P<day>\d{1,2})_')
    bed_files = [x for x in os.listdir(bed_dir) if x.endswith('.bed')]
    bed_files.sort(key=lambda x: int(p.search(x).group('day')))
    for x in bed_files:
        m = p.search(x)
        bed_file = os.path.join(bed_dir, x)
        snp_file = '{}_impSNPs.bed'.format(m.group('initials')) # remove this line if running on sentinel SNPs
        bedtools = subprocess.run(
            'module load bedtools\nbedtools intersect -a {} -b {}'.format(
                snp_file,
                bed_file),
            shell=True, stdout=subprocess.PIPE)
        out = bedtools.stdout.decode('utf-8')
        #data = np.array([x.split('\t') for x in out.strip().split('\n')])
        data = np.array([x.split('\t') for x in out.split('\n')][:-1])
        #print(data)
        if len(data)>0:
            atac = np.where([any(x==data[:,3]) for x in df_dict[m.group('initials')]['SNP name']],
                            'Yes', 'No')
            #print('high')
        else:
            atac = ['No']*df_dict[m.group('initials')].shape[0]
            #print('low')
        new_col = 'ATAC_d{}'.format(m.group('day'))
        df_dict[m.group('initials')][new_col] = atac
        
    return df_dict
    
def write_df(df, out_ext):
    for x in df:
        df[x].to_csv(x+out_ext, index=False, sep='\t')
    
    return

