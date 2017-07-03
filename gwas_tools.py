'''This module contains functions for parsing GWAS data text files'''

import os
import re
import subprocess

import pandas as pd
import numpy as np

jh = pd.read_table('JH_vDH-SNPinfo_c1.txt', header=1)
cd = pd.read_table('CD_vDH-SNPinfo_c1.txt', header=1)
st = pd.read_table('ST_vDH-SNPinfo_c1.txt', header=1)
df_dict = {'JH': jh, 'CD': cd, 'ST': st}

def get_snp_coor(snp_file, bed_dir):
    bed_files = [x for x in os.listdir(bed_dir) if x.endswith('.bed')]
    p = re.compile('_(?P<initials>[CJS][DHT])_*_(?P<day>d\d{1,2})_')
    for x in bed_files:
        m = p.search(x)
        bed_file = os.path.join(bed_dir, x)
        bedtools = subprocess.run(
            'module load bedtools\nbedtools intersect -a {} -b {}'.format(
                snp_file,
                bed_file),
            shell=True, stdout=subprocess.PIPE)
        out = bedtools.stdout.decode('utf-8')
        data = np.array([x.split('\t') for x in out.strip().split('\n')])
        if len(data)>0:
            atac = np.where([any(x==data[:,3]) for x in jh['SNP name']],
                            'Yes', 'No')
        else:
            atac = ['No']*df_dict[m.group('initials')].shape[0]
        new_col = 'ATAC_{}'.format(m.group('day'))
        df_dict[m.group('initials')][new_col] = atac
        
    return df_dict
    
def write_df(df):
    for x in df:
        file_name = '{}_ATAC_intersect.txt'.format(x)
        df[x].to_csv(file_name, index=False, sep='\t')
    
    return

