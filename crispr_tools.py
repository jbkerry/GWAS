#!/usr/bin/env python

from __future__ import print_function, division

import argparse
import math
import os
import re
import shutil
import subprocess

silents = {'1A': ('7578495-C-G', '7578498-C-T', '7578501-C-G', '7578504-A-T',
                  '7578507-G-T', '7578510-G-C', '7578513-C-T', '7578516-G-A',
                  '7578519-C-A', '7578534-C-T'),
           '1C': ('7578495-C-G', '7578498-C-T', '7578501-C-G', '7578504-A-T',
                  '7578510-G-C', '7578513-C-T', '7578516-G-A', '7578519-C-A',
                  '7578534-C-T'),
           '1G': ('7578495-C-G', '7578498-C-T', '7578501-C-G', '7578504-A-T',
                  '7578507-G-C', '7578510-G-C', '7578513-C-T', '7578516-G-A',
                  '7578519-C-A', '7578534-C-T'),
           '1T': ('7578495-C-G', '7578498-C-T', '7578501-C-G', '7578504-A-T',
                  '7578507-G-A', '7578510-G-C', '7578513-C-T', '7578516-G-A',
                  '7578519-C-A', '7578534-C-T'),
           '2': ('7578498-C-T', '7578504-A-T', '7578510-G-C', '7578516-G-A',
                 '7578519-C-A', '7578523-TGG-TGGG', '7578524-GG-GGG',
                 '7578534-C-T'),
           '3': ('7578495-C-G', '7578501-C-G', '7578507-G-A', '7578513-C-T',
                 '7578519-C-A', '7578534-C-T'),
           'A-wt': ('31022633-T-C', '31022636-C-T', '31022642-C-T',
                    '31022643-G-T', '31022657-C-A', '31022660-A-G'),
           'A-mut': ('31022633-T-C', '31022636-C-T', '31022642-C-T',
                     '31022657-C-A', '31022660-A-G')}
vep_cmd = 'module load ensembl-api/20170130\nperl /package/ensembl/20170130/' \
          'ensembl-tools/scripts/variant_effect_predictor/variant_effect_' \
          'predictor.pl --database --assembly GRCh37 --port 3337 --polyphen b ' \
          '--sift b --force_overwrite'

def collapse_vcfs(oligo, v_dir='.'):
    vcf_list = [x for x in os.listdir(v_dir) if (x.endswith('.vcf')) & (x!='test.vcf')]
    vcf_dict = {}
    bar_dict = {}
    bar_collapse = [0, 0]
    vcf_collapse = [0, 0]
    
    if not os.path.exists('./BARCODES'): os.mkdir('./BARCODES')
    
    for vcf_file in vcf_list:
        snp_list = []
        with open(os.path.join(v_dir, vcf_file)) as f:
            for x in f:
                i = x.strip().split('\t')
                snp_list.append('-'.join([i[1], i[3], i[4]]))
        key = '_'.join(snp_list)
        if all([x in silents[oligo] for x in snp_list]):
            if key in bar_dict:
                bar_dict[key][0]+=1
                bar_collapse[1]+=1
                #os.remove(vcf_file)
            else:
                bar_dict[key] = [1, vcf_file]
                bar_collapse[0]+=1; bar_collapse[1]+=1
                shutil.copyfile(os.path.join(v_dir, vcf_file),
                                os.path.join('./BARCODES', vcf_file))
        else:
            if key in vcf_dict:
                vcf_dict[key][0]+=1
                vcf_collapse[1]+=1
                #os.remove(vcf_file)
            else:
                vcf_dict[key] = [1, vcf_file]
                vcf_collapse[0]+=1; vcf_collapse[1]+=1
    for bar_key, bar_item in bar_dict.items():
        os.rename('./BARCODES/{}'.format(bar_item[1]),
                  './BARCODES/{}_{}.vcf'.format(bar_item[1][:-4], bar_item[0]))
    
    print('{0} barcode-only vcfs were collapsed down to {1} unique vcfs. These'
          ' {1} were copied to the ./BARCODES/ directory'.format(
            bar_collapse[1],
            bar_collapse[0]))
    print('{0} new-mutant vcfs were collapsed down to {1} unique vcfs. These '
          '{1} will now be run through VEP.'.format(
            vcf_collapse[1],
            vcf_collapse[0]))
    return vcf_dict

def run_vep(vcf_dict, oligo, v_dir='.', verbose=False, keep_vep=False):
    p = re.compile('(?P<file_num>\d*\d)')
    p2 = re.compile('SIFT=(?P<sift>.*)\((?P<sift_score>.*)\);PolyPhen=(?P<pp>.*)\((?P<pp_score>.*)\)')
    
    total = len(vcf_dict)
    counter = 0
    del_count = nondel_count = del_total_count = nondel_total_count = 0
    print('Starting VEP runs...')
    
    if not os.path.exists('./DELETERIOUS'): os.mkdir('./DELETERIOUS')
    if not os.path.exists('./NON-DELETERIOUS'): os.mkdir('./NON-DELETERIOUS')
    
    for i, j in vcf_dict.items():
        m = p.search(j[1])
        outfile = '{}_out.txt'.format(m.group('file_num'))
        cmd = '{} -i {} -o {}'.format(vep_cmd,
                                      os.path.join(v_dir, j[1]),
                                      outfile)
        if verbose: print('Running VEP on {}...'.format(j[1]))
        #subprocess.call(cmd, shell=True)
        vep = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
        vep_stdout = vep.stdout.read()
        #if verbose: print(vep_stdout.strip())
        if verbose: print('\t...complete')
        
        outstring = ''
        deleterious = False
        ambiguous = False
        score_dict = {}
        with open(outfile) as f:
            for x in f:
                if x.startswith('#'): continue
                cols = x.strip().split('\t')
                extras = cols[-1]
                m2 = p2.search(extras)
                if m2:
                    ss = float(m2.group('sift_score'))
                    pps = float(m2.group('pp_score'))
                    s_term = m2.group('sift')
                    pp_term = m2.group('pp')
                    if (ss<0.05) | (pps>=0.15):
                        effect = 'DELETERIOUS'
                        deleterious = True
                    else:
                        effect = 'NON-DELETERIOUS'
                else:
                    ss = s_term = '-'
                    pps = pp_term = '-'
                    effect = 'NON-DELETERIOUS'
                chr_num, base, var = cols[0].split('_')
                try:
                    ref, alt = var.split('/')
                except ValueError:
                    print('Ambigous mutation detected in {}. Skipping.'.format(
                        j[1]))
                    ambiguous = True
                    break
                list_id = '-'.join([base, ref, alt])
                if (ref == '-') | (alt == '-') | (len(ref) != len(alt)):
                    if list_id not in silents[oligo]:
                        effect = 'INDEL'
                        deleterious = True
                if list_id in silents[oligo]: effect = 'BARCODE'
                if cols[0] not in score_dict:
                    score_dict[cols[0]] = (cols[0], ss, s_term, pps, pp_term, effect)
                elif (effect == 'DELETERIOUS') & (score_dict[cols[0]][5] != 'DELETERIOUS'):
                    score_dict[cols[0]] = (cols[0], ss, s_term, pps, pp_term, effect)
                elif (ss != '-') & (pps != '-'):
                    if (score_dict[cols[0]][1]=='-') & (score_dict[cols[0]][3]=='-'):
                        score_dict[cols[0]] = (cols[0], ss, s_term, pps, pp_term, effect)
                    elif (ss < score_dict[cols[0]][1]) & (pps > score_dict[cols[0]][3]):
                        score_dict[cols[0]] = (cols[0], ss, s_term, pps, pp_term, effect)
        if not ambiguous:
            for key, item in score_dict.items():
                outstring = '\n'.join([outstring, '\t'.join(map(str, item))])
            out_file = 'mut_id{}_{}.txt'.format(m.group('file_num'), j[0])
            if deleterious:
                out_dir = './DELETERIOUS'
                del_count += 1
                del_total_count += j[0]
            else:
                out_dir = './NON-DELETERIOUS'
                nondel_count += 1
                nondel_total_count += j[0]
                
            with open(os.path.join(out_dir, out_file), 'w') as f_out:
                f_out.write('Variant\tSIFT_score\tSIFT_term\tPolyPhen_score\tPolyPhen_term\tErica_term\n')
                f_out.write('{}\n'.format(outstring.lstrip()))
            if verbose:
                print('Collapsed idential variants and wrote results to {} in '
                      '{}\n'.format(out_file, out_dir))
        if not keep_vep:    
            os.remove(outfile)
            os.remove('{}_summary.html'.format(outfile))
        
        counter+=1
        for y in range(1, 21, 1):    
            if counter == math.ceil((total/20)*y):
                #print('total = {}, counter = {}; y = {}'.format(total, counter, y))
                print('** {}% complete **'.format(y*5))
                if verbose: print()
                break
    
    if not keep_vep: shutil.rmtree('_Inline')
    print('{1} VEP output files had only non-deleterious mutations and were '
          'stored in the ./NON-DELETERIOUS/ directory ({1} unique events, {2} '
          'total events)\n{3} VEP output files had at least one deleterious '
          'mutation and were stored in the ./DELETERIOUS/ directory ({3} '
          'unique events, {4} total events).'.format(nondel_count,
                nondel_total_count, del_count, del_total_count))
    return None

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-d',
        '--dir',
        type = str,
        default = '.',
        help = 'Directory containing vcf files. Default = current directory',
        required = False,
    )
    parser.add_argument(
        '-o',
        '--oligo',
        type = str,
        choices = ['1A', '1C', '1G', '1T', '2', '3', 'A-wt', 'A-mut'],
        help = 'The number of the oligo',
        required = True,
    )
    parser.add_argument(
        '-v',
        '--verbose',
        action = 'store_true',
        help = 'Print progress output for each individual VEP run',
        required = False,
    )
    parser.add_argument(
        '-k',
        '--keep_vep',
        action = 'store_true',
        help = 'Keep the Ensembl summary html and txt files from VEP runs. ' \
               'Deleted by default.',
        required = False,
    )
    
    args = parser.parse_args()
    
    vcf_dict = collapse_vcfs(str(args.oligo), v_dir=args.dir)
    run_vep(vcf_dict,
            str(args.oligo),
            v_dir=args.dir,
            verbose=args.verbose,
            keep_vep=args.keep_vep)