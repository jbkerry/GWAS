#!/usr/bin/env python

from __future__ import print_function, division

import argparse
import os

def check_file(x):
    if not os.path.exists(x):
        print('Warning, could not locate {}. Skipping.'.format(x))
        return False
    else:
        return True

def store_mutants(dirs):
    """Stored positional mutation counts from VEP output
    
    Parameters
    ----------
    dirs : array-like
        directories containing output files from crispr_tools.py
    
    Returns
    -------
    del_dict : dict
        contains deleterious mutation counts per position
    nondel_dict : dict
        contains non-deleterious mutantion counts per position
        
    """
    
    del_dict = {}
    nondel_dict = {}
    for curr_dir in dirs:
        file_list = (x for x in os.listdir(curr_dir))
        for x in file_list:
            file_name = os.path.join(curr_dir, x)
            with open(file_name) as f:
                next(f)
                for line in f:
                    cols = line.strip().split('\t')
                    chr_num, coor, mut = cols[0].split('_')
                    if (cols[5] == 'DELETERIOUS') | (cols[5] == 'INDEL'):
                        del_dict[coor] = del_dict.get(coor, 0) + 1
                    elif cols[5] == 'NON-DELETERIOUS':
                        nondel_dict[coor] = nondel_dict.get(coor, 0) +1
                        
    return del_dict, nondel_dict

def write_mutants(dels, nondels, gene, sample, outfile='positional_mutants.txt', wigs=True):
    """Outputs number of deleterious and non-deleterious mutations per 
    position as a text file, and optionally as separate wig files
    
    Parameters
    ----------
    dels : dict
        contains deleterious mutation counts per position
    nondels : dict
        contains non-deleterious mutation counts per position
    outfile : str
        name of output text file, default='positional_mutants.txt'
    wigs : bool
        also output count data as wig tracks, default=True
        
    """
    
    outfile = '_'.join([sample, outfile])
    out = open(outfile, 'w')
    out.write('Position\tDeleterious_mutations\tNon-deleterious_mutations\n')
    if wigs:
        print('Generating wig tracks')
        header = 'track type=wiggle_0 name="{} wig" description="{} mutants"' \
                 ' visibility=full color={} offset=0\n'
        chr_name = 'chr17'
        if gene == 'ASXL1': chr_name = 'chr20'
        chr_line = 'variableStep chrom={}\n'.format(chr_name)
        delwig = open('deleterious_{}.wig'.format(sample), 'w')
        delwig.write(header.format('del', 'deleterious', '177,0,0'))
        delwig.write(chr_line)
        nondelwig = open('non-deleterious_{}.wig'.format(sample), 'w')
        nondelwig.write(header.format('nondel', 'non-deleterious', '0,0,177'))
        nondelwig.write(chr_line)
    for pos, count in sorted(dels.items()):
        out.write('{}\t{}\t{}\n'.format(pos, count, nondels.get(pos, 0)))
        if wigs:
            delwig.write('{} {}\n'.format(pos, count))
    for pos, count in sorted(nondels.items()):
        if pos not in dels:
            out.write('{}\t{}\t{}\n'.format(pos, 0, count))
        if wigs:
            nondelwig.write('{} {}\n'.format(pos, count))
    
    out.close()
    if wigs:
        delwig.close()
        nondelwig.close()
    
    return None

def generate_track(gene, sample):
    w2bw_cmd = 'module load ucsctools\n' \
               'wigToBigWig {0}_{1}.wig /databank/igenomes/Homo_sapiens/UCSC/' \
               'hg19/Sequence/WholeGenomeFasta/chr_sizes.txt ' \
               '/public/{2}/{0}_{1}.bw'
    print('Converting wigs to bigWigs')
    os.system(w2bw_cmd.format('non-deleterious', sample, os.environ['LOGNAME']))
    os.system(w2bw_cmd.format('deleterious', sample, os.environ['LOGNAME']))
    print('Generating Track Hub for UCSC')
    hub_path = '/public/{}/crispr_hub'.format(os.environ['LOGNAME'])
    if not os.path.exists(hub_path): os.mkdir(hub_path)
    if not os.path.exists(hub_path): os.mkdir(hub_path)
    
    hubtxt_path = os.path.join(hub_path, '{}_hub.txt'.format(gene))
    genometxt_path = os.path.join(hub_path, '{}_genomes.txt'.format(gene))
    tracktxt_path = os.path.join(hub_path, '{}_tracks.txt'.format(gene))
    
    if not os.path.exists(hubtxt_path):
        with open(hubtxt_path, 'w') as f:
            f.write('hub {0}\nshortLabel {0}\nlongLabel {0} CRISPR mutations\n'
                    'genomesFile {0}_genomes.txt\nemail jon.kerry@imm.ox.ac.uk\n'.format(gene))
    
    if not os.path.exists(genometxt_path):
        with open(genometxt_path, 'w') as f:
            f.write('genome hg19\ntrackDb {}_tracks.txt'.format(gene))
    
    trackStr = """
#--------------------------------------

track {1}_{0}
container multiWig
shortLabel {1}, {0}
longLabel {1}, {0}, deleterious = orange, non-deleterious = blue
type bigWig
visibility full
aggregate transparentOverlay
showSubtrackColorOnUi on
windowingFunction maximum
configurable on
autoScale on
alwaysZero on
dragAndDrop subtracks

track {0}_del
parent {1}_{0}
bigDataUrl http://userweb.molbiol.ox.ac.uk/public/{2}/deleterious_{0}.bw
shortLabel {0} del
longLabel {1} {0} B deleterious
type bigWig
color 204,77,0
priority 100

track {0}_nondel
parent {1}_{0}
bigDataUrl http://userweb.molbiol.ox.ac.uk/public/{2}/non-deleterious_{0}.bw
shortLabel {0} nondel
longLabel {1} {0} non-deleterious
type bigWig
color 0,204,255
priority 110

#--------------------------------------

    
    """.format(sample, gene, os.environ['LOGNAME'])
    
    with open(tracktxt_path, 'a') as f:
        f.write(trackStr)
    
    hub_address = ''.join(['http://sara.molbiol.ox.ac.uk', hubtxt_path])
    print("""
Add this URL to the My Data -> Track Hubs -> My Hubs tab in any UCSC hg19 genome browser session: {}
You only need to do this once per gene. If you have already added the URL for {} and this was just a different sample, the session will update automatically upon refereshing.
    """.format(hub_address, gene))
    
if __name__ == '__main__':
    
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-d',
        '--dir',
        type = str,
        default = '.',
        help = 'Colon-separated list of directories e.g. dir1:dir2:dir3',
        required = True,
    )
    parser.add_argument(
        '-s',
        '--sample',
        type = str,
        help = 'Sample name, e.g. 12d_B_allOligos',
        required = True,
    )
    parser.add_argument(
        '-g',
        '--gene',
        type = str,
        help = 'Name of gene e.g. TP53',
        required = True,
    )
    
    args = parser.parse_args()
    
    top_list = args.dir.split(':')
    dir_list = []
    for oligo in top_list:
        del_dir = os.path.join(oligo, 'DELETERIOUS')
        nondel_dir = os.path.join(oligo, 'NON-DELETERIOUS')
        if check_file(del_dir): dir_list.append(del_dir)
        if check_file(nondel_dir): dir_list.append(nondel_dir)
            
    del_dict, nondel_dict = store_mutants(dir_list)
    write_mutants(del_dict, nondel_dict, args.gene, args.sample)
    generate_track(args.gene, args.sample)
    
    
    
    