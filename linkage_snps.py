from __future__ import print_function, division

import re
import subprocess
import math
import sys

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import argparse
import pysam
from scipy import stats

def group_imputations():
    """Generates dictionary with every proximal SNP as a key referencing
    a dictionary with every linked imputed SNP as key with that SNP's
    location and base change as a value
    
    Returns
    -------
    imp_dict : dict
        Contains imputed SNPs linked to proximal SNPs
    
    """
    
    imp_dict = {}
    with open('/t1-data1/WTSA_Dev/jkerry/BloodATAC/Ron_imputed_SNPs.txt') as f:
        for x in f:
            (chr_name, start, stop, imp_snp, ref_al, alt_al,
             pr_snp, desc, snp_type) = x.strip().split('\t')
            if pr_snp not in imp_dict:
                imp_dict[pr_snp] = {}
            loc = ':'.join((chr_name, stop))
            ref_change = '/'.join((ref_al, alt_al))
            if (len(ref_al)==1) & (len(alt_al)==1):  # only include imputed SNP if it really is a SNP and not an indel
                imp_dict[pr_snp].update({imp_snp: '_'.join((loc, ref_change))})
    
    return imp_dict
    

def match_vdh_to_ngs(initials, length):
    """Determines if any of the van De Harst 75 SNPs are in the
    specified NGS sequencing file
    
    Parameters
    ----------
    initials : str
        Initials of the individual who has been sequenced
    length : int
        Read length of the sequences
        
    """
    
    cnames=('SNP', 'GWAS', 'A', 'C', 'T', 'G')
    df = pd.DataFrame(columns=cnames)
    
    bam = '/t1-data1/WTSA_Dev/jkerry/BloodATAC/ATAC_K4_{}.sorted.bam'.format(
        initials
    )
    
    snp_dict = {}
    info_dict = {}
    
    with open('/t1-data1/WTSA_Dev/jkerry/BloodATAC/' \
              'vanDeHarst2012_SNPs_hg19.txt') as f:
        for x in f:
            chr_name, snp, loc, base, genes = x.strip().split('\t')
            snp_dict[snp] = {}
            coor = ':'.join((chr_name, loc))
            info_dict[snp] = coor
            ea, oa = base.split('/')
            chr_num = int(chr_name[3:])
            loc = int(loc)
            start = loc - length
            stop = loc
            
            sf = pysam.AlignmentFile(bam, 'rb')
            for r in sf.fetch(chr_name, start, stop):
                pos = loc-(r.reference_start+1)
                seq = r.query_sequence
                if r.cigartuples[0][0]==0:
                    m_len = r.cigartuples[0][1]-1
                    if (pos>=0) & (pos<=m_len):
                        seq_base = seq[pos]
                        snp_dict[snp][seq_base] = snp_dict[snp].get(
                            seq_base, 0) + 1
                        
            if snp_dict[snp]:
                GWASCount = 0
                nGWASCount = 0
                base_dict = {'A': 0, 'C': 0, 'T': 0, 'G': 0}
                for snp_base in snp_dict[snp]:
                    if snp_base in base_dict:
                        base_dict[snp_base]+=snp_dict[snp][snp_base]
                df = df.append({'SNP': snp,
                                'GWAS': ea,
                                'A': base_dict['A'],
                                'C': base_dict['C'],
                                'T': base_dict['T'],
                                'G': base_dict['G']}, ignore_index=True)
    
    return df


## at this point the script has read to each read and matched it to a SNP location where possible. For a given SNP it has calculated how many times the nucleotide matches that in the vanDeHarst paper or one of the other 3 nucleotides.
## the dataframe also contains the 'replicate' i.e. the different day of the experiment (not really a replicate so this should be fixed later when taking averages) and the initials of the person

#GroupList = ['SNP']
#grouped_df = df.groupby(GroupList)
#Counter=0
#SNPList = []
#df_2 = pd.DataFrame(columns=["SNP","GWAS","A","C","T","G"])
#for key, item in grouped_df:
#    SNPList.append(grouped_df.get_group(key).reset_index()['GWAS'][0])
#    ASum = sum(grouped_df.get_group(key)['A'])
#    CSum = sum(grouped_df.get_group(key)['C'])
#    TSum = sum(grouped_df.get_group(key)['T'])
#    GSum = sum(grouped_df.get_group(key)['G'])
#    InsertList = [key,grouped_df.get_group(key).reset_index()['GWAS'][0],ASum,CSum,TSum,GSum]
#    df_2.loc[Counter] = InsertList
#    Counter+=1
#

def get_genotype(initials, cutoff, length, data=df):

    imp_snp_out = open('{}_Hap_c{}.bed'.format(initials, cutoff), 'w')
    pr_snp_out = open('{}_vDH_SNPs_c{}.bed'.format(initials, cutoff), 'w')
    ColumnList = ['A','C','T','G']
    ATots = []
    CTots = []
    TTots = []
    GTots = []
    UserDict = {}
    code_dict = {}
    for i, row in data.iterrows():
        
        binary = np.where(row[['A', 'C', 'T', 'G']]>1, 1, 0)  #if read cut-off a nucleotide must appear in at least 2 sequences to be counted thereby reducing the chance of a sequencing error being included
        total = sum(binary)
        percentage = [x/total*100 for x in binary]
        ATots.append(Percentage[0])
        CTots.append(Percentage[1])
        TTots.append(Percentage[2])
        GTots.append(Percentage[3])
        #print(Percentage)
        #print(BinaryList)
        Status = "N"
        imp_snp_num = len(imp_dict[row['SNP']])
        imp_yes = 0
        imp_no = 0
        #if df_2.iloc[i][df_2.iloc[i]['GWAS']]>ReadCutoff: ##if read cut-off a nucleotide must appear in at least 2 sequences to be counted thereby reducing the chance of a sequencing error being included
            #Status = "pos"
        pr_snp_chr, pr_snp_loc = info_dict[row['SNP']].split(':')
        pr_snp_out.write("{}({})\n".format('\t'.join((pr_snp_chr,
                                                      int(pr_snp_loc)-1,
                                                      pr_snp_loc,
                                                      row['SNP'])),
                                           row['GWAS']))
        
        ##Check imputed SNPs here
        pr_snp_change = imp_dict[row['SNP']][row['SNP']].split('_')
        pr_ref,pr_alt = pr_snp_change.split('/')
        code = 0
        #print("For proxy SNP "+df_2.iloc[i]['SNP']+":")
        #print("\tReference = {0}, status = ".format(PrImpRefAl),end="")
        if row[pr_ref]>cutoff:
            #print("present")
            code += 1
        #else:
            #print("absent")
        #print("\tImputedSNP = {0}, status = ".format(PrImpAltAl),end="")
        if row[pr_alt]>cutoff:
            #print("present")
            code+=10
        #else:
            #print("absent")
        #print("\tEA = {0}, status = ".format(df_2.iloc[i]['GWAS']),end="")
        #if df_2.iloc[i][df_2.iloc[i]['GWAS']]>ReadCutoff:
            #print("present")
        #else:
            #print("absent")
        code_dict[row['SNP']] = code
        if row[row['GWAS']]>cutoff: # Will only check imputed SNPs if vDH EA is present. Alter this IF statement to include other groups like the imputed SNP change, if different   
            status = "Y"
            for imp_snp in imp_dict[row['SNP']]:
                imp_loc, imp_change = imp_dict[row['SNP']][imp_snp].split('_')
                imp_chr, imp_coor = imp_loc.split(':')
                imp_ref, imp_alt = imp_change.split('/')
                alt_lengths = [len(x) for x in imp_alt.split(',')]
                #print("SVs = {0}, seq lengths = {1}".format(Alt_Al,newSeqLengths))
                loc = int(imp_coor)
                start = loc - length
                stop = loc
                #BaseDict = {'A': 0, 'C': 0, 'T': 0, 'G': 0}
                base_dict = {}
                
                sf = pysam.AlignmentFile(bam, 'rb')
                for r in sf.fetch(imp_chr, start, stop):
                    pos = loc-(r.reference_start+1)
                    seq = r.query_sequence
                    if r.cigartuples[0][0]==0:
                        m_len = r.cigartuples[0][1]-1
                        if (pos>=0) & (pos<=m_len):
                            for x in alt_lengths:
                                end_pos = pos + x
                                seq_base = seq[pos:end_pos]
                                base_dict[seq_base] = base_dict.get(
                                    base_dict[seq_base], 0) + 1
                all_alts = imp_alt.split(',')
                running_yes = 0
                running_no = 0
                for alt_al in all_alts:
                    #print(This_Alt_Al)
                    if alt_al not in base_dict:
                    #print("\tFor imputed SNP "+impSNP+", this code isn't smart enough yet to check if this alternative SNP allele exists")
                        running_no += 1
                    else:
                        if base_dict[alt_al] > cutoff:
                            #print("\tFor imputed SNP "+impSNP+", imputed SNP exists!")
                            running_yes += 1
                        else:
                            #print("\tFor imputed SNP "+impSNP+", imputed SNP does not exist")
                            running_no += 1
                            
                if running_yes>=1:
                    #print("\tFor imputed SNP "+impSNP+", imputed SNP exists!")
                    value = 1000
                    imp_yes += 1
                else:
                    #print("\tFor imputed SNP "+impSNP+", imputed SNP does not exist")
                    value = 250
                    imp_no += 1
                imp_pr_snp = '{}({})'.format(imp_snp, row['SNP'])
                imp_snp_out.write('{}\n'.format('\t'.join(imp_chr,
                                                          loc-1,
                                                          loc,
                                                          imp_pr_snp,
                                                          value)))
            collate = '_'.join(imp_snp_num, imp_yes, imp_no)
            
            ### to here
            
            LDSNP_num[df_2.iloc[i]['SNP']] = Collate
            #print("\tOut of a total of {0} imputed SNP, {1} were present and {2} were not".format(ImpSNPNum,ImpYes,ImpNo))
                    
                #if bool(SNPdict[SNP])==False:
                    #x=1
                #else:
                    #GWASCount = 0
                    #nGWASCount = 0
                    #BaseDict = {'A': 0, 'C': 0, 'T': 0, 'G': 0}
                    #for SNPbase in SNPdict[SNP].keys():
                        #BaseDict[SNPbase]+=SNPdict[SNP][SNPbase]
                    #InsertList = [SNP,EA,BaseDict['A'],BaseDict['C'],BaseDict['T'],BaseDict['G']]
                    #df.loc[RowCounter] = InsertList
                    #RowCounter+=1
            
            #print(ImpDict[df_2.iloc[i]['SNP']])
            
        Geno = "Hom"
        if Total>1:
            Geno = "Het"
        elif Total==0:
            Geno = "Ambiguous"
        String = Geno+","+Status
        if Geno=="Ambiguous":
            String = Geno
        UserDict[df_2.iloc[i]['SNP']] = String
    impSNPout.close()
    prSNPout.close()
#output = open(InputInitial+"_vDH-SNPinfo_c"+str(ReadCutoff)+".txt","w")
#output.write(InputInitial+", read cutoff>"+str(ReadCutoff)+"\n")
#output.write("SNP name\tHet/Hom\tref/alt\tvDH EA\tLD SNPs\n")
#
#
#for ThisSNP in sorted(InfoDict.keys()):
#    output.write(ThisSNP+"\t")
#    if ThisSNP not in UserDict.keys():
#        output.write("not detected\t-\t-\t-\n")
#    else:
#        if UserDict[ThisSNP]!="Ambiguous":
#            Geno,Status = UserDict[ThisSNP].split(',')
#            Comp = "n/a"
#            if Geno=="Hom":
#                if CodeDict[ThisSNP]==1:
#                    Comp = "ref"
#                elif CodeDict[ThisSNP]==10:
#                    Comp = "alt"
#                else:
#                    Comp = "error"
#            LD = "n/a"
#            if Status=="Y":
#                TLD,YLD,NLD = LDSNP_num[ThisSNP].split('_')
#                LD = YLD+"/"+TLD
#        else:
#            Geno = "Ambiguous"
#            Comp = "-"
#            Status = "-"
#            LD = "-"
#        
#        output.write("{0}\t{1}\t{2}\t{3}\n".format(Geno,Comp,Status,LD))
#output.close()
#
##print(df_2[:20])
#
#N = len(ATots)
#ind = np.arange(N)
#width = 0.7
##F1911E
#AColours = np.where([x=='A' for x in SNPList],'#B1D9F4','white')
##AColours = np.where([x=='A' for x in SNPList],0.9,0.6)
#CColours = np.where([x=='C' for x in SNPList],'#B1D9F4','white')
#TColours = np.where([x=='T' for x in SNPList],'#B1D9F4','white')
#GColours = np.where([x=='G' for x in SNPList],'#B1D9F4','white')
#
#pA = plt.bar(ind,ATots,width,color=AColours,hatch='/')
#pC = plt.bar(ind,CTots,width,bottom=ATots,color=CColours,hatch='.')
#pT = plt.bar(ind,TTots,width,bottom=[x+y for x,y in zip(ATots,CTots)],color=TColours,hatch='\\')
#pG = plt.bar(ind,GTots,width,bottom=[x+y+z for x,y,z in zip(ATots,CTots,TTots)],color=GColours,hatch='X')
#
##pA = plt.bar(ind,ATots,width,color='blue',alpha=AColours)
##pC = plt.bar(ind,CTots,width,bottom=ATots,color='red',alpha=0.5)
##pT = plt.bar(ind,TTots,width,bottom=[x+y for x,y in zip(ATots,CTots)],color='green',alpha=0.5)
##pG = plt.bar(ind,GTots,width,bottom=[x+y+z for x,y,z in zip(ATots,CTots,TTots)],color='yellow',alpha=0.5)
#
#plt.ylabel('Percentage of locus (%)',fontsize=18)
#plt.xlabel('Van de Harst SNPs',fontsize=18)
#plt.title(InputInitial+", Nucleotide base at Van de Harst SNP locations",fontsize=20,fontweight='bold',y=1.03)
#plt.xticks([i+(width/2) for i in ind], df_2['SNP'],rotation='vertical')
#plt.yticks(np.arange(0,101,50))
#plt.ylim(0,105)
#plt.margins(0.01)
#plt.tick_params(axis='both',direction='out',length=5,top='off',right='off')
#
#snp_patch = mpatches.Patch(color='#B1D9F4',label='GWAS\nSNP base')
#second_legend = plt.legend(handles=[snp_patch],loc=2,bbox_to_anchor=(1.01,0.8),fontsize=18)
#plt.gca().add_artist(second_legend)
#ax = plt.gca()
#leg = ax.get_legend()
#leg.legendHandles[0].set_edgecolor('black')
#
#plt.legend((pA[0],pC[0],pT[0],pG[0]),('A','C','T','G'),loc=2,bbox_to_anchor=(1.01,1),fontsize=18)
#ax = plt.gca()
#leg = ax.get_legend()
#i=0
#while i<=3:
#    leg.legendHandles[i].set_color('white')
#    leg.legendHandles[i].set_edgecolor('black')
#    i+=1
#plt.subplots_adjust(bottom=0.15,right=0.87)
##plt.savefig("ST_SNP.png")
##plt.show()

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-i',
        '--initials',
        type = str,
        help = 'Initials of person.',
        required = True,
    )
    parser.add_argument(
        '-r',
        '--read_length',
        type = int,
        help = 'Read length (bp)',
        required = True,
    )
    parser.add_argument(
        '-c',
        '--cutoff',
        type = int,
        help = 'Sequencing depth minimum threshold',
        required = True,
    )
    
    args = parser.parse_args()
