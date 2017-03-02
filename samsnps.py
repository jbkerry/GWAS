#!/usr/bin/env python

from __future__ import print_function
from scipy import stats
import re,subprocess,pandas as pd, numpy as np, matplotlib.pyplot as plt, math, matplotlib.patches as mpatches

ReadLength = 40
cnames=['SNP','GWAS','Non-GWAS','Rep','Person']
df = pd.DataFrame(columns=cnames)
RowCounter=0
RepDict = {7:1,10:2,13:3,16:4}
#BAMfile = "/t1-data/user/milne_group/ChIPseq/SEM/MLL_N/MLL_sorted.bam"
BAMfiles = ["/t1-data1/WTSA_Dev/jkerry/BloodATAC/ATAC_CD34_d7_CD_rep1.sorted.bam","/t1-data1/WTSA_Dev/jkerry/BloodATAC/ATAC_CD34_d10_CD_rep1.sorted.bam","/t1-data1/WTSA_Dev/jkerry/BloodATAC/ATAC_CD34_d13_CD_rep1.sorted.bam","/t1-data1/WTSA_Dev/jkerry/BloodATAC/ATAC_CD34_d16_CD_rep1.sorted.bam","/t1-data1/WTSA_Dev/jkerry/BloodATAC/ATAC_CD34_d7_ST_rep1.sorted.bam","/t1-data1/WTSA_Dev/jkerry/BloodATAC/ATAC_CD34_d10_ST_rep1.sorted.bam","/t1-data1/WTSA_Dev/jkerry/BloodATAC/ATAC_CD34_d13_ST_rep1.sorted.bam","/t1-data1/WTSA_Dev/jkerry/BloodATAC/ATAC_CD34_d16_ST_rep1.sorted.bam","/t1-data1/WTSA_Dev/jkerry/BloodATAC/ATAC_CD34_d7_JH_rep1.sorted.bam","/t1-data1/WTSA_Dev/jkerry/BloodATAC/ATAC_CD34_d10_JH_rep1.sorted.bam","/t1-data1/WTSA_Dev/jkerry/BloodATAC/ATAC_CD34_d13_JH_rep1.sorted.bam","/t1-data1/WTSA_Dev/jkerry/BloodATAC/ATAC_CD34_d16_JH_rep1.sorted.bam"]
#BAMfiles = ["/t1-data1/WTSA_Dev/jkerry/BloodATAC/ATAC_CD34_d7_CD_rep1.sorted.bam","/t1-data1/WTSA_Dev/jkerry/BloodATAC/ATAC_CD34_d10_CD_rep1.sorted.bam","/t1-data1/WTSA_Dev/jkerry/BloodATAC/ATAC_CD34_d13_CD_rep1.sorted.bam","/t1-data1/WTSA_Dev/jkerry/BloodATAC/ATAC_CD34_d16_CD_rep1.sorted.bam"]
for BAMfile in BAMfiles:
    Rep = re.findall('d[0-9]*',BAMfile[-30:])[0]
    DayNum = int(Rep[1:])
    Replicate = "Rep"+str(RepDict[DayNum])
    PersonIni = str(BAMfile)[-18:-16]
    SNPdict = {}
    Lines = [x.rstrip('\n') for x in open('/t1-data1/WTSA_Dev/jkerry/HiC/tests/vanDeHarst2012_SNPs.txt')]
    
    for i in Lines:
        Chr, SNP, Loc, Base, Genes = i.split('\t')
        SNPdict[SNP] = {}
        EA, OA = Base.split('/')
        chrNum = int(re.findall('[0-9]*', Chr)[0])
        Loc = int(re.sub(',', '', Loc[1:-1]))
        start=Loc-ReadLength
        end=Loc
        GetReads = subprocess.Popen("samtools view "+BAMfile+" chr"+str(chrNum)+":"+str(start)+"-"+str(end), shell=True, stdout=subprocess.PIPE)
        SAMlines = GetReads.stdout.read().rstrip('\n').split('\n')
        for j in SAMlines:
            if j!="":
                parts = j.split('\t')
                ReadLoc = parts[3]
                CIGAR = parts[5]
                Align = re.findall('[0-9]*[A-Z]', CIGAR)
                Pos = Loc-int(ReadLoc)
                Seq = parts[9]
                if Align[0][-1:]=='M':
                    CheckLen = int(Align[0][:-1])-1
                    if (Pos>=0) & (Pos<=CheckLen):
                        SeqBase = Seq[Pos]
                        if SeqBase not in SNPdict[SNP].keys():
                            SNPdict[SNP].update({SeqBase: 1})
                        else:
                            SNPdict[SNP][SeqBase]+=1
        if bool(SNPdict[SNP])==False:
            x =1
        else:
            GWASCount = 0
            nGWASCount = 0
            for SNPbase in SNPdict[SNP].keys():
                if SNPbase==EA:
                    GWASCount+=SNPdict[SNP][SNPbase]
                else:
                    nGWASCount+=SNPdict[SNP][SNPbase]
            InsertList = [SNP,GWASCount,nGWASCount,Replicate,PersonIni]
            df.loc[RowCounter] = InsertList
            RowCounter+=1

## at this point the script has read to each read and matched it to a SNP location where possible. For a given SNP it has calculated how many times the nucleotide matches that in the vanDeHarst paper or one of the other 3 nucleotides.
## the dataframe also contains the 'replicate' i.e. the different day of the experiment (not really a replicate so this should be fixed later when taking averages) and the initials of the person

GroupList = ['Person','SNP']
grouped_df = df.groupby(GroupList)
width = 0.5
namelist = ['GWAS', 'Non-GWAS']
Counter=0
df_2 = pd.DataFrame(columns=["Person","SNP","LogP","GWAS"])
for key, item in grouped_df:
    GWASMean = grouped_df.get_group(key)['GWAS'].mean()
    GWASstdev = grouped_df.get_group(key)['GWAS'].std()
    nGWASMean = grouped_df.get_group(key)['Non-GWAS'].mean()
    nGWASstdev = grouped_df.get_group(key)['Non-GWAS'].std()
    mannWhit = stats.mannwhitneyu(grouped_df.get_group(key)['GWAS'],grouped_df.get_group(key)['Non-GWAS'])
    LogP = math.log(float(mannWhit[1]),2)*-1
    data = [GWASMean,nGWASMean]
    error = [GWASstdev,nGWASstdev]
    xlocations = np.array(range(len(data)))+0.5
    valuesList = [grouped_df.get_group(key)['GWAS'].values,grouped_df.get_group(key)['Non-GWAS'].values]
    if GWASMean>=nGWASMean:
        GWASvalue=1
    else:
        GWASvalue=0
    InsertList = [key[0],key[1],LogP,GWASvalue]
    df_2.loc[Counter] = InsertList
    Counter+=1

## at this point the df groups the replicates of the SNPs together (in a person-dependent manner), takes the mean and SD of the GWAS-matching hits for a given SNP, across the 4 replicates, and compares this to mean of the non-GWAS matching hit across the 4 replicates.
## The data points are compared against each other using the Mann-Whitney U test and the direction (i.e. GWAS SNP or not) is noted by checking which mean is higher and putting a '1' in the GWAS column of the df if it is the GWAS nucleotide, or '0' if it isn't
    
grouped_df2 = df_2.groupby('Person')
###fig, axs = plt.subplots(nrows=2,ncols=1)
###Counter=0
###CounterList=[(0,0),(0,1)]
###for key, item in grouped_df2:
###    ax = axs[Counter]
###    initials = key
###    colours = np.where(grouped_df2.get_group(key)['GWAS']==1,'r','b')
###    sizes = np.where(grouped_df2.get_group(key)['GWAS']==1,50,25)
###    xtickRange = range(1,len(grouped_df2.get_group(key)['SNP'])+1)
###    ax.scatter(xtickRange,grouped_df2.get_group(key)['LogP'],c=colours,s=sizes,marker='x')
###    ax.set_xticks(xtickRange,minor=False)
###    ax.set_xticklabels(grouped_df2.get_group(key)['SNP'],rotation='vertical')
###    ax.set_title(initials)
###    ax.set_ylim(0,8)
###    ax.set_xlim(0,len(grouped_df2.get_group(key)['SNP'])+1)
###    ax.set_xlabel('SNP')
###    ax.set_ylabel('-log2(p-value)')
###    Counter+=1
###plt.tight_layout() 
###plt.subplots_adjust(bottom=0.15)   
###plt.show()

## below, the dataframe is pivoted to consists of two major columns: the log p-value and whether the nucleotide is GWAS or not. These two columns are subdivided by person, with the SNP name being the key for each row so that each SNP only has one record.
## where a particular SNP location is not present in a person's sequence file, the NaN for log p-value is filled in with a false '0'

bySNP = pd.pivot_table(df_2, index='SNP', columns='Person', values=['LogP','GWAS'])
bySNP = bySNP.fillna(0)
print(bySNP[:20])
fig, ax = plt.subplots()
xtickRange = range(1,len(bySNP.index.values)+1)
ax.set_xticks(xtickRange,minor=False)
ax.set_xticklabels(bySNP.index.values,rotation='vertical')
ax.set_xlim(0,len(bySNP.index.values)+1)
ax.set_yticklabels(range(0,9),fontsize=14)
ax.set_ylim(0,8)
ax.set_xlabel('SNP',fontsize=20)
ax.set_ylabel('-log2(p-value)',fontsize=20)
markerList = ['o','x','^']
counter = 0
#ax.scatter(xtickRange,bySNP['LogP']['CD'])
for initials in bySNP['LogP'].columns.values:
    colours = np.where(bySNP['GWAS'][initials]==1,'r','b')
    ax.scatter(xtickRange,bySNP['LogP'][initials],c=colours,s=35,label=initials,marker=markerList[counter])
    counter+=1
    #CDcolours = np.where(bySNP['GWAS']['CD']==1,'r','b')
    #STcolours = np.where(bySNP['GWAS']['ST']==1,'r','b')
    #ax.scatter(xtickRange,bySNP['LogP']['CD'],c=CDcolours,label='CD')
    #ax.scatter(xtickRange,bySNP['LogP']['ST'],c=STcolours,marker='x',label='ST')
red_patch = mpatches.Patch(color='red',label='GWAS SNP base')
blue_patch = mpatches.Patch(color='blue',label='Other base')
first_legend = plt.legend(handles=[red_patch,blue_patch])
plt.gca().add_artist(first_legend)
plt.legend(bbox_to_anchor=(0.8,1))
plt.subplots_adjust(bottom=0.15)
#plt.savefig('scatter.png')
#plt.show()

#print(bySNP['LogP'].columns.values)
#print(bySNP.index.values)
#print(bySNP['LogP']['CD'])
#print(bySNP['LogP'])