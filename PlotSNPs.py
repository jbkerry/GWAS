#!/usr/bin/env python

from __future__ import print_function,division
from scipy import stats
import re,subprocess,pandas as pd, numpy as np, matplotlib.pyplot as plt, math, matplotlib.patches as mpatches, sys, getopt

def usage():
    print("usage: samsnps_2.py -i <initials of person> -r <read length (bp)>")

try:
    opts, args = getopt.getopt(sys.argv[1:], 'i:r:h',)
except getopt.GetoptError:
    usage()
    sys.exit(2)

InputInitial = ""
ReadLength = 0
    
if not opts:
    usage()
    sys.exit(2)
else:
    for opt, arg in opts:
        if opt=='-h':
            usage()
            sys.exit(2)
        elif opt=='-i':
            InputInitial = arg
        elif opt=='-r':
            ReadLength = int(arg)
        else:
            usage()
            sys.exit(2)
            
cnames=['SNP','GWAS','A','C','T','G','Rep','Person']
df = pd.DataFrame(columns=cnames)
RowCounter=0
RepDict = {7:1,10:2,13:3,16:4}
BAMfiles = ["/t1-data1/WTSA_Dev/jkerry/BloodATAC/ATAC_CD34_d7_"+InputInitial+"_rep1.sorted.bam","/t1-data1/WTSA_Dev/jkerry/BloodATAC/ATAC_CD34_d10_"+InputInitial+"_rep1.sorted.bam","/t1-data1/WTSA_Dev/jkerry/BloodATAC/ATAC_CD34_d13_"+InputInitial+"_rep1.sorted.bam","/t1-data1/WTSA_Dev/jkerry/BloodATAC/ATAC_CD34_d16_"+InputInitial+"_rep1.sorted.bam"]
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
            BaseDict = {'A': 0, 'C': 0, 'T': 0, 'G': 0}
            for SNPbase in SNPdict[SNP].keys():
                BaseDict[SNPbase]+=SNPdict[SNP][SNPbase]
            InsertList = [SNP,EA,BaseDict['A'],BaseDict['C'],BaseDict['T'],BaseDict['G'],Replicate,PersonIni]
            df.loc[RowCounter] = InsertList
            RowCounter+=1
            


## at this point the script has read to each read and matched it to a SNP location where possible. For a given SNP it has calculated how many times the nucleotide matches that in the vanDeHarst paper or one of the other 3 nucleotides.
## the dataframe also contains the 'replicate' i.e. the different day of the experiment (not really a replicate so this should be fixed later when taking averages) and the initials of the person

GroupList = ['Person','SNP']
grouped_df = df.groupby(GroupList)
Counter=0
SNPList = []
df_2 = pd.DataFrame(columns=["Person","SNP","GWAS","A","C","T","G"])
for key, item in grouped_df:
    SNPList.append(grouped_df.get_group(key).reset_index()['GWAS'][0])
    ASum = sum(grouped_df.get_group(key)['A'])
    CSum = sum(grouped_df.get_group(key)['C'])
    TSum = sum(grouped_df.get_group(key)['T'])
    GSum = sum(grouped_df.get_group(key)['G'])
    InsertList = [key[0],key[1],grouped_df.get_group(key).reset_index()['GWAS'][0],ASum,CSum,TSum,GSum]
    df_2.loc[Counter] = InsertList
    Counter+=1

ColumnList = ['A','C','T','G']
ATots = []
CTots = []
TTots = []
GTots = []

for i in df_2.index.values:
    
    BinaryList = np.where(df_2.iloc[i][ColumnList]>0,1,0)
    Total = sum(BinaryList)
    Percentage = [x/Total*100 for x in BinaryList]
    ATots.append(Percentage[0])
    CTots.append(Percentage[1])
    TTots.append(Percentage[2])
    GTots.append(Percentage[3])
    #print(Percentage)

N = len(ATots)
ind = np.arange(N)
width = 0.7
AColours = np.where([x=='A' for x in SNPList],'#F1911E','white')
CColours = np.where([x=='C' for x in SNPList],'#F1911E','white')
TColours = np.where([x=='T' for x in SNPList],'#F1911E','white')
GColours = np.where([x=='G' for x in SNPList],'#F1911E','white')

pA = plt.bar(ind,ATots,width,color=AColours,hatch='/')
pC = plt.bar(ind,CTots,width,bottom=ATots,color=CColours,hatch='.')
pT = plt.bar(ind,TTots,width,bottom=[x+y for x,y in zip(ATots,CTots)],color=TColours,hatch='\\')
pG = plt.bar(ind,GTots,width,bottom=[x+y+z for x,y,z in zip(ATots,CTots,TTots)],color=GColours,hatch='X')

plt.ylabel('Percentage of locus (%)',fontsize=18)
plt.xlabel('Van de Harst SNPs',fontsize=18)
plt.title(InputInitial+", Nucleotide base at Van de Harst SNP locations",fontsize=20,fontweight='bold',y=1.03)
plt.xticks([i+(width/2) for i in ind], df_2['SNP'],rotation='vertical')
plt.yticks(np.arange(0,101,50))
plt.ylim(0,105)
plt.margins(0.01)
plt.tick_params(axis='both',direction='out',length=5,top='off',right='off')

snp_patch = mpatches.Patch(color='#F1911E',label='GWAS\nSNP base')
second_legend = plt.legend(handles=[snp_patch],loc=2,bbox_to_anchor=(1.01,0.8),fontsize=18)
plt.gca().add_artist(second_legend)
ax = plt.gca()
leg = ax.get_legend()
leg.legendHandles[0].set_edgecolor('black')

plt.legend((pA[0],pC[0],pT[0],pG[0]),('A','C','T','G'),loc=2,bbox_to_anchor=(1.01,1),fontsize=18)
ax = plt.gca()
leg = ax.get_legend()
i=0
while i<=3:
    leg.legendHandles[i].set_color('white')
    leg.legendHandles[i].set_edgecolor('black')
    i+=1
plt.subplots_adjust(bottom=0.15,right=0.87)
plt.show()
