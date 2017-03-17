#GWAS
`Description`
The main functioning script is PlotSNPs_inputBam.py. This checks aligned read in a sorted bam file for the presence of the effect allele (EA) SNP from the van de Harst 75 loci for different blood phenotypes (van der Harst et al., Nature 2012. PMID: 23222517).
For each EA SNP that is present in the sequencing file, the script then checks how many imputed SNPs from linkage association are also present and outputs the SNP locations as a bed file.
