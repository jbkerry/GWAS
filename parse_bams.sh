Counter=1
paireds=($(samtools view $sam_file | cut -f 1 | sort | uniq -c | grep '^\s*2' | awk -F' ' {'print $2'}))
for read in ${paireds[@]}; do
    samtools view -H $sam_file >temp.sam
    samtools view $sam_file | grep $read >>temp.sam 
    samtools view -h -b temp.sam | samtools sort -o temp.sorted.bam
    samtools index temp.sorted.bam
    samtools mpileup -f /databank/igenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa -u -v -A -o temp.vcf temp.sorted.bam
    grep -v \# temp.vcf | awk '$5!="<*>" {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t.\t.\t."}' | sed 's/,<[*]>//' >Pair_$Counter.vcf
    let Counter=Counter+1
done
