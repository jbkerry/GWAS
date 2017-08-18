sam_lines=$(wc -l < $sam_file)
head_lines=$(samtools view -H $sam_file | wc -l)
LowCount=$(($head_lines+1))
HighCount=$(($LowCount+1))
MaxCount=$(($sam_lines+1))
while [ $HighCount -lt $MaxCount ]; do
    readName=$(awk 'FNR=='$LowCount' {print $1}' $sam_file)
    readNum=$(awk 'FNR>='$LowCount' && FNR<='$HighCount $sam_file | grep $readName | wc -l)
    if [ $readNum -eq 1 ]
    then
        let LowCount=LowCount+1
        let HighCount=HighCount+1
    elif [ $readNum -eq 2 ]
    then
        samtools view -H $sam_file >test.sam
        awk 'FNR>='$LowCount' && FNR<='$HighCount' {print}' $sam_file >>test.sam
        samtools view -h -b test.sam | samtools sort -o test.sorted.bam
        samtools index test.sorted.bam
        samtools mpileup -f /databank/igenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa -u -v -A -o test.vcf test.sorted.bam
        grep -v \# test.vcf | awk '$5!="<*>" {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t.\t.\t."}' | sed 's/,<[*]>//' >test_parse_$LowCount.vcf
        let HighCount=HighCount+2
        let LowCount=LowCount+2
    else
	echo "error"
    fi
done
