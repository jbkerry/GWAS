### Given an input sam file this will count the total number of paired and single reads ###

paireds=$(samtools view $1 | cut -f 1 | sort | uniq -c | grep -c '^\s*2')
paired_num=$(($paireds*2))
singles=$(samtools view $1 | cut -f 1 | sort | uniq -c | grep -c '^\s*1')
echo "$paireds read pairs ($paired_num reads found in pairs). $singles single reads."
