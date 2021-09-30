#!/bin/bash

cd data
for i in $(ls *.fastq)
do
	echo $i >> stat_fq_temp.txt
	cat $i | echo $(($(wc -l)/4)) >> stat_fq_temp.txt
done

perl -pe 's/fastq\n/fastq\t/' stat_fq_temp.txt > stat_fq.txt
                                                                    
rm *_temp.txt

#Variable that count the total number of reads
t=0
#Variable that count the total number of samples
c=0
cat stat_fq.txt | ( while read line
do
	id=$(echo $line | cut -d' ' -f1)
	nbseq=$(echo $line | cut -d' ' -f2)

	t=$(( $nbseq + t ))
	c=$(( c + 1 ))
done

#Compute the average number of reads from the sample pool
m=$((t / c))
echo "Average number of reads: " $m

#Compute 10% of the average
a=`bc <<< "scale=0; $m*0.1"`

#Round to integer the preceeding number
gf=$(echo "$a" | awk '{printf("%d\n",$1 + 0.5)}')

echo "10% of the average: " $gf

cat stat_fq.txt | ( while read line
do
	id=$(echo $line | cut -d' ' -f1)
	nbseq=$(echo $line | cut -d' ' -f2)
	
	test=""
	if [ "$gf" -gt "$nbseq" ]
		then test="Remove this sample from the analysis"
		mv $id ../reject/
	fi
	
	echo $id $nbseq $gf $test
done


)
)

mv stat_fq.txt ../stat
