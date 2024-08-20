                                                                  bam2fastq.sh                                                                                  
while IFS= read -r line
do
    # Split the line into array to get paired-end files
    files=($line)
    #s0="${files%_*}" 
 
s1=$(basename $files .l1.AACK27LHV.1.cleaned.bam) 
echo $s1
samtools fastq $files -1 ${s1}_R1.fastq -2 ${s1}_R2.fastq
echo $s1>>samplenames0.txt
done < "bamfiles.txt"
sort samplenames0.txt|uniq>bamfilestrim1.txt
rm samplenames0.txt
