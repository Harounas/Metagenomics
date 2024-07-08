while IFS= read -r line
do
    # Split the line into array to get paired-end files
    files=($line)
    s0="${files%_*}" 
 s1="${files%%_R[1-9]_*}"
 
echo $s1

#bowtie2 -p 16 -x GRCh38_noalt_as -1 ${s1}_R1_001.fastq.gz -2 ${s1}_R2_001.fastq.gz  --un-conc-gz ${s1}_unmapped_reads --out ${s1}_mapped_and_unmapped.sam
#rm ${s1}_mapped_and_unmapped.sam
#mv ${s1}_unmapped_reads.fastq.1.gz  ${s1}_R1_001_unmapped.fastq.gz
#mv ${s1}_unmapped_reads.fastq.2.gz  ${s1}_R2_001_unmapped.fastq.gz

#mkdir -p ${s1}_kaiju
#mv ${s1}_unmapped_reads.1 ${s1}_unmapped_reads.1.fastq.gz
#mv ${s1}_unmapped_reads.2 ${s1}_unmapped_reads.2.fastq.gz
#metaspades.py -1 ${s1}_R1_001_unmapped.fastq.gz -2 ${s1}_R2_001_unmapped.fastq>
#metaspades.py  -1 ${s1}_R1_001.fastq.gz  -2 ${s1}_R2_001.fastq.gz -o ${s1} -t 32
kaiju -t /home/harouna/Kaiju/nodes.dmp -f /home/harouna/Kaiju/kaiju_db_nr.fmi -i ${s1}_R1_001.fastq.gz -j ${s1}_R2_001.fastq.gz -o ${s1}_kaiju/${s1}_output  
kaiju2table -t /home/harouna/Kaiju/nodes.dmp -n /home/harouna/Kaiju/names.dmp -r genus -o ${s1}_kaiju/${s1}_summary1.tsv ${s1}_kaiju/${s1}_output 
#rm ${s1}_unmapped_reads.1.fastq.gz 
#rm ${s1}_unmapped_reads.2.fastq.gz
done < "fastq_files1.txt"
