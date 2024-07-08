#!/bin/bash

# Specify the file containing the list of filenames
file_list="samplename.txt"

# Iterate over each line in the file
while IFS= read -r file
do
    echo "Processing $file"
sequences_file="${file}/${file}_contigs.txt"  # File containing sequence identifiers or file names
mkdir -p  Alignment
# Loop over each line in sequences_file
while IFS= read -r sequence_id; do
    input_fasta="${file}/scaffolds.fasta"  # Path to your input FASTA file
    output_fasta="${file}/${file}_${sequence_id}.fasta"  # Path to output subsequence FASTA file
   
    # Use seqtk subseq to extract subsequence based on range (e.g., 1-500)
    seqtk subseq "$input_fasta" <(echo "$sequence_id") > "$output_fasta"
     
    #echo "Extracted subsequence for $sequence_id"
    bwa index $output_fasta
    bwa mem -t 16 $output_fasta /home/harouna/240405_VH01476_5_AAC3HYYHV_batch2/${file}_R1_001.fastq.gz /home/harouna/240405_VH01476_5_AAC3HYYHV_batch2/${file}_R2_001.fastq.gz |samtools view  -u -@ 3 - | samtools sort -@ 16 -o Alignment/${file}_${sequence_id}_mapped.bam
samtools depth Alignment/${file}_${sequence_id}_mapped.bam -o Alignment/${file}_${sequence_id}_coverage.tsv
samtools sort Alignment/${file}_${sequence_id}_mapped.bam -o Alignment/${file}_${sequence_id}_mapped_sorted.bam
samtools index Alignment/${file}_${sequence_id}_mapped_sorted.bam
#samtools mpileup -A -d 1000 -f $output_fasta Alignment/AE_24_169_S20_${sequence_id}_mapped_sorted.bam > pileup.txt

#ivar consensus -p Alignment/AE_24_169_S20_${sequence_id}_consensus -t 0.5 -m 10 -n N -i <(samtools mpileup -A -d 1000 -f $output_fasta  Alignment/AE_24_169_S20_${sequence_id}_mapped_sorted.bam)
samtools mpileup -aa -A -d 0 -Q 0 Alignment/${file}_${sequence_id}_mapped_sorted.bam | ivar consensus -p  Alignment/${file}_${sequence_id}_consensus 
done < "$sequences_file"
    # Add your commands to process each file here
done < "$file_list"
