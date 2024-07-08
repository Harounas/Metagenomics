# Example script to link bacteria with AMR genes across multiple samples
import os
lines = []
with open('samplename.txt', 'r') as f:
    for line in f:
        lines.append(line.strip())

print(lines)
# List of samples (adjust as per your directory structure)
samples = lines

# Dictionary to store linked results
linked_results = {}

for sample in samples:
    # Read Kaiju output
    kaiju_summary = {}
    with open(f'{sample}_kaiju.names.out', 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            contig_id = fields[1]
            taxonomy = fields[3]
            kaiju_summary[contig_id] = taxonomy
           
    # Read AMRFinderPlus output
    amr_genes = {}
    with open(f'/Users/haroun/AMR/{sample}_amr.txt', 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            contig_id = fields[1]
            gene_symbol = fields[6]
            if contig_id in amr_genes:
                amr_genes[contig_id].append(gene_symbol)
            else:
                amr_genes[contig_id] = [gene_symbol]
    print(amr_genes)
    # Link bacteria with AMR genes for the current sample
    for contig_id in kaiju_summary:
        if contig_id in amr_genes:
            bacteria = kaiju_summary[contig_id]
            if bacteria in linked_results:
                linked_results[bacteria].extend(amr_genes[contig_id])
            else:
                linked_results[bacteria] = amr_genes[contig_id]

# Output linked results
for bacteria, genes in linked_results.items():
    print(f"Bacteria: {bacteria}")
    print(f"Associated AMR Genes: {', '.join(set(genes))}")  # Use set to remove duplicates
   

