from pathlib import Path

input_paths = sorted(Path("./kraken_summary_files").glob("*/contigs.fasta"))
output_file = Path("merged_contigs.fasta")

def process_and_merge(contigs_files, output_path):
    with open(output_path, "w") as outfile:
        for path in contigs_files:
            sample_name = path.parent.name  # extract sample name from path like kraken_summary_files/{sample}/contigs.fasta
            with open(path, "r") as infile:
                for line in infile:
                    if line.startswith(">"):
                        contig_name = line.strip()[1:]  # remove ">"
                        new_header = f">{sample_name}|{contig_name}\n"
                        outfile.write(new_header)
                    else:
                        outfile.write(line)

process_and_merge(input_paths, output_file)
print(f"Merged contigs written to: {output_file.resolve()}")

