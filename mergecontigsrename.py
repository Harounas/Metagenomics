from pathlib import Path

input_dir = Path("./kraken_summary_files")

for fasta_path in input_dir.glob("*/contigs.fasta"):
    sample_name = fasta_path.parent.name
    output_path = fasta_path.with_name(f"{sample_name}_renamed.fasta")

    with open(fasta_path, "r") as infile, open(output_path, "w") as outfile:
        for line in infile:
            if line.startswith(">"):
                contig_name = line.strip()[1:]  # remove ">"
                new_header = f">{sample_name}|{contig_name}\n"
                outfile.write(new_header)
            else:
                outfile.write(line)

    print(f"Renamed contigs written to: {output_path.resolve()}")





