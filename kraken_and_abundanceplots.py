import os
import glob
import subprocess
import argparse
import pandas as pd
import random
from collections import defaultdict
import plotly.express as px

def run_trimmomatic(forward, reverse, base_name, output_dir, threads):
    trimmed_forward = os.path.join(output_dir, f"{base_name}_trimmed_R1.fastq.gz")
    trimmed_reverse = os.path.join(output_dir, f"{base_name}_trimmed_R2.fastq.gz")
    unpaired_forward = os.path.join(output_dir, f"{base_name}_unpaired_R1.fastq.gz")
    unpaired_reverse = os.path.join(output_dir, f"{base_name}_unpaired_R2.fastq.gz")

    if reverse:
        trimmomatic_cmd = [
            "trimmomatic", "PE", "-threads", str(threads),
            forward, reverse,
            trimmed_forward, unpaired_forward,
            trimmed_reverse, unpaired_reverse,
            "LEADING:3", "TRAILING:3", "SLIDINGWINDOW:4:15", "MINLEN:36"
        ]
    else:
        trimmed_forward = os.path.join(output_dir, f"{base_name}_trimmed_R1.fastq.gz")
        trimmomatic_cmd = [
            "trimmomatic", "SE", "-threads", str(threads),
            forward, trimmed_forward,
            "LEADING:3", "TRAILING:3", "SLIDINGWINDOW:4:15", "MINLEN:36"
        ]

    print("Running Trimmomatic command:", " ".join(trimmomatic_cmd))  # Debug
    subprocess.run(trimmomatic_cmd, check=True)
    
    return trimmed_forward, trimmed_reverse if reverse else None

def run_bowtie2(forward, reverse, base_name, bowtie2_index, output_dir, threads):
    unmapped_r1 = os.path.join(output_dir, f"{base_name}_unmapped_1.fastq.gz")
    unmapped_r2 = os.path.join(output_dir, f"{base_name}_unmapped_2.fastq.gz") if reverse else None

    bowtie2_cmd = [
        "bowtie2", "--threads", str(threads),
        "-x", bowtie2_index,
        "-1", forward, "-2", reverse,
        "--un-conc-gz", os.path.join(output_dir, f"{base_name}_unmapped_%.fastq.gz"),
        "-S", "/dev/null"
    ] if reverse else [
        "bowtie2", "--threads", str(threads),
        "-x", bowtie2_index,
        "-U", forward,
        "--un-gz", unmapped_r1,
        "-S", "/dev/null"
    ]

    print("Running Bowtie2 command:", " ".join(bowtie2_cmd))  # Debug
    subprocess.run(bowtie2_cmd, check=True)
    
    return unmapped_r1, unmapped_r2

def run_kraken2(forward, reverse, base_name, kraken_db, output_dir, threads):
    kraken_report = os.path.join(output_dir, f"{base_name}_report.txt")
    kraken_output = os.path.join(output_dir, f"{base_name}_kraken.txt")

    kraken_cmd = [
        "kraken2", "--db", kraken_db,
        "--threads", str(threads),
        "--report", kraken_report,
        "--output", kraken_output,
    ]

    if reverse:
        kraken_cmd.extend(["--paired", "--gzip-compressed", forward, reverse])
    else:
        kraken_cmd.extend(["--gzip-compressed", forward])

    print("Running Kraken2 command:", " ".join(kraken_cmd))  # Debug
    subprocess.run(kraken_cmd, check=True)

    return kraken_report

def process_sample(forward, reverse, base_name, bowtie2_index, kraken_db, output_dir, threads):
    # Step 1: Run Trimmomatic
    trimmed_forward, trimmed_reverse = run_trimmomatic(forward, reverse, base_name, output_dir, threads)
    
    # Step 2: Run Bowtie2 with trimmed reads to deplete host genome reads
    unmapped_r1, unmapped_r2 = run_bowtie2(trimmed_forward, trimmed_reverse, base_name, bowtie2_index, output_dir, threads)
    
    # Step 3: Use the unmapped reads from Bowtie2 as input for Kraken2
    kraken_report = run_kraken2(unmapped_r1, unmapped_r2, base_name, kraken_db, output_dir, threads)

    return kraken_report

def aggregate_kraken_results(kraken_dir, metadata_file, read_count, top_N, virus, bacteria):
    metadata = pd.read_csv(metadata_file, sep=",")
    sample_id_col = metadata.columns[0]  # Assume the first column is the sample ID

    # Dictionary to store aggregated results
    aggregated_results = {}

    # Iterate over each Kraken report file
    for file_name in os.listdir(kraken_dir):
        if file_name.endswith("_report.txt"):
            with open(os.path.join(kraken_dir, file_name), 'r') as f:
                for line in f:
                    fields = line.strip().split('\t')
                    perc_frag_cover = fields[0]
                    nr_frag_cover = fields[1]
                    nr_frag_direct_at_taxon = int(fields[2])
                    rank_code = fields[3]
                    ncbi_ID = fields[4]
                    scientific_name = fields[5]
                    parts = file_name.split('_')
                    extracted_part = '_'.join(parts[:-1])
                    sampleandtaxonid = extracted_part + str(ncbi_ID)

                    if rank_code == 'S' and nr_frag_direct_at_taxon >= read_count:
                        if extracted_part in metadata[sample_id_col].unique():
                            sample_metadata = metadata.loc[metadata[sample_id_col] == extracted_part].iloc[0].to_dict()
                            aggregated_results[sampleandtaxonid] = {
                                'Perc_frag_cover': perc_frag_cover,
                                'Nr_frag_cover': nr_frag_cover,
                                'Nr_frag_direct_at_taxon': nr_frag_direct_at_taxon,
                                'Rank_code': rank_code,
                                'NCBI_ID': ncbi_ID,
                                'Scientific_name': scientific_name,
                                'SampleID': extracted_part,
                                **sample_metadata
                            }

    # Output aggregated results to a TSV file
    merged_tsv_path = os.path.join(kraken_dir, "merged_kraken1.tsv")
    with open(merged_tsv_path, 'w') as f:
        # Write headers dynamically
        headers = ['Perc_frag_cover', 'Nr_frag_cover', 'Nr_frag_direct_at_taxon', 'Rank_code', 'NCBI_ID', 'Scientific_name', 'SampleID'] + metadata.columns[1:].tolist()
        f.write("\t".join(headers) + "\n")
        for sampleandtaxonid, data in aggregated_results.items():
            f.write("\t".join(str(data[col]) for col in headers) + "\n")

    return merged_tsv_path

def generate_abundance_plots(merged_tsv_path, virus, bacteria, top_N):
    df = pd.read_csv(merged_tsv_path, sep="\t")
    df.columns = df.columns.str.replace('/', '_').str.replace(' ', '_')
    df = df.apply(lambda col: col.map(lambda x: x.strip() if isinstance(x, str) else x))
    df = df[df['Scientific_name'] != 'Homo sapiens']

    if virus:
        df = df[df['Scientific_name'].str.contains('Virus', case=False, na=False)]
        df = df.rename(columns={'Scientific_name': 'Virus_Type'})
    elif bacteria:
        df = df[~df['Scientific_name'].str.contains('Virus', case=False, na=False)]
        df = df.rename(columns={'Scientific_name': 'Bacteria_Type'})

    if top_N:
        target_column = 'Virus_Type' if virus else 'Bacteria_Type'
        top_N_categories = df[target_column].value_counts().head(top_N).index
        df = df[df[target_column].isin(top_N_categories)]

    target_column = 'Virus_Type' if virus else 'Bacteria_Type'
    categorical_cols = df.select_dtypes(include=['object']).columns.tolist()
    categorical_cols.remove(target_column)

    for col in categorical_cols:
        grouped_sum = df.groupby([target_column, col])['Nr_frag_direct_at_taxon'].mean().reset_index()

        colordict = defaultdict(int)
        random_colors = ["#{:06x}".format(random.randint(0, 0xFFFFFF)) for _ in range(len(grouped_sum[col].unique()))]
        for target, color in zip(grouped_sum[target_column].unique(), random_colors):
            colordict[target] = color

        plot_width = 1100 + 5 * len(grouped_sum[col].unique())
        plot_height = 800 + 5 * len(grouped_sum[col].unique())
        font_size = max(10, 14 - len(grouped_sum[col].unique()) // 10)

        title_prefix = "Viral" if virus else "Bacterial"
        fig = px.bar(
            grouped_sum,
            x=col,
            y='Nr_frag_direct_at_taxon',
            color=target_column,
            color_discrete_map=colordict,
            title=f"{title_prefix} Abundance by {col}"
        )

        fig.update_layout(
            xaxis=dict(tickfont=dict(size=font_size), tickangle=45),
            yaxis=dict(tickfont=dict(size=font_size)),
            title=dict(text=f'Average {title_prefix} Abundance by {col}', x=0.5, font=dict(size=16)),
            bargap=0.5,
            legend=dict(
                font=dict(size=font_size),
                x=1,
                y=1,
                traceorder='normal',
                orientation='v',
                itemwidth=30,
                itemsizing='constant',
                itemclick='toggleothers',
                itemdoubleclick='toggle'
            ),
            width=plot_width,
            height=plot_height
        )

        fig.write_image(f"{title_prefix}Abundance_by_{col}.png", format='png', scale=3)

def main():
    parser = argparse.ArgumentParser(description="Pipeline for Trimmomatic trimming, Bowtie2 host depletion, and Kraken2 taxonomic classification.")
    parser.add_argument("--kraken_db", required=True, help="Path to Kraken2 database.")
    parser.add_argument("--bowtie2_index", required=True, help="Path to Bowtie2 index.")
    parser.add_argument("--output_dir", required=True, help="Directory to save output files.")
    parser.add_argument("--input_dir", required=True, help="Directory containing input FASTQ files.")
    parser.add_argument("--threads", type=int, default=8, help="Number of threads to use for Trimmomatic, Bowtie2, and Kraken2.")
    parser.add_argument("--metadata_file", required=True, help="Path to the metadata CSV file.")
    parser.add_argument("--read_count", type=int, default=0, help="Minimum read count threshold.")
    parser.add_argument("--top_N", type=int, default=None, help="Select the top N most common viruses or bacteria.")
    parser.add_argument("--virus", action='store_true', help="Focus on viruses.")
    parser.add_argument("--bacteria", action='store_true', help="Focus on bacteria.")

    args = parser.parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    # Step 1: Process each sample
    for forward in glob.glob(os.path.join(args.input_dir, "*_R1.fastq*")):
        base_name = os.path.basename(forward).replace("_R1.fastq.gz", "").replace("_R1.fastq", "")
        reverse = os.path.join(args.input_dir, f"{base_name}_R2.fastq.gz") if forward.endswith(".gz") else os.path.join(args.input_dir, f"{base_name}_R2.fastq")

        if not os.path.isfile(reverse):
            reverse = None

        process_sample(forward, reverse, base_name, args.bowtie2_index, args.kraken_db, args.output_dir, args.threads)

    # Step 2: Aggregate Kraken results
    merged_tsv_path = aggregate_kraken_results(args.output_dir, args.metadata_file, args.read_count, args.top_N, args.virus, args.bacteria)

    # Step 3: Generate abundance plots
    generate_abundance_plots(merged_tsv_path, args.virus, args.bacteria, args.top_N)

if __name__ == "__main__":
    main()

