import os
import subprocess
import argparse
import glob

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
    if reverse:
        unmapped_r1 = os.path.join(output_dir, f"{base_name}_unmapped_1.fastq.gz")
        unmapped_r2 = os.path.join(output_dir, f"{base_name}_unmapped_2.fastq.gz")

        bowtie2_cmd = [
            "bowtie2", "--threads", str(threads),
            "-x", bowtie2_index,
            "-1", forward, "-2", reverse,
            "--un-conc-gz", os.path.join(output_dir, f"{base_name}_unmapped_%.fastq.gz"),
            "-S", "/dev/null"
        ]
        print("Running Bowtie2 command:", " ".join(bowtie2_cmd))  # Debug
        subprocess.run(bowtie2_cmd, check=True)
        
        return unmapped_r1, unmapped_r2
    else:
        unmapped_r1 = os.path.join(output_dir, f"{base_name}_unmapped_1.fastq.gz")

        bowtie2_cmd = [
            "bowtie2", "--threads", str(threads),
            "-x", bowtie2_index,
            "-U", forward,
            "--un-gz", unmapped_r1,
            "-S", "/dev/null"
        ]
        print("Running Bowtie2 command:", " ".join(bowtie2_cmd))  # Debug
        subprocess.run(bowtie2_cmd, check=True)
        
        return unmapped_r1, None

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

def process_sample(forward, reverse, base_name, bowtie2_index, kraken_db, output_dir, threads):
    # Step 1: Run Trimmomatic
    trimmed_forward, trimmed_reverse = run_trimmomatic(forward, reverse, base_name, output_dir, threads)
    
    # Step 2: Run Bowtie2 with trimmed reads to deplete host genome reads
    unmapped_r1, unmapped_r2 = run_bowtie2(trimmed_forward, trimmed_reverse, base_name, bowtie2_index, output_dir, threads)
    
    # Step 3: Use the unmapped reads from Bowtie2 as input for Kraken2
    run_kraken2(unmapped_r1, unmapped_r2, base_name, kraken_db, output_dir, threads)

def main():
    parser = argparse.ArgumentParser(description="Pipeline for Trimmomatic trimming, Bowtie2 host depletion, and Kraken2 taxonomic classification.")
    parser.add_argument("--kraken_db", required=True, help="Path to Kraken2 database.")
    parser.add_argument("--bowtie2_index", required=True, help="Path to Bowtie2 index.")
    parser.add_argument("--output_dir", required=True, help="Directory to save output files.")
    parser.add_argument("--input_dir", required=True, help="Directory containing input FASTQ files.")
    parser.add_argument("--threads", type=int, default=8, help="Number of threads to use for Trimmomatic, Bowtie2, and Kraken2.")
    
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    # Loop through all forward read files
    for forward in glob.glob(os.path.join(args.input_dir, "*_R1.fastq*")):
        extension = forward.split('.')[-1]
        if extension == "gz":
            base_name = os.path.basename(forward).replace("_R1.fastq.gz", "")
            reverse = os.path.join(args.input_dir, f"{base_name}_R2.fastq.gz")
        else:
            base_name = os.path.basename(forward).replace("_R1.fastq", "")
            reverse = os.path.join(args.input_dir, f"{base_name}_R2.fastq")

        if os.path.isfile(reverse):
            process_sample(forward, reverse, base_name, args.bowtie2_index, args.kraken_db, args.output_dir, args.threads)
        else:
            process_sample(forward, None, base_name, args.bowtie2_index, args.kraken_db, args.output_dir, args.threads)

if __name__ == "__main__":
    main()
