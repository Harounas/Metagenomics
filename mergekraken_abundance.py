import os
import pandas as pd
import plotly.express as px
import argparse
import random
from collections import defaultdict

# Argument parser to handle inputs
parser = argparse.ArgumentParser(description="Process Kraken summary files with metadata.")
parser.add_argument("kraken_dir", type=str, help="Path to the directory containing Kraken summary files.")
parser.add_argument("metadata_file", type=str, help="Path to the metadata CSV file.")
parser.add_argument("--virus", action='store_true', help="Focus on viruses.")
parser.add_argument("--bacteria", action='store_true', help="Focus on bacteria.")
parser.add_argument("--read_count", type=int, default=0, help="Minimum read count threshold.")
parser.add_argument("--top_N", type=int, default=None, help="Select the top N most common viruses or bacteria.")
args = parser.parse_args()

# Load metadata
metadata = pd.read_csv(args.metadata_file, sep=",")
sample_id_col = metadata.columns[0]  # Assume the first column is the sample ID

# Dictionary to store aggregated results
aggregated_results = {}

# Iterate over each summary output file
for file_name in os.listdir(args.kraken_dir):
    if file_name.endswith("_krakent.txt"):
        with open(os.path.join(args.kraken_dir, file_name), 'r') as f:
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
                
                if rank_code == 'S' and nr_frag_direct_at_taxon >= args.read_count:  # Only consider species-level entries with threshold
                    if extracted_part in metadata[sample_id_col].unique(): 
                        metadata1 = metadata.loc[metadata[sample_id_col] == extracted_part]
                        sample_metadata = metadata1.iloc[0].to_dict()

                        # Aggregate the results
                        aggregated_results[sampleandtaxonid] = {
                            'Perc_frag_cover': perc_frag_cover,
                            'Nr_frag_cover': nr_frag_cover,
                            'Nr_frag_direct_at_taxon': nr_frag_direct_at_taxon,
                            'Rank_code': rank_code,
                            'NCBI_ID': ncbi_ID,
                            'Scientific_name': scientific_name,
                            'SampleID': extracted_part,
                            **sample_metadata  # Unpack all metadata fields dynamically
                        }

# Output aggregated results to a TSV file
merged_tsv_path = os.path.join(args.kraken_dir, "merged_kraken1.tsv")
with open(merged_tsv_path, 'w') as f:
    # Write headers dynamically
    metadata_columns = metadata.columns.tolist()
    headers = ['Perc_frag_cover', 'Nr_frag_cover', 'Nr_frag_direct_at_taxon', 'Rank_code', 'NCBI_ID', 'Scientific_name', 'SampleID'] + metadata_columns[1:]
    f.write("\t".join(headers) + "\n")
    for sampleandtaxonid, data in aggregated_results.items():
        f.write("\t".join(str(data[col]) for col in headers) + "\n")

# Load the merged TSV file into a DataFrame
df = pd.read_csv(merged_tsv_path, sep="\t")
df.columns = df.columns.str.replace('/', '_').str.replace(' ', '_')
df = df.apply(lambda col: col.map(lambda x: x.strip() if isinstance(x, str) else x))
df = df[df['Scientific_name'] != 'Homo sapiens']  
# Conditionally filter the DataFrame for viruses or bacteria based on provided options
if args.virus:
    #df = df[df['Scientific_name'] != 'Homo sapiens']  # Remove human data
    df = df[df['Scientific_name'].str.contains('Virus', case=False, na=False)]  # Focus on viruses
    df = df.rename(columns={'Scientific_name': 'Virus_Type'})
elif args.bacteria:
    #df = df[df['Scientific_name'] != 'Homo sapiens']  # Remove human data
    df = df[~df['Scientific_name'].isin(
    df[df['Scientific_name'].str.contains('Virus', case=False, na=False)]['Scientific_name'].unique()
)]  # Focus on bacteria
    df = df.rename(columns={'Scientific_name': 'Bacteria_Type'})

# Apply top N filtering if the flag is set
if args.top_N:
    target_column = 'Virus_Type' if args.virus else 'Bacteria_Type'
    category_counts = df[target_column].value_counts()
    top_N_categories = category_counts.head(args.top_N).index
    df = df[df[target_column].isin(top_N_categories)]

# Identify categorical columns excluding 'Virus_Type' or 'Bacteria_Type'
target_column = 'Virus_Type' if args.virus else 'Bacteria_Type'
categorical_cols = df.select_dtypes(include=['object']).columns.tolist()
categorical_cols.remove(target_column)  # Exclude the target column as it's already used in the plots

# Function to generate a random color
def generate_random_color():
    return "#{:06x}".format(random.randint(0, 0xFFFFFF))

# Iterate over each categorical column and create a plot
for col in categorical_cols:
    grouped_sum = df.groupby([target_column, col])['Nr_frag_direct_at_taxon'].mean().reset_index()

    # Generate a color dictionary for the target categories
    colordict = defaultdict(int)
    N = len(grouped_sum[col].unique())
    random_colors = [generate_random_color() for _ in range(N)]
    for target, color in zip(grouped_sum[target_column].unique(), random_colors):
        colordict[target] = color

    # Adjust plot size and font size based on the number of unique target values
    plot_width = 1100 + 5 * N  # Increase width for better readability
    plot_height = 800 + 5* N
    font_size = max(10, 14 - N // 10)

    # Create the bar chart
    title_prefix = "Viral" if args.virus else "Bacterial"
    fig = px.bar(
        grouped_sum,
        x=col,
        y='Nr_frag_direct_at_taxon',
        color=target_column,
        color_discrete_map=colordict,
        title=f"{title_prefix} Abundance by {col}"
    )
    
    # Update layout for axis font size, legend, and figure dimensions
    fig.update_layout(
        xaxis=dict(
            tickfont=dict(size=font_size),
            tickangle=45  # Rotate x-axis labels by 45 degrees
        ),
        yaxis=dict(
            tickfont=dict(size=font_size)
        ),
        title=dict(
            text=f'Average {title_prefix} Abundance by {col}',
            x=0.5,
            font=dict(size=16)
        ),
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
        width=plot_width,  # Increased width
        height=plot_height
    )

    # Save the figure as an image file
    #fig.write_image(f"{title_prefix}Abundance_by_{col}.svg", format='svg')
    fig.write_image(f"{title_prefix}Abundance_by_{col}.png", format='png', scale=3)
    
