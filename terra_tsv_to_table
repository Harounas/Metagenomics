import io
import os
import pandas as pd
from google.cloud import storage
df=pd.read_csv('/Users/haroun/Downloads/sample_Terra.tsv', sep="\t")
storage_client = storage.Client()

# 2. Get list of GCS URLs from your DataFrame
# Assuming you have a DataFrame 'df' with 'assembly_stats_by_taxon_tsv' column
urls = df['assembly_stats_by_taxon_tsv'].dropna().astype(str).tolist()
def download_and_concat_tsvs(gcs_urls):
    all_data = []
    header = None
    
    for idx, url in enumerate(gcs_urls):
        if not url or url.lower() == "nan":
            continue
            
        if not url.startswith('gs://'):
            print(f"Skipping invalid URL: {url}")
            continue

        # Parse GCS URL
        _, _, bucket_name, *blob_parts = url.split('/')
        blob_path = '/'.join(blob_parts)
        
        # Get blob content
        bucket = storage_client.bucket(bucket_name)
        blob = bucket.blob(blob_path)
        content = blob.download_as_text()
        
        # Process TSV content
        lines = content.split('\n')
        if not lines:
            continue
            
        # Handle header
        if idx == 0:
            header = lines[0]
            all_data.append(content)
        else:
            if lines[0] != header:
                print(f"Header mismatch in {url} - using original header")
            all_data.append('\n'.join(lines[1:]))
    
    # Create combined TSV using io.StringIO
    combined_tsv = '\n'.join(all_data)
    return pd.read_csv(io.StringIO(combined_tsv), sep='\t')

# Usage remains the same
urls = df['assembly_stats_by_taxon_tsv'].dropna().astype(str).tolist()
storage_client = storage.Client()
combined_df = download_and_concat_tsvs(urls)
combined_df.to_csv("combined_assembly_stats.tsv", sep='\t', index=False)
