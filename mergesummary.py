import os
import pandas as pd
"""
# Directory containing Kaiju summary output files
output_dir = "/Users/haroun/Soladata/output_kaiju"

# Dictionary to store aggregated results
aggregated_results = {}

# Iterate over each summary output file
for file_name in os.listdir(output_dir):
    if file_name.endswith("_kaiju.names.out"):
        with open(os.path.join(output_dir, file_name), 'r') as f:
            for line in f:
                #if line.startswith('C'):
                    #continue
                fields = line.strip().split('\t')
                taxon_id = fields[1]
                taxon_name = fields[3]
                read_count = int(fields[2])
                parts = file_name.split('_')
                if len(parts) > 4:
                  extracted_part = '_'.join(parts[:3])
           
      
                elif len(parts) ==4 :
                   extracted_part = '_'.join(parts[:2])
                elif len(parts) ==3 :
                   extracted_part = '_'.join(parts[:1])
     
                if taxon_id in aggregated_results:
                    aggregated_results[taxon_id]['read_count'] += read_count
                else:
                    aggregated_results[taxon_id] = {
                        'taxon_name': taxon_name,
                        'read_count': read_count,
                        'Sample ID': extracted_part
                    }

# Output aggregated results
with open(os.path.join(output_dir, "aggregated_kaiju_summary.out"), 'w') as f:
    f.write("Taxon ID\tTaxon Name\tRead Count\tSample ID\n")
    for taxon_id, data in aggregated_results.items():
        f.write(f"{taxon_id}\t{data['taxon_name']}\t{data['read_count']}\t{data['Sample ID']}\n")
"""
df=pd.read_csv("aggregated_kaiju_summary.out",sep="\t")
print(df)
cov=[]
Len=[]
for col in df['Taxon ID']:
  parts = col.split('_')
  cov.append(parts[5])
  Len.append(parts[3])


df['Length']=Len
df['Cov']=cov
#print(max(Len))
df=df[df['Taxon Name']=='Sphingomonas sp. TF3;']
#max_row = df.loc[df['Length'].idxmax()]
df['Length'] = pd.to_numeric(df['Length'], errors='coerce')
df = df.dropna(subset=['Length'])
max_row = df.loc[df['Length'].idxmax()]
print(df['Read Count'].max())
print(max_row)
