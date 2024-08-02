import os
import pandas as pd
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
   # kaiju_summary = {}
    #with open(f'{sample}/{sample}_kaiju.names.out', 'r') as f:
    #with open(f'{sample}/{sample}_kaiju.names.out', 'r') as f:
      #  for line in f:
     #       fields = line.strip().split('\t')
       #     #contig_id = fields[1]
        #    taxonomy = fields[3]
         #   print(
  df1=pd.read_csv(f'{sample}/{sample}_kaiju.names.out',sep="\t",header=None)
  df1.columns=['Class','contigs','taxon_id','specie']
  df2=pd.read_csv(f'{sample}/{sample}_summary2.tsv',sep="\t")
  mylist=[]
  for i in range(10):
   mylist.append(df2['taxon_id'].iloc[i])
    #df1=df1[df1['taxon_id']==df2['taxon_id'].iloc[i]]
 #int(df2['taxon_id'].iloc[i])]
    #print(df2['taxon_id'].iloc[i])
    #print(df1)
  print(mylist)  
  with open(f'{sample}/{sample}_contigs.txt', 'w') as file:
   for i in mylist:
    df2=df1[df1['taxon_id']==i]
    df3=df2['contigs'].iloc[0]
    file.write(f"{df3}\n")
