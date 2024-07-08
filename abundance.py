import numpy as np
import matplotlib.pyplot as plt
from matplotlib import style 
import pyreadstat
import scipy.optimize as opt
plt.style.use('ggplot')
import rpy2.robjects as robjects
import pandas as pd
import statsmodels.api as sm
from scipy.stats import pearsonr
import csv
import glob
import os

# Define the directory containing the .tsv files
directory = '/Users/haroun/Soladata/output_kaiju'  # Change this to your directory
metadata=pd.read_csv('/Users/haroun/Soladata/SENTINEL_LASV_SEQ_SIJU_OA.csv', sep=",")

# Find all .tsv files in the specified directory
tsv_files = glob.glob(os.path.join(directory, '*.tsv'))

#df=pd.read_csv('/Users/haroun/Downloads/Coverage/')


dfs = []
samplename=[]
for tsv_file in tsv_files:
    
    print(f"Processing file: {tsv_file}")
    if os.path.getsize(tsv_file) == 0:
        print(f"Skipping {tsv_file} because it is empty.")
        continue

    # Read the .tsv file into a pandas DataFrame
    df= pd.read_csv(tsv_file,sep='\t',header=None)
    print(df)
    df.columns=['Taxon_id','counts']
    path= tsv_file.split('/')
    myfile=path[-1]
    parts = myfile.split('_')
    if len(parts) > 4:
      extracted_part = '_'.join(parts[:3])
      extracted_part0 = '_'.join(parts[:4])
      
      print(f"Extracted part after the second underscore: {extracted_part}")
    elif len(parts) ==4 :
      extracted_part = '_'.join(parts[:2])
      extracted_part0 = '_'.join(parts[:3]) 
    elif len(parts) ==3 :
      extracted_part = '_'.join(parts[:1])
      extracted_part0 = '_'.join(parts[:2])
    samplename.append(extracted_part0)
    with open('samplename.txt', 'w') as f:
      for line in samplename:
        f.write(f"{line}\n")
    print(f"Extracted part after the second underscore: {extracted_part}")  
    
   # print(metadata.columns)
    #metadata=metadata.loc[metadata["ID"] ==extracted_part, 'Sample Type'] 
 
       #df['SampleID']=extracted_part
      # df['Sample Type']=metadata['Sample Type']
      # df['ID']=metadata['ID']
    #print(df.head())
        
    metadata1=metadata.loc[metadata["ID"]==extracted_part]
     
    df['Sample Type']=metadata1['Sample Type'].values[0]
    df['Sample ID']=metadata1['ID'].values[0]
    df['Index']=metadata1['Index'].values[0]
    df['Location']=metadata1['Location'].values[0]
    df['Specie']=metadata1['Specie'].values[0]
    df['Sample ID']=metadata1['ID'].values[0]
    df=df[df['Taxon_id']!='unclassified Sphingomonas']
    print(df['Taxon_id'].values.tolist())
    dfs.append(df)
    

merged_df = pd.concat(dfs, ignore_index=True)
print(merged_df,'hello')


merged_df.to_csv("/Users/haroun/Soladata/merged_df.csv", index=False)

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import plotly.express as px
#px.bar(melt, x='value', y='Time', color='variable', orientation='h', text='value')

import seaborn as sns

"""
#df=pd.read_csv('Downloads/bacterialcounts.csv',sep=',')
df=pd.read_csv('/Users/haroun/Downloads/virals_samples.csv',sep='\t')

df =df.rename(columns={'STNL018':'Sample_ID','Plasma':'Sample_type','Human':'Specie', 'Ikorodu':'Location', 'counts':'Counts',  'species':'Viral_type'})
df['Counts']= df['Counts'].apply(pd.to_numeric, errors='coerce')

# Drop rows with NaN values (non-numerical)
df= df.dropna(subset=['Counts'])
print(df['Counts'].unique)
#df=df[df['Counts']!='counts']
#df=df[df['Counts']!='countss']
#df=df[df['Counts']!='count']
#df = df.dropna(subset=['Counts'])
print(df.head())
print(df['Counts'].sum())
grouped_sum = df.groupby(['Viral_type', 'Location'])['Counts'].sum().reset_index()
grouped_sum=grouped_sum[grouped_sum['Counts']>100]
print(grouped_sum)

 
fig=px.bar(grouped_sum,
    x='Counts', y='Viral_type', color='Location', orientation='h',title="Automatic Labels Based on Data Frame Column Names",width=1200, height=1200)
 
fig.update_layout(title_text='Viral abundance according to Location', title_x=0.5)
fig.show()
#fig.savefig('Viral abundance according to location',dpi=100)
"""
merged_df=merged_df.rename(columns={'Taxon_id':'Bacteria_Type'})
merged_df['Bacteria_Type'] =merged_df['Bacteria_Type'].str.replace(';', '', regex=True)
#merged_df=merged_df[merged_df['Location']=="Osun"]
grouped_sum = merged_df.groupby(['Bacteria_Type', 'Sample ID'])['counts'].sum().reset_index()
print(grouped_sum['counts'].sum())
grouped_sum=grouped_sum[grouped_sum['counts']>100]
print(grouped_sum)
grouped_sum.to_csv('grouped_sum.cst')

import random
from collections import defaultdict
import colorsys
import colorama

def generate_random_color():
    return "#{:06x}".format(random.randint(0, 0xFFFFFF))

def generate_random_color_list(n):
    return [generate_random_color() for _ in range(n)]
colordict=defaultdict(int)



N=len(grouped_sum['Bacteria_Type'].unique())
#random_colors = generate_random_color_list(N)

random_colors=generate_random_color_list(N)
for (bac,col) in zip(grouped_sum['Bacteria_Type'].unique(),random_colors):
  
   colordict[bac]=col
"""  
# Save the dictionary to a text file
with open('my_dict.txt', 'w') as f:
    for key, value in colordict.items():
        f.write(f"{key}: {value}\n")
"""
colordict = {}
with open('my_dict.txt', 'r') as f:
    for line in f:
        key, value = line.strip().split(': ')
        colordict[key] = value

#fig=px.bar(grouped_sum,
   # x='counts', y='Bacteria_Type', color='Sample ID', orientation='h',title="Automatic Labels Based on Data Frame Column Names")#width=1000, height=800
fig=px.bar(grouped_sum,
    x='Sample ID', y='counts', color='Bacteria_Type', color_discrete_map=colordict,title="Automatic Labels Based on Data Frame Column Names")#width=1000, height=800
fig.update_layout(xaxis=dict(tickfont=dict(size=10)),  # Adjust size as needed
                  yaxis=dict(tickfont=dict(size=9))) 
#fig.update_traces(width=0.8) 
"""
fig.update_layout(
    title='Bar Chart with Adjusted Bar Width',
    xaxis=dict(title='Category'),
    yaxis=dict(title='Values'),
    bargap=0.8  # This controls the gap between bars in the same category group
)
"""
fig.update_layout(title_text='Bacterial abundance in all samples', title_x=0.3)
fig.write_image("sampleloc.png")
fig.show()




    
    
    
  
