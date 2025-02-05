#!/usr/bin/env python
# coding: utf-8

# Python 3.12.8, pandas 2.2.3

import pandas as pd

# Read in bed
df = pd.read_csv("catLiftOffGenesV1.bed", sep='\t', header=None, low_memory=False)

# Split csv string into list of exon start sites and lengths
df['exons_list'] = df[11].apply(lambda y: [int(x) for x in y.split(',')])
df['exons_length'] = df[10].apply(lambda y: [int(x) for x in y.split(',')])

# Keep the relevant columns
df2 = df[[3,0,5,1,2,'exons_list','exons_length',19,31]].copy()

# Convert transcript based exon positions to genome based exon positions by adding the transcript start to each exon start
df2['exon_start'] = df2.apply(lambda row: [row[1] + x for x in row['exons_list']], axis=1)

# Get exon end positions by adding the exon length to the exon start
df2['exon_end'] = df2.apply(lambda row: [x + y for x, y in zip(row['exon_start'], row['exons_length'])], axis=1)

# Remove uneeded things
df2[0] = df2[0].str.replace("chr","")
df2=df2.drop(['exons_list','exons_length'],axis=1)

# Rename columns and convert lists to strings
df2.columns = ['#NAME', 'CHROM', 'STRAND', 'TX_START', 'TX_END', 'type','CAT_class', 'EXON_START', 'EXON_END']
df2['EXON_START'] = df2['EXON_START'].apply(lambda x: ','.join(map(str, x)))
df2['EXON_END'] = df2['EXON_END'].apply(lambda x: ','.join(map(str, x)))

# Output to file
df2.to_csv('chm13_spliceAI_annotation_unfiltered.txt', sep='\t', index=False)
df2[df2['type']=='protein_coding'].drop(['type','CAT_class'],axis=1).to_csv('chm13_spliceAI_annotation_protein_coding.txt', sep='\t', index=False)
df2.drop(['type','CAT_class'],axis=1).to_csv('chm13_spliceAI_annotation_all.txt', sep='\t', index=False)

