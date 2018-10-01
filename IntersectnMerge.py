#!/usr/bin/env python
import argparse
import requests
import pandas as pd
#import numpy as np
import sys

parser = argparse.ArgumentParser(description='Input Files')
parser.add_argument('--assembly', metavar='A', type=str, nargs=1,
                   help='The mapped assembly in gtf format (filepath from pwd)')
parser.add_argument('--genome', metavar='G', type=str, nargs=1,
                   help='The annotated reference genome in gff3 format (filepath from pwd)')
parser.add_argument('--distance', metavar='D', type=int, nargs=1, default=10,
                   help='Distance of kept feature to annotated features')
parser.add_argument('--output_path', metavar='O', type=str, nargs=1,
                   help='Path of the output files')
parser.add_argument('--project_prefix', metavar='N', type=str, nargs=1,
                   help='Project name, all output files will be prefixed with this')

args = parser.parse_args()
assembly=args.assembly[0]
genome=args.genome[0]
distance=args.distance[0]+1
output=args.output_path[0]
projectname=args.project_prefix[0]
col_names = ['seqid', 'source', 'type', 'start', 'stop', 'score', 'strand', 'phase', 'attributes']

#Read in gff files to be compared as pandas dataframes
df1 = pd.read_csv(assembly, sep='\t', comment='#', low_memory=False, header=None, names=col_names)
df1 = df1.sort_values(["strand", "start", "stop"], axis=0).reset_index(drop=True)
df1 = df1.drop_duplicates(subset=["start", "stop", "strand"])
df1 = df1.reset_index(drop=True)
df1.loc[df1['type'] != "exon", 'type'] = "exon"

df2 = pd.read_csv(genome, sep='\t', comment='#', low_memory=False, header=None, names=col_names)
df2 = df2.sort_values(["strand", "start", "stop"], axis=0).reset_index(drop=True)
df2 = df2.drop_duplicates(subset=["start", "stop", "strand"])
df2 = df2.reset_index(drop=True)
df2.loc[df2['type'] != "exon", 'type'] = "exon"

def iterate_df2(df1,df2, df1_chrm, df1_range, df1_row, intersecting_df, non_intersecting_df):
    number_intersections = 0
    intersection = set()
    for index, row in df2.iterrows():
        df2_chrm = row["strand"]
#Only compare, if the two reads are on the same chromosome
        df2_range = range(row["start"],row["stop"])
        xs = set(df1_range)
        intersection = intersection.union(xs.intersection(df2_range))
        if intersection:
            number_intersections = 1
            return number_intersections
        else:
            continue
    return number_intersections

def iterate_df1(df1,df2,distance):
    intersecting_df = pd.DataFrame(columns=df1.columns, index=[])
    non_intersecting_df = pd.DataFrame(columns=df1.columns, index=[])
    for index, row in df1.iterrows():
        sys.stdout.write('\r' + "Processing Assembly Row: " + str(row.name+1) + "/" + str(len(df1.index)))
        sys.stdout.flush()
#Create a range of the start and end point of the feature
        df1_range = range(row["start"] - distance,row["stop"] + distance)
        df1_row = row.name
        df1_chrm = row["strand"]
#Create Dataset for intersections and iterate over df1 Dataframe
        number_intersections = iterate_df2(df1,df2, df1_chrm, df1_range, df1_row,intersecting_df,non_intersecting_df)
        if int(number_intersections) == 0:
            non_intersecting_df = non_intersecting_df.append(df1.iloc[df1_row].copy())
        elif int(number_intersections) == 1:
            intersecting_df = intersecting_df.append(df1.iloc[df1_row].copy())
    return non_intersecting_df, intersecting_df
#Function to merge two Dataframes
def merge(merge1, merge2):
    merged_df = pd.DataFrame(columns=merge1.columns, index=[])
    lenm1 = len(merge1.index)
    for index, row in merge1.iterrows():
        trans_start = [row["start"]]
        trans_end = [row["stop"]]
        trans_row = row.name
        trans_strand=row["strand"]
        transcript_range=range(row["start"],row["stop"])
        sys.stdout.write('\r'+ "Processing feature: " + str(trans_row+1) + "/" + str(lenm1))
        sys.stdout.flush()
#        print("Processing non-intersecting feature: " + str(trans_row+1) + "/" + str(lenm1))
        for index, row in merge2.iterrows():
            gene_strand=row["strand"]
            if gene_strand == trans_strand:
                gene_range=range(row["start"],row["stop"])
                xs = set(transcript_range)
                intersection = set()
                intersection=xs.intersection(gene_range)
                if intersection:
#                    print("Inside for loop: " + str(trans_start))
                    trans_start.append(row["start"])
                    trans_end.append(row["stop"])
                    del intersection
#        merge1.set_value(trans_row, "start", min(trans_start))
        merge1.at[trans_row, "start"] = min(trans_start)
#        print("Outside for loop: " + str(trans_start))
#        merge1.set_value(trans_row, "stop", max(trans_end))
        merge1.at[trans_row, "stop"] = max(trans_end)
        merged_df = merged_df.append(merge1.iloc[trans_row].copy())
    return merged_df



#Create an empty dataframe for intersecting and nonintersecting reads each
#Put out the lengths of the files read in
print("Number of Features in Assembly File: " + str(len(df1.index)))
print("Number of Features in Genome File: " + str(len(df2.index)))

#Iterate over RAD21 dataframe
non_intersecting_df,intersecting_df = iterate_df1(df1,df2,distance)
del df1
del df2
#Remove duplicates and save results as csv
non_intersecting_df = non_intersecting_df.drop_duplicates(subset=["start", "stop", "strand", "attributes"])
non_intersecting_df = non_intersecting_df.reset_index(drop=True)
intersecting_df = intersecting_df.drop_duplicates(subset=["start", "stop", "strand", "attributes"])
intersecting_df = intersecting_df.reset_index(drop=True)
non_intersecting_df.to_csv(output + projectname + "non-intersecting.csv", sep='\t', index=False, header=False)
intersecting_df.to_csv(output + projectname + "_intersecting_discarded.csv", sep='\t', index=False, header=False)
del intersecting_df

####BEGIN MERGING ####

merge1 = non_intersecting_df
merge2 = non_intersecting_df
del non_intersecting_df
print("\nMerging remaining " + str(len(merge1.index)) + " Features")
merged_df = merge(merge1, merge2)
duplicates_merged = merged_df
duplicates_merged.to_csv(output + projectname + "_merged_with-duplicates.csv", sep='\t', index=False, header=False)

merged_df = merged_df.drop_duplicates(subset=["start", "stop", "strand"])
merged_df = merged_df.reset_index(drop=True)
print("\n" + str(len(merged_df.index)) + " Features remaining. Searching for duplicated Features...")

merge1 = merged_df
merge2 = merged_df
merged_df = merge(merge1, merge2)
del merge1
del merge2
merged_df = merged_df.drop_duplicates(subset=["start", "stop", "strand"])
merged_df = merged_df.reset_index(drop=True)
merged_df = merged_df.sort_values(["start", "stop"], axis=0).reset_index(drop=True)
merged_df.loc[merged_df['type'] != "exon", 'type'] = "exon"
merged_df.to_csv(output + projectname + "_final_merged.gff", sep='\t', index=False, header=False)
print("\nDone. " + str(len(merged_df.index)) + " Features remaining. Final output file was written to: " + output + projectname + "_final_merged.csv")
