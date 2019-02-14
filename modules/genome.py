import numpy as np
import pandas as pd
import sys
import csv
import time

class Genome():
    def __init__(self):
        self.df = pd.DataFrame()

    def load_genome(self, genome_path, col_names):
        start_time = time.time()
        df = pd.read_csv(genome_path, sep='\t', comment='#', low_memory=False, header=None, names=col_names)
        df = df.sort_values(["strand", "start", "stop"], axis=0).reset_index(drop=True)
        df = df.drop_duplicates(subset=["start", "stop", "strand"])
        df = df.reset_index(drop=True)
        df.loc[df['type'] != "exon", 'type'] = "exon"
        self.df = df
        print(f"Finished reading in Genome in {(time.time() - start_time) * 1000} milliseconds.\nFound {len(self.df.index)} Features in Genome\n")
        return self

    def calc_intersect(self):
        start_time = time.time()
        strand_pos = set()
        strand_neg = set()
        for index, row in self.df.iterrows():
            if row['strand'] == "+":
                for elem in range(row['start'], row['stop']):
                    strand_pos.add(elem)
            elif row['strand'] == "-":
                for elem in range(row['start'], row['stop']):
                    strand_neg.add(elem)
#        max_len = max(strand_pos, strand_neg)
#        strand_pos = strand_pos + [0] * (len(max_len) - len(strand_pos))
#        strand_pos = strand_pos + [0] * (len(max_len) - len(strand_neg))
#        occupied = pd.DataFrame(np.column_stack([strand_pos, strand_neg]), columns=['+', '-'])
#        intersect.df[row["strand"]] = df.apply(lambda row: (row['start']))
#        intersection_places.append(range(row['start'], row['stop'])
#        self.df['intersections'] = self.df.apply(lambda row: (list(range(row['start'], row['stop']))), axis=1)
        print("Finished calculating intersections in "+ str((time.time() - start_time)*1000) + " milliseconds. \nStarting Intersection Analysis...")

        return strand_pos, strand_neg

    def max_size(self):
        size_calculation = set()
        df = self.df
        for index, row in df.iterrows():
            size_calculation.add(row["stop"] - row["start"])
        longest_feature = max(size_calculation)
        print("Longest Feature (search radius is length + 10): " + str(longest_feature) + " bp")
        return longest_feature
