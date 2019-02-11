import numpy as np
import pandas as pd
import sys
import csv
import time
import more_itertools as mit

class IntersTable():
    def __init__(self, strand_pos, strand_neg):
        self.pos = strand_pos
        self.neg = strand_neg


    def intersect_ns(self, assembly, distance):
        self.intersecting = []
        self.non_intersecting = []

        for index, row in assembly.df.iterrows():
            intersection = set()
            xs = set(range(row["start"] - distance, row["stop"] + distance))
            if not intersection.union(xs.intersection(self.neg)) and not intersection.union(xs.intersection(self.pos)):
                self.non_intersecting.append((row['seqid'], row['source'], row['type'], row['start'], row['stop'], row['score'], row['strand'], row['phase'], row['attributes']))
                continue


            else:
                self.intersecting.append((row['seqid'], row['source'], row['type'], row['start'], row['stop'], row['score'], row['strand'], row['phase'], row['attributes']))
                continue


        return self

    def intersect_ss(self, assembly, distance):
        self.intersecting = []
        self.non_intersecting = []
        for index, row in assembly.df.iterrows():
            intersection = set()
            xs = set(range(row["start"] - distance, row["stop"] + distance))
            if row["strand"] == "+" and not intersection.union(xs.intersection(self.pos)):
                self.non_intersecting.append((row['seqid'], row['source'], row['type'], row['start'], row['stop'], row['score'], row['strand'], row['phase'], row['attributes']))

            elif row["strand"] == "-" and not intersection.union(xs.intersection(self.neg)):
                self.non_intersecting.append((row['seqid'], row['source'], row['type'], row['start'], row['stop'], row['score'], row['strand'], row['phase'], row['attributes']))

            else:
                self.intersecting.append((row['seqid'], row['source'], row['type'], row['start'], row['stop'], row['score'], row['strand'], row['phase'], row['attributes']))

        return self

    def merge(self):
        self.occup_pos = set()
        self.occup_neg = set()
        for elem in self.non_intersecting:
            if elem[6] == "+":
                self.occup_pos.update(range(elem[3], elem[4]))
            elif elem[6] == "-":
                self.occup_neg.update(range(elem[3], elem[4]))


        self.merged_pos = []
        self.merged_neg = []
        i = 0
        for group in mit.consecutive_groups(sorted(self.occup_pos)):
            intermediate = list(group)
            self.merged_pos.append(["---", "TIM_Analysis", "exon", min(intermediate), max(intermediate), "0","+",".","IGR_"+str(i)])
            i+=1

        i = 0
        for group in mit.consecutive_groups(sorted(self.occup_neg)):
            #self.merged_neg.append(list(group))
            self.merged_neg.append(["---", "TIM_Analysis", "exon", min(intermediate), max(intermediate), "0", "-", ".", "IGR_" + str(i)])
            i += 1
        self.merged_final = self.merged_pos
        self.merged_final.extend(self.merged_neg)




















