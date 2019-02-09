import pandas as pd
import sys

def intersect(assembly, occupied_pos, occupied_neg):
    count = 0
    for index, row in assembly.df.iterrows():
        occur_count = 0
        if row["strand"] == "+":
            for elem in range(row['start']-10, row['stop']+10):
                if elem in occupied_pos or occupied_neg:
                    occur_count += 1
                else:
                    continue
        elif row["strand"] == "-":
            for elem in range(row['start'] - 10, row['stop'] + 10):
                if elem in occupied_neg or occupied_pos:
                    occur_count += 1
                else:
                    continue

        if occur_count > 0:
            continue
        elif occur_count == 0:
            count += 1



    print(count)

