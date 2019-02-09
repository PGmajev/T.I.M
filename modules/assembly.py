import pandas as pd
import sys
import csv
import time

class assembly():
    def __init__(self):
        self.df = pd.DataFrame()

    def load_assembly(self, assembly_path):
        start_time = time.time()
        col_names = ['seqid', 'source', 'type', 'start', 'stop', 'score', 'strand', 'phase', 'attributes']
        df = pd.read_csv(assembly_path, sep='\t', comment='#', low_memory=False, header=None, names=col_names)
        df = df.sort_values(["strand", "start", "stop"], axis=0).reset_index(drop=True)
        df = df.drop_duplicates(subset=["start", "stop", "strand"])
        df = df.reset_index(drop=True)
        df.loc[df['type'] != "exon", 'type'] = "exon"
        self.df = df
        print(self.df.index)
        print("Finished reading in Assembly in " + str((time.time() - start_time) * 1000) + " milliseconds")
        return self