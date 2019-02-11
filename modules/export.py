import csv
import pandas as pd

def exportpandas(export_df, output, projectname, status):
    export_df = export_df.drop_duplicates(subset=["start", "stop", "strand", "attributes"])
    export_df = export_df.reset_index(drop=True)
    export_df = export_df.sort_values(["start", "stop"], axis=0).reset_index(drop=True)
    export_df.to_csv(str(output) + str(projectname) + str(status), sep='\t', index=False, header=False)

def export(export_list, output, projectname, status):
    with open(str(output) + str(projectname) + str(status), 'w') as output_file:
        csv_out = csv.writer(output_file, delimiter='\t')
        csv_out.writerow(['seqid', 'source', 'type', 'start', 'stop', 'score', 'strand', 'phase', 'attributes'])
        csv_out.writerows(export_list)