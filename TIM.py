#!/usr/bin/env python
import argparse
from modules import *
import time

print("\nHello, I am T.I.M. Let's get started!")
total_start_time = time.time()
#####Define Input Arguments#####
parser = argparse.ArgumentParser(description='Input Files')
parser.add_argument('--assembly', metavar='A', default=["./test_data/Assembly.gff3"], required=False, type=str, nargs=1,
                   help='The mapped assembly in gtf format (filepath from pwd)')

parser.add_argument('--genome', metavar='G', default=["./test_data/Genome.gff3"], required=False, type=str, nargs=1,
                   help='The annotated reference genome in gff3 format (filepath from pwd)')

parser.add_argument('--distance', metavar='D', default=["10"], type=int, nargs=1,
                   help='Distance of kept feature to annotated features')

parser.add_argument('--output_path', metavar='O', default=["./"], type=str, nargs=1,
                   help='Path of the output files')

parser.add_argument('--project_prefix', metavar='N', default=["TIM"], type=str, nargs=1,
                   help='Project name, all output files will be prefixed with this')

parser.add_argument('--ss', metavar='N', default=["0"], type=int, nargs=1,
                   help='Only consider intersects, if they occur on the same strand as the respective feature? 1= Yes 0=No (Default:0)')

parser.add_argument('--merge', metavar='M', default=["0"], type=int, nargs=1,
                   help='Should non-intersecting features be merged if they overlap? 0= Yes 1=No (Default:0)')




#####Load Input Arguments into variables#####
args = parser.parse_args()
assembly_path=args.assembly[0]
genome_path=args.genome[0]
distance=int(args.distance[0])+1
output=str(args.output_path[0])
projectname=str(args.project_prefix[0])
ss=int(args.ss[0])
mergeyn=int(args.merge[0])

#####Column Names of Input FIles#####
col_names = ['seqid', 'source', 'type', 'start', 'stop', 'score', 'strand', 'phase', 'attributes']


#####Load Assembly into Object#####
assey = Assembly()
assey = assey.load_assembly(assembly_path, col_names)

#####Load Genome into Object#####
geno = Genome()
geno = geno.load_genome(genome_path, col_names)


#####Calculate all genomic positions occupied by annotations in the Genome for each strand#####
strand_pos, strand_neg = geno.calc_intersect()

###Look for intersections between the assembly and the genome annotation#####
inters = IntersTable(strand_pos, strand_neg)
if ss == 0:
    inters.intersect_ns(assey, distance)
elif ss == 1:
    inters.intersect_ss(assey, distance)
else:
    print("Incompatible parameter given!")

print("\nNumber of non-intersecting features found: " + str(len(inters.non_intersecting)))
print("\nIntersection Analysis concluded within "+ str((time.time() - total_start_time)) + " seconds. \nStarting with feature merging...")

#####Merge remaining features#####
if mergeyn == 0:
    inters.merge()
else:
    pass

print("\nNumber of newly merged features: " + str(len(inters.merged_final)))
print("\nMerging concluded within "+ str((time.time() - total_start_time)) + " seconds. \nStarting with data export...")


#####Exporting all non-intersecting features premerging#####
export(inters.non_intersecting, output, projectname,"non-intersecting_premerge.tsv")

#####Exporting all intersecting/discarded features premerging#####
export(inters.intersecting, output, projectname,"intersecting_discarded.tsv")

#####Outputing newly created features from merges#####
export(inters.merged_final, output, projectname,"newly_merged_features.tsv")

print("\nTotal Analysis concluded within "+ str((time.time() - total_start_time)) + " seconds. Good Bye!")


