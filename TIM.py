#!/usr/bin/env python
import argparse
from modules import *
import time

print("Hello, I am T.I.M, lets get started...")
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

assey = assembly()
geno = genome()

assey = assey.load_assembly(assembly_path)

geno = geno.load_genome(genome_path)
strand_pos, strand_neg = geno.calc_intersect()
print("Length of pos: " + str(len(strand_pos)))
print("Length of neg: " + str(len(strand_neg)))
max_size_tc = geno.max_size()
intersect(assey, strand_pos, strand_neg)


print("Analysis concluded within "+ str((time.time() - total_start_time)) + " seconds")
