# Determines whether every candidate system in a group of candidates is nuclease dead,
# and if so, prints the filename of the gzipped CSV file to stdout


import os
import sys

from operon_analyzer import load

nuclease_dead = set(load.load_gzipped_operons(sys.argv[1]))
cluster_directory = sys.argv[2]


for cluster_file in os.listdir(cluster_directory):
    if not cluster_file.endswith('csv.gz'):
        continue
    cluster = tuple(load.load_gzipped_operons(os.path.join(cluster_directory, cluster_file)))
    dead_count = sum([1 for operon in cluster if operon in nuclease_dead])
    if dead_count == len(cluster):
        print(cluster_file)
