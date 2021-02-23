# We've re-BLASTed many systems, but only are manually selected for further review. In order to isolate the
# gene diagrams of those systems, this script reads the list that defines them and copies those images to
# another directory.

import sys
import os
import shutil


interesting_file = sys.argv[1]
plot_directory = sys.argv[2]
output_dir = sys.argv[3]

pngs = []
for group in 'cas9', 'cas12', 'cas13', 'class1':
    directory = os.path.join(plot_directory, group)
    group_pngs = [os.path.join(directory, path) for path in os.listdir(directory) if path.endswith('.png')]
    pngs.extend(group_pngs)

with open(interesting_file) as f:
    for line in f:
        directory, cluster_number = line.strip().split()
        prefix = f"{directory}/cluster{cluster_number:>03}-"
        output_group_dir = os.path.join(output_dir, directory)
        os.makedirs(output_group_dir, exist_ok=True)
        for png in pngs:
            if prefix in png:
                shutil.copy(png, output_group_dir)
