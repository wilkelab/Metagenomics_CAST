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

clusters = set()
found = set()

with open(interesting_file) as f:
    for line in f:
        directory, cluster_number = line.strip().split()
        clusters.add((directory, cluster_number))
        prefix = f"{directory}/cluster{cluster_number:>03}-"
        for png in pngs:
            if prefix in png:
                shutil.copy(png, os.path.join(output_dir, directory))
                found.add((directory, cluster_number))

print(clusters-found)
print(len(clusters))
print(len(clusters-found))
