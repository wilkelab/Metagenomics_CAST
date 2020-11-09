from operon_analyzer import genes, visualize
import flabpal
import gzip
from Bio import SeqIO
import subprocess


FEATURE_COLORS = {'transposase': flabpal.blue,
                  'cas9': flabpal.red,
                  'cas12': flabpal.green,
                  'cas5': flabpal.pink,
                  'cas6': flabpal.brown,
                  'cas7': flabpal.purple,
                  'cas8': flabpal.yellow,
                  'CRISPR array': flabpal.orange,
                  '': flabpal.gray}


def explore_operon(operon: genes.Operon):
    fig = visualize.create_operon_figure(operon, plot_ignored=False, feature_colors=FEATURE_COLORS)
    if fig is None:
        raise ValueError("Could not create image")
    visualize.save_operon_figure(fig, '/tmp/explore-operon.png')
    with gzip.open(operon.contig_filename, 'rt') as f:
        records = SeqIO.parse(f, 'fasta')
        for record in records:
            if record.id == operon.contig:
                with open("/tmp/explore-contig.fasta", "w") as f:
                    SeqIO.write(record, f, 'fasta')
                break

        command = ["pilercr", "-in", '/tmp/explore-contig.fasta',
                   "-out", '/tmp/explore-pilerout.txt',
                   "-minarray", "2",
                   "-quiet",
                   "-noinfo"]
        subprocess.call(command)
