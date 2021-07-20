import sys
import os
from filters import fs
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

from operon_analyzer import load


def get_feature(filtered_features, name):
    features = [feature for feature in filtered_features if name in feature.name]
    if len(features) != 1:
        return None
    return features[0]

output_dir = sys.argv[1]


for operon in load.load_operons(sys.stdin):
    fs.evaluate(operon)
    features = list(sorted(set(operon), key=lambda feature: feature.start))
    filtered_features = []
    for feature in features:
        for keyword in ('Repeat Spacer', 'Target from', 'Target for', 'Self-targeting spacer'):
            if keyword in feature.name:
                break
        else:
            filtered_features.append(feature)

    topo = get_feature(filtered_features, 'type IA DNA topoisomerase')
    viroplasmin = get_feature(filtered_features, 'viroplasmin')
    if topo is None and viroplasmin is None:
        continue

    atp_helicase = get_feature(filtered_features, 'ATP-dependent helicase')
    deadbox_helicase = get_feature(filtered_features, 'DEAD/DEAH box helicase')
    if atp_helicase is None and deadbox_helicase is None:
        continue

    upfeature = topo if topo is not None else viroplasmin
    downfeature = deadbox_helicase if deadbox_helicase is not None else atp_helicase
    upindex = filtered_features.index(upfeature)
    downindex = filtered_features.index(downfeature)
    if upindex < downindex:
        upfeature, downfeature = downfeature, upfeature
        upindex, downindex = downindex, upindex

    cas12 = get_feature(filtered_features, 'Cpf1')
    assert cas12 is not None

    cas12_index = filtered_features.index(cas12)
    dnapol_index = filtered_features.index(cas12)

    assert downindex < cas12_index < upindex
    assert downindex < dnapol_index < upindex

    good_features = filtered_features[max(0, downindex - 2): upindex + 3]
    assert len(good_features) > 4
    for feature in cas12, upfeature, downfeature:
        assert feature in good_features

    dnapols = [feature for feature in good_features if 'DNA polymerase III' in feature.name]
    if len(dnapols) == 0:
        continue

    coords = []
    for feature in good_features:
        coords.append(feature.start)
        coords.append(feature.end)
    start, end = min(coords), max(coords)

    sequence = load.load_sequence(operon)
    nucleotide = sequence[start - 200: end + 200]
    fasta_filename = os.path.join(output_dir, 'nucleotide', f"{operon.contig}_{operon.start}_{operon.end}.fasta")
    if os.path.exists(fasta_filename):
        continue

    seqrecord = SeqRecord(sequence, f"{operon.contig}_{operon.start}_{operon.end}", name="", description="")
    with open(fasta_filename, "w") as f:
        SeqIO.write(seqrecord, f, 'fasta')

    protein_records = []
    for feature in good_features:
        seqrecord = SeqRecord(seq=Seq(feature.sequence), id=feature.name.replace(" ", "_"), name=feature.name, description=feature.description)
        protein_records.append(seqrecord)

    protein_fasta_filename = os.path.join(output_dir, 'protein', f"{operon.contig}_{operon.start}_{operon.end}.fasta")
    if os.path.exists(protein_fasta_filename):
        continue
    with open(protein_fasta_filename, "w") as f:
        SeqIO.write(protein_records, f, 'fasta')
