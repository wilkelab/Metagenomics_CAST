# Type V Tn7 CAST Output

  - `canonical-pams.txt` PAM sequences from bona fide self-targeting systems
  - `cas12k-tns-minced-array.csv.gz` Putative Type V CASTs from the initial `gene_finder` run.
  - `other-pams.txt` PAM sequences from systems where we suspect the self-target annotation is an artifact and cannot be trusted
  - `reblasted-cas12k-with-sts-rules-deduped.csv.gz` Putative Type V CASTs after re-BLASTing with Cas12k as the seed.
  - `self-target-blastn-nt.tsv` Results of running blastn with the NCBI nt database on target sequences.
  - `self-target-blastx-nr.tsv` Results of running blastp with the NCBI nr database on target sequences.
  - `self-target-context-seqs.fa` DNA sequences 50 bp up/downstream of each target.
  - `unique-filenames.txt` Paths to raw sequence files for each CAST.

### TSV formats

`self-target-blastn-nt.tsv` and `self-target-blastx-nr.tsv` are created by `blastn-2.10` and `blastx-2.10`, respectively, with the `outfmt` parameter set to `6 qseqid sseqid evalue qseq sseq pident stitle`. These values mean:

| value | meaning |
| ----- | ------- | 
| qseqid | `<contig accession>_<contig_start>..<contig end>_<gene start>..<gene_end>` |
| sseqid | accession of nearest BLAST result |
| qseq | aligned part of query sequence |
| sseq | aligned part of BLAST result sequence |
| pident | percentage of identical matches |
| stitle | title of BLAST hit | 
