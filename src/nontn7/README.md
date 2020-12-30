# Non-Tn7 CRISPR-transposons

A pipeline for discovering non-Tn7 CRISPR transposons.

### Software Dependencies

  - `Python 3.6+`  
  - `Opfi` (https://github.com/alexismhill3/Opfi) and its dependencies  
  - `GNU parallel verion: 20161222` [https://www.gnu.org/software/parallel/](https://www.gnu.org/software/parallel/)  
  - `MAFFT v7.310 (2017/Mar/17)` [https://mafft.cbrc.jp/alignment/software/](https://mafft.cbrc.jp/alignment/software/)  
  - `MMseqs2 version: 45c4de7f1daefa06b45688195305eadedaea4d97` [https://github.com/soedinglab/MMseqs2](https://github.com/soedinglab/MMseqs2)  
  - `NCBI Blast 2.10+` [https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)  
  - `ripgrep` [https://github.com/BurntSushi/ripgrep](https://github.com/BurntSushi/ripgrep)  
  - `fd` [https://github.com/sharkdp/fd](https://github.com/sharkdp/fd)  

### Data Dependencies

  - **Protein BLAST database:** combine the TrEMBL and Swissprot databases ([https://www.uniprot.org/downloads](https://www.uniprot.org/downloads)) into a single protein database by concatenating their FASTA files and run `makeblastdb` (version 2.10.1+) with the default options.
  - **Nucleotide database:** download the bacterial and archaeal tRNA sequences from ([http://lowelab.ucsc.edu/GtRNAdb/download.html](http://lowelab.ucsc.edu/GtRNAdb/download.html)) and run `makeblastdb` (version 2.10.1+) with the default options.

### Inputs and outputs

Data directories are defined at the top of `main.sh`. The `$INPUT` and `$STORE` directories will need to be changed if you want to run the pipeline on your own system. `$INPUT` is the location of the `gene_finder` CSVs, and `$STORE` should be a filesystem with the capacity to store all the inputs. `$OUTPUT` and `$DATA` refer to directories within this repo.

The pipeline requires a set of CSVs produced by `gene_finder` as its input. It will search through those putative operons and produce these outputs:
  - `$crispr.fully-analyzed`: `gene_finder`-formatted CSVs of each candidate operon, containing details on all proteins and CRISPR arrays, self-targeting spacers and inverted repeats. `$crispr` here is one of: `cas9`, `cas12`, `cas13` or `class1`.
  - `reblast`: the pipeline will take the operons in the `fully-analyzed` directories and run them through `gene_finder` with the TrEMBL, Swissprot and tRNA databases. The results are stored in the `reblast` directory, separated by type.
  - `plots`: for each candidate operon, a diagram will be generated containing the gene annotations from the curated database (on top) and the TrEMBL/Swissprot/tRNA annotations (on bottom).

Optionally, after running the pipeline and manually picking systems that merit further examination, create a file named `plausible-operon-ids.txt` in the output directory with the first two values of the operons CSV entries (i.e., formatted like: `contig_accession,start_coordinate..end_coordinate`). Running `./main.sh` again will now generate `plausible.csv.gz` which contains detailed information on the genes and other features of each candidate system.

### Running the pipeline

Issue the command: `./main.sh`  

Tests can be run with `python -m pytest` in this directory.  

### Software Dependencies

`Python 3.6+`
`Opfi` (https://github.com/alexismhill3/Opfi) and its dependencies
`GNU parallel verion: 20161222` [https://www.gnu.org/software/parallel/](https://www.gnu.org/software/parallel/)
`MAFFT v7.310 (2017/Mar/17)` [https://mafft.cbrc.jp/alignment/software/](https://mafft.cbrc.jp/alignment/software/)
`MMseqs2 version: 45c4de7f1daefa06b45688195305eadedaea4d97` [https://github.com/soedinglab/MMseqs2](https://github.com/soedinglab/MMseqs2)
`NCBI Blast 2.10+` [https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
`ripgrep` [https://github.com/BurntSushi/ripgrep](https://github.com/BurntSushi/ripgrep)
`fd` [https://github.com/sharkdp/fd](https://github.com/sharkdp/fd)

### Data Dependencies

Protein BLAST database: combine the TrEMBL and Swissprot databases [https://www.uniprot.org/downloads](https://www.uniprot.org/downloads) into a single protein database by concatenating their FASTA files and run `makeblastdb` (version 2.10.1+) with the default options.

Nucleotide database: download the bacterial and archaeal tRNA sequences from [http://lowelab.ucsc.edu/GtRNAdb/download.html](http://lowelab.ucsc.edu/GtRNAdb/download.html) and run `makeblastdb` (version 2.10.1+) with the default options.

### Inputs and outputs

Data directories are defined at the top of `main.sh`. The `$INPUT` and `$STORE` directories will need to be changed if you want to run the pipeline on your own system. `$INPUT` is the location of the `gene_finder` CSVs, and `$STORE` should be a filesystem with the capacity to store all the inputs. `$OUTPUT` and `$DATA` refer to directories within this repo.


The pipeline requires a set of CSVs produced by `gene_finder` as its input. It will search through those putative operons and produce these outputs:
  - `$crispr.fully-analyzed`: `gene_finder`-formatted CSVs of each candidate operon, containing details on all proteins and CRISPR arrays, self-targeting spacers and inverted repeats. $crispr here is one of: cas9, cas12, cas13 or class1.
  - `reblast`: the pipeline will take the operons in the `fully-analyzed` directories and run them through `gene_finder` with the TrEMBL, Swissprot and tRNA databases. The results are stored here, separated by type.
  - `plots`: for each candidate operon, a diagram will be generated containing the gene annotations from the curated database (on top) and the TrEMBL/Swissprot/tRNA annotations (on bottom).

The input CSVs have full paths to their FASTA file. This will need to be updated to a valid path on your system, and the value `KEEP_PATHS` set to `YES` in `main.sh`.

After running the pipeline and manually picking systems that merit further examination, create a file named `plausible-operon-ids.txt` in the output directory with the first two values of the operon's CSV entry (i.e., formatted like: `contig_accession,start_coordinate..end_coordinate`). Running `./main.sh` again will now generate `plausible.csv.gz` which contains detailed information on the genes and other features of each candidate system.

### Running the pipeline

Issue the command: `./main.sh` 

Tests can be run with `python -m pytest` in this directory.

### Scripts

`count_operons.py` Counts the number of unique putative operons in a CSV 
`count_passing_operons.sh` Counts the number of systems that passed our initial filters
`dedup.py` Deduplicates operons  
`extract_candidates.py` Serializes manually-selected candidate systems
`find_inverted_repeats.py` Finds inverted repeats in contigs that contain putative operons
`find_nuclease_dead_clusters.py` Determine which Class 2 systems have inactivated catalytic residues
`fix-paths.py` Updates paths to those on our local cluster  
`get_clustered_operons.py` Parses output from mmseqs2 to help cluster candidate systems
`identify-catalytic-residues.py` Finds residues that are canonically nuclease-active  
`load-nuclease-dead.py` Filters operons with nuclease-dead effector proteins  
`make-effector-fasta.py` Creates a FASTA file containing every instance of a given protein  
`paper-stats.py` Performs calculations for the manuscript. Used in `paper-values.sh`  
`paper-values.sh` Performs calculations for the manuscript  
`plot_reblast.py` Plots putative operons with annotations from our custom database, alongside annotations from TrEMBL/Swissprot and a tRNA database  
`reblast.py` Runs BLAST on a random subset of operons using given protein and nucleotide databases  
`rules-all.py` Filters out operons that cannot possibly contain non-Tn7 CRISPR-transposons  
`rules-cas12.py` Selects operons with cas12  
`rules-cas13.py` Selects operons with cas13  
`rules-cas9.py` Selects operons with cas9  
`rules-class1.py` Selects operons with Class 1 systems  
`rules-nocas1-2.py` Filters out operons with cas1 and/or cas2  
`rules-nocas3-10.py` Filters out operons with cas3 and/or cas10  
`self-targeting.py` Finds operons with CRISPR arrays that target a location in its contig  
`size-select.py` Filters operons by the size of a given protein  
`validate-plausible-candidates.py` Prints detailed information about a manually-curated list of putative operons  
`tools/`  
  `colors.py` Contains the color palette for operon diagrams  
  `filters.py` Contains a universal filter  
