# Non-Tn7 CASTs

A pipeline for discovering non-Tn7 CASTs.

### Software Dependencies

  - `Python 3.6+`  
  - `Opfi` (https://github.com/alexismhill3/Opfi) and its dependencies  
  - `GNU parallel verion: 20161222` [https://www.gnu.org/software/parallel/](https://www.gnu.org/software/parallel/)  
  - `MAFFT v7.310 (2017/Mar/17)` [https://mafft.cbrc.jp/alignment/software/](https://mafft.cbrc.jp/alignment/software/)  
  - `MMseqs2 version: 45c4de7f1daefa06b45688195305eadedaea4d97` [https://github.com/soedinglab/MMseqs2](https://github.com/soedinglab/MMseqs2)  
  - `NCBI Blast 2.10+` [https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)  
  - `ripgrep` [https://github.com/BurntSushi/ripgrep](https://github.com/BurntSushi/ripgrep)  
  - `fd` [https://github.com/sharkdp/fd](https://github.com/sharkdp/fd)  
  - `seqkit 0.14.0` [https://bioinf.shenwei.me/seqkit/](https://bioinf.shenwei.me/seqkit/)

### Data Dependencies

  - **Protein BLAST database:** combine the TrEMBL and Swissprot databases ([https://www.uniprot.org/downloads](https://www.uniprot.org/downloads)) into a single protein database by concatenating their FASTA files and run `makeblastdb` (version 2.10.1+) with the default options.
  - **Nucleotide database:** download the bacterial and archaeal tRNA sequences from ([http://lowelab.ucsc.edu/GtRNAdb/download.html](http://lowelab.ucsc.edu/GtRNAdb/download.html)) and run `makeblastdb` (version 2.10.1+) with the default options.

### Inputs and Outputs

Data directories are defined at the top of `main.sh`. The `$INPUT` and `$STORE` directories will need to be changed if you want to run the pipeline on your own system. `$INPUT` is the location of the `gene_finder` CSVs, and `$STORE` should be a filesystem with the capacity to store all the inputs. `$OUTPUT` and `$DATA` refer to directories within this repo.

The pipeline requires a set of CSVs produced by `gene_finder` as its input. It will search through those putative systems and produce these outputs:
  - `$crispr.fully-analyzed`: `gene_finder`-formatted CSVs of each candidate system, containing details on all proteins and CRISPR arrays, self-targeting spacers and inverted repeats. `$crispr` here is one of: `cas9`, `cas12`, `cas13` or `class1`.
  - `reblast`: the pipeline will take the systems in the `fully-analyzed` directories and run them through `gene_finder` with the TrEMBL, Swissprot and tRNA databases. The results are stored in the `reblast` directory, separated by type.
  - `plots`: for each candidate system, a diagram will be generated containing the gene annotations from the curated database (on top) and the TrEMBL/Swissprot/tRNA annotations (on bottom).

Optionally, after running the pipeline and manually picking systems that merit further examination, create a file named `plausible-system-ids.txt` in the output directory with the first two values of the systems CSV entries (i.e., formatted like: `contig_accession,start_coordinate..end_coordinate`). Running `./main.sh` again will now generate `plausible.csv.gz` which contains detailed information on the genes and other features of each candidate system.

### Running the Pipeline

Issue the command: `./main.sh`  

Tests can be run with `python -m pytest` in this directory.  

### Scripts

  - `cas12k-nuclease-dead-control.sh` Control to ensure that Cas12k genes from NCBI are all assigned as nuclease-dead
  - `count_systems.py` Counts the number of unique putative systems in a CSV 
  - `dedup.py` Deduplicates systems  
  - `extract_candidates.py` Serializes manually-selected candidate systems
  - `find-cas12-sts.py` Looks for Cas12-containing systems and tries to identify their self-targeting spacers, if any exist
  - `find_inverted_repeats.py` Finds inverted repeats in contigs that contain putative systems
  - `find_nuclease_dead_clusters.py` Determine which Class 2 systems have inactivated catalytic residues
  - `fix-paths.py` Updates paths to those on our local cluster  
  - `get_clustered_systems.py` Parses output from mmseqs2 to help cluster candidate systems
  - `identify-catalytic-residues.py` Finds residues that are canonically nuclease-active  
  - `load-nuclease-dead.py` Filters systems with nuclease-dead effector proteins  
  - `make-effector-fasta.py` Creates a FASTA file containing every instance of a given protein  
  - `paper-stats.py` Performs calculations for the manuscript. Used in `paper-values.sh`  
  - `paper-values.sh` Performs calculations for the manuscript  
  - `plot_reblast.py` Plots putative systems with annotations from our custom database, alongside annotations from TrEMBL/Swissprot and a tRNA database  
  - `reblast.py` Runs BLAST on a random subset of systems using given protein and nucleotide databases  
  - `reblast_interesting_candidates.py` Selectively runs BLAST on a given list of candidate systems 
  - `rules-array.py` Selects systems with CRISPR arrays
  - `rules-cas9.py` Selects systems with cas9  
  - `rules-cas12.py` Selects systems with cas12  
  - `rules-cas13.py` Selects systems with cas13  
  - `rules-class1.py` Selects systems with Class 1 CRISPR systems  
  - `rules-minimal.py` Selects systems with one transposase, one cas gene, and no Tn7 genes
  - `rules-nocas1-2.py` Filters out systems with cas1 and/or cas2  
  - `rules-nocas3-10.py` Filters out systems with cas3 and/or cas10  
  - `self_targeting.py` Finds systems with CRISPR arrays that target a location in its contig  
  - `simple-reblast.py` Re-BLASTs systems with given nucleotide and protein databases
  - `size_select.py` Filters systems by the size of a given protein  
  - `validate-plausible-candidates.py` Prints detailed information about a manually-curated list of putative systems  
  - `tools/`  
    `colors.py` Contains the color palette for system diagrams  
    `filters.py` Contains a universal filter 
