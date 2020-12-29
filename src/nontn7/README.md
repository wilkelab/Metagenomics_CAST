# Non-Tn7 CRISPR-transposons

`main.sh` Runs the entire pipeline.  
Tests can be run with `python -m pytest` in this directory.

`find_nuclease_dead_clusters.py` Determine which Class 2 systems have inactivated catalytic residues
`count_operons.py` Counts the number of unique putative operons in a CSV 
`count_passing_operons.sh` Counts the number of systems that passed our initial filters
`dedup.py` Deduplicates operons  
`extract_candidates.py` Serializes manually-selected candidate systems
`find_inverted_repeats.py` Finds inverted repeats in contigs that contain putative operons
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
`summarize.py` Provides an overview of the results of re-BLASTing operons  
`validate-plausible-candidates.py` Prints detailed information about a manually-curated list of putative operons
