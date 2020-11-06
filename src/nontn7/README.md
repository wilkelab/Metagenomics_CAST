# Non-Tn7 CRISPR-transposons

`main.sh` Runs the entire pipeline.  

`dedup.py` Deduplicates operons
`fix-paths.py` Updates paths to those on our local cluster
`identify-catalytic-residues.py` Finds residues that are canonically nuclease-active
`load-nuclease-dead.py` Filters operons with nuclease-dead effector proteins
`make-effector-fasta.py` Creates a FASTA file containing every instance of a given protein
`reblast.py` Runs BLAST on a random subset of operons using given protein and nucleotide databases. 
`rules-all.py` Filters out operons that cannot possibly contain non-Tn7 CRISPR-transposons.
`rules-cas12.py` Selects operons with cas12.
`rules-cas13.py` Selects operons with cas13.
`rules-cas9.py` Selects operons with cas9.
`rules-class1.py` Selects operons with Class 1 systems.
`rules-nocas1-2.py` Filters out operons with cas1 and/or cas2.
`rules-nocas3-10.py` Filters out operons with cas3 and/or cas10.
`self-targeting.py` Finds operons with CRISPR arrays that target a location in its contig.
`size-select.py` Filters operons by the size of a given protein.
`summarize.py` Provides an overview of the results of re-BLASTing operons.
