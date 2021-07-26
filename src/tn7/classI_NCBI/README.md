# tn7 scripts to identify Tn7 CRISPR-transposons in NCBI genebank database#

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

  - **Protein BLAST database:** cas, tns, attsite fasat in "data/Tn7" folder, run `makeblastdb` (version 2.10.1+) with the default options
  - **Nucleotide database:** ffs fasta in "data/Tn7" folder, run `makeblastdb` (version 2.10.1+) with the default options


### Scripts

  - `typeI-F_ffs.py` Selects CAST IF systems with ffs attachment site.
  - `typeI-F_yciA.py` Selects CAST IF systems with yciA attachment site.
  - `typeI-F_rsmJ.py` Selects CAST IF systems with rsmJ attachment site.
  - `typeI-F_guaC.py` Selects CAST IF systems with guaC attachment site.
  - `gene_finder_att_site.py` Runs BLAST on contig using given protein and nucleotide databases.  
  - `png_generator.py` Plots putative systems with annotations from our custom database.
  - `tools/`  
    `colors.py` Contains the color palette for system diagrams 

