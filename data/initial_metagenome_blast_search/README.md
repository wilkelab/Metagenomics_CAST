This directory contains the following:

`databases/` Blast databases used for the intial metagenomic search for putative CRISPR-Transposons

`cas*.fasta` Representative amino acid sequences for each cas protein type/subtype

`transposase.fasta` Representative amino acid sequences for all available transposases, including the heteromeric Tn7 transposase proteins TnsA and TnsB

`tn7-accessory.fasta` Representative amino acid sequences for Tn7 accessory proteins, i.e TnsC, TnsD/TniQ, and TnsE

**Blast database construction**

`databases/tn7-accessory` and `databases/transposase` were built using the corresponding fasta files as inputs. For `databases/cas_all`, all cas protein sequence files were concatenated before piping to `makeblastdb`.

The specific `makeblastdb` command is:

`cat transposase.fasta | makeblastdb -title "transposase" -out transposase -dbtype prot -hash_index`

**Protein sequence collection/curation**

**Transposases:**

All available bacterial and archaeal transposase sequences were downloaded from UniRef50. All transposases associated with transposons listed in the Transposon Registry (https://transposon.lstmed.ac.uk/tn-registry) were downloaded from NCBI. Finally, 100 transposases associated with each of the major families of insertion sequences were downloaded from NCBI (sorted by relevance).

Separately, sequences for the Tn7 heteromeric transposase proteins TnsA and TnsB were downloaded from UniRef50. 

All transposase sequences were combined and then deduplicated with the `MMseqs2` (https://github.com/soedinglab/MMseqs2) easy-cluster workflow, using default parameters. Representative sequences from each cluster were used for the final `transposase.fasta` collection and corresponding blast database.

**Accessory Tn7 proteins:**

Sequences for the Tn7 accessory proteins TnsC, TnsD/TniQ, and TnsE were downloaded from UniRef50 and combined in `tn7-accessory.fasta`. 

**Cas proteins:**

Sequences for Cas proteins 1-11 were downloaded from UniRef50. For Cas12, Cas13, and Casphi sequence collection, a combination of UniProt (UniRef), NCBI (protein), and primary literature sources was necessary to ensure that newly discovered variants be included. Cas12 and Cas13 sequences were deduplicated using `CD-HIT` (http://weizhongli-lab.org/cd-hit/) with a 50% sequence identity threshold and 80% alignment overlap. 

To ensure that the representative Cas protein sequences be as comprehensive as possible, we compiled (for each major Cas type, i.e 1-13) a list of contemporary and historic protein names to search for (summarized in the following tables). 

Cas1-11:
|Common Name|Alternate Name(s) / Query|Source|
|---|---|---|
|Cas1|cas1|UniRef50|
|Cas2|cas2|UniRef50|
|Cas3|cas3|UniRef50|
|Cas4|cas4|UniRef50|
|Cas5|pbprb1992 OR csf3 OR csc1 OR csy2 OR cmr3 OR csx10 OR csm4 OR cas5|UniRef50|
|Cas6|cas6 OR cas6f OR cas6e|UniRef50|
|Cas7|cas7 OR csc2 OR csf2 OR cmr1 OR cmr4 OR csm5 OR cmr6 OR csm3 OR pbprb1993 OR csy3|UniRef50|
|Cas8|cas8 OR cse1 OR csy1 OR cas8b OR csh1 OR cas8c OR csd1 OR cst1 OR cmx1 OR csf1 OR cas8a2 OR csa4 OR csx9|UniRef50|
|Cas9|cas9|UniRef50|
|Cas10|cas10 OR cas10d OR csx11 OR csm1 OR cmr2|UniRef50|
|Cas11|cas11 OR csa5 OR cmr5 OR cse2 OR csm2|UniRef50|

Cas12:
|Subtype|Source|Alternate Name(s) / Query|DOI|
|---|---|---|---|
|Cas12 (all)|UniRef50|cas12||
|Cas12a|UniRef50|cas12a OR cpf1||
|Cas12b|UniRef50|cas12b OR c2c1||
|Cas12c|UniRef50|cas12c||
|Cas12d|NCBI|cas12d OR casY||
|Cas12e|NCBI|cas12e OR casX||
|Cas12f|UniRef50|cas12f OR cas14||
|Cas12g|literature||https://doi.org/10.1126/science.aav7271|
|Cas12h|literature||https://doi.org/10.1126/science.aav7271|
|Cas12i|NCBI|cas12i||
|Cas12k|NCBI|cas12k OR c2c5||
|Cas12V-U1|literature||https://doi.org/10.1038/nrmicro.2016.184|
|Cas12V-U2|literature||https://doi.org/10.1038/nrmicro.2016.184|
|Cas12V-U3|literature||https://doi.org/10.1038/nrmicro.2016.184|
|Cas12V-U4|literature||https://doi.org/10.1038/nrmicro.2016.184|
|Cas12V-U5|literature||https://doi.org/10.1038/nrmicro.2016.184|

Cas13:
|Common Name|Source|Alternate Name(s) / Query|
|---|---|---|
|Cas13|NCBI|cas13 OR cas13a OR cas13b OR cas13c OR cas13d|

Casphi:
|Common Name|Source|DOI|
|---|---|---|
|Casphi|literature|https://doi.org/10.1126/science.abb1400|
