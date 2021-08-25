# Non-Tn7 CRISPR transposon data

  - `cas12-rpn-candidates.csv` Cas12-Rpn systems described in Figure 6
  - `interesting-candidates.txt` manually-selected clusters from the CRISPR-array-optional datasets that were prioritized for re-BLASTing
  - `reblast` contains `gene_finder` data for interesting systems that have been BLASTed with the Trembl/Swissprot/tRNA databases. 
  - `systems-with-crispr-arrays` datasets where CRISPR arrays were required, among other rules
  - `systems-with-or-without-crispr-arrays` datasets where CRISPR arrays were optional, among other rules

### The following files are currently excluded from git due to their large size but will be regenerated when running the pipeline from scratch:

  - `minimal_passing_systems.csv.gz` all systems that met our minimal criteria and optionally contained CRISPR arrays
  - `minimal_passing_systems_with_crispr_arrays.csv.gz` all systems that met our minimal criteria and definitely contained CRISPR arrays
