# Non-Tn7 CRISPR transposon data

  - `reblast` contains `gene_finder` data for interesting systems that have been BLASTed with the Trembl/Swissprot/tRNA databases. 
  - `systems-with-crispr-arrays` datasets where CRISPR arrays were required, among other rules
  - `systems-with-or-without-crispr-arrays` datasets where CRISPR arrays were optional, among other rules
  - `interesting-candidates.csv` manually-selected clusters from the CRISPR-array-optional datasets that were prioritized for re-BLASTing

### The following files are currently excluded from git due to their large size but will be regenerated when running the pipeline from scratch:

  - `minimal_passing_systems.csv.gz` all systems that met our minimal criteria and optionally contained CRISPR arrays
  - `minimal_passing_systems_with_crispr_arrays.csv.gz` all systems that met our minimal criteria and definitely contained CRISPR arrays
