# Non-Tn7 CRISPR transposon data

  - `reblast` contains `gene_finder` data for interesting systems that have been BLASTed with the Trembl/Swissprot/tRNA databases. 
  - `systems-with-crispr-arrays` datasets where CRISPR arrays were required, among other rules
  - `systems-with-or-without-crispr-arrays` datasets where CRISPR arrays were optional, among other rules
  - `interesting-candidates.csv` manually-selected clusters from the CRISPR-array-optional datasets that were prioritized for re-BLASTing

### The following files are currently excluded from git due to their large size:

  - `minimal-composite-with-arrays.csv.gz` potential composite transposons with CRISPR arrays
  - `minimal-composite.csv.gz` potential composite transposons with or without CRISPR arrays
  - `minimal-tn3-with-arrays.csv.gz` systems with Tn3-family transposases and CRISPR arrays
  - `minimal-tn3.csv.gz` systems with Tn3-family transposases
  - `minimal_passing_systems.csv.gz` all systems that met our minimal criteria and optionally contained CRISPR arrays
  - `minimal_passing_systems_with_crispr_arrays.csv.gz` all systems that met our minimal criteria and definitely contained CRISPR arrays
