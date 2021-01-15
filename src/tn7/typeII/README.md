A pipeline for discovering Tn7 CRISPR transposons.
### Software Dependencies
  - `Python 3.6+`  
  - `Opfi` (https://github.com/alexismhill3/Opfi) and its dependencies  
  - `GNU parallel verion: 20161222` [https://www.gnu.org/software/parallel/](https://www.gnu.org/software/parallel/)  
  - `NCBI Blast 2.10+` [https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)  
  
### Inputs and Outputs

Data directories are defined at the top of `main.sh`. The `$INPUT` directories will need to be changed if you want to run the pipeline on your own system. `$INPUT` is the location of the `gene_finder` CSVs.

The pipeline requires a set of CSVs produced by `gene_finder` as its input. It will search through those putative operons and generate figure for these outputs.

### Running the Pipeline

Issue the command: `./main.sh`  
