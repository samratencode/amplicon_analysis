# amplicon_analysis
This repository contains R & shell scripts used for the analysis of amplicon based microbiome data (16S & ITS sequence).


Author: Samrat Ghosh (samratencode) Contact: samrat.ghosh2010@gmail.com


### sampling


$ ./amplicon_qiime2_sampling_wf.sh


### import database

$./amplicon_import_silva_db_wf.sh

$./amplicon_import_unite_db_wf.sh

### run QIIME2

$./amplicon_qiime2_bacteria_wf.sh

$./amplicon_qiime2_fungi_wf.sh


### export QIIME2 output
´´´
$./amplicon_exporting_wf.sh
´´´
### run R scripts
Open this script in Rstudio and run step by step:

amplicon_analysis_wf.R


