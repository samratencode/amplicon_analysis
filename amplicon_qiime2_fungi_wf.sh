###DataImport
qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-path manifest_fungi.csv --output-path paired_end_fungi_demux.qza --input-format PairedEndFastqManifestPhred33


###Dereplication-DADA2
qiime dada2 denoise-paired --i-demultiplexed-seqs paired_end_fungi_demux.qza --p-trim-left-f 0 --p-trim-left-r 0 --p-trunc-len-f 0 --p-trunc-len-r 0 --p-n-threads 30 --o-table 18S_table.qza --o-representative-sequences 18S_rep_seqs.qza --o-denoising-stats 18S_denoising_stats.qza --verbose


###Fungi-UNITE-OTU-referene-Phylogeny
qiime feature-classifier classify-sklearn --i-classifier /home/xxxx/SOFTWARES/UNITE_9_QIIME_fungi/unite-ver9-dynamic-classifier-fungi_29.11.2022.qza --i-reads 18S_rep_seqs.qza --o-classification 18S_classified_rep_seqs.qza
qiime phylogeny align-to-tree-mafft-fasttree --i-sequences 18S_rep_seqs.qza --output-dir phylogenetic_tree_fun --p-n-threads 30 --verbose &> phylogenetic_fun_tree_generation.log 

