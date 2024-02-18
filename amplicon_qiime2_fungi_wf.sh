###DataImport
#qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-path manifest_fungi.csv --output-path paired_end_fungi_demux.qza --input-format PairedEndFastqManifestPhred33


##PrimersTrim 
#qiime cutadapt trim-paired --i-demultiplexed-sequences paired_end_fungi_demux.qza --p-front-f TAGAGGAAGTAAAAGTCGTAA --p-front-r TTCAAAGATTCGATGATTCA --p-error-rate 0 --o-trimmed-sequences paired_end_fungi_demux_trimmed.qza

##Check quality plot
#qiime demux summarize --i-data paired_end_fungi_demux_trimmed.qza --o-visualization paired_end_fungi_demux_trimmed.qzv


###Dereplication-DADA2
qiime dada2 denoise-paired --i-demultiplexed-seqs paired_end_fungi_demux_trimmed.qza --p-trim-left-f 0 --p-trim-left-r 0 --p-trunc-len-f 200 --p-trunc-len-r 200 --p-n-threads 30 --o-table 18S_table.qza --o-representative-sequences 18S_rep_seqs.qza --o-denoising-stats 18S_denoising_stats.qza --verbose


###Fungi-UNITE-OTU-referene-Phylogeny
qiime feature-classifier classify-sklearn --i-classifier /home/xxxx/SOFTWARE/UNITE9/fungi_unite_insdc_classifier_latest.qza --i-reads 18S_rep_seqs.qza --o-classification 18S_classified_rep_seqs.qza
qiime phylogeny align-to-tree-mafft-fasttree --i-sequences 18S_rep_seqs.qza --output-dir phylogenetic_tree_fun --p-n-threads 30 --verbose &> phylogenetic_fun_tree_generation.log 

