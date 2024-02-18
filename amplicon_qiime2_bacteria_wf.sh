###DataImport 
#qiime tools import --type 'SampleData[SequencesWithQuality]' --input-path manifest_bacteria.csv --output-path single_end_bacteria_demux.qza  --input-format SingleEndFastqManifestPhred33V2
#qiime tools import --type SampleData[JoinedSequencesWithQuality] --input-path manifest_bacteria.csv  --output-path joined_end_bacteria_demux.qza  --input-format SingleEndFastqManifestPhred33V2
qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-path manifest_bacteria.csv --output-path paired_end_bacteria_demux.qza --input-format PairedEndFastqManifestPhred33

##PrimersTrim 
qiime cutadapt trim-paired --i-demultiplexed-sequences paired_end_bacteria_demux.qza --p-front-f AACMGGATTAGATACCCKG --p-front-r AGGGTTGCGCTCGTTRC --p-error-rate 0 --o-trimmed-sequences paired_end_bacteria_demux_trimmed.qza

##Check quality plot
qiime demux summarize --i-data paired_end_bacteria_demux_trimmed.qza --o-visualization paired_end_bacteria_demux_trimmed.qzv

###Dereplication-DADA2
#qiime dada2 denoise-single --i-demultiplexed-seqs single_end_bacteria_demux.qza --p-trim-left 0 --p-trunc-len 0 --p-n-threads 30 --o-table 16S_table.qza --o-representative-sequences 16S_rep_seqs.qza  --o-denoising-stats 16S_denoising_stats.qza --verbose
qiime dada2 denoise-paired --i-demultiplexed-seqs paired_end_bacteria_demux_trimmed.qza --p-trim-left-f 0 --p-trim-left-r 0 --p-trunc-len-f 190 --p-trunc-len-r 190 --p-n-threads 30 --o-table 16S_table.qza --o-representative-sequences 16S_rep_seqs.qza --o-denoising-stats 16S_denoising_stats.qza --verbose

###Bacteria-SILVA-OTU-referene-Phylogeny
qiime feature-classifier classify-sklearn --i-classifier /home/xxxx/SOFTWARE/SILVA138_1/silva_ver138_1_99_classifier_bac.qza --i-reads 16S_rep_seqs.qza --o-classification 16S_classified_rep_seqs.qza
qiime phylogeny align-to-tree-mafft-fasttree --i-sequences 16S_rep_seqs.qza --output-dir phylogenetic_tree_bac --p-n-threads 30 --verbose &> phylogenetic_bac_tree_generation.log 
