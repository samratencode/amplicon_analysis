#qiime tools import --type 'FeatureData[Sequence]' --input-path ./rep_set/rep_set_16S_only/99/silva_132_99_16S.fna --output-path silva_ver132_99_seqs_bac.qza
#qiime tools import --type 'FeatureData[Taxonomy]' --input-format HeaderlessTSVTaxonomyFormat --input-path ./taxonomy/16S_only/99/majority_taxonomy_7_levels.txt --output-path silva_ver132_99_tax_bac.qza
qiime feature-classifier fit-classifier-naive-bayes --i-reference-reads ref-seqs.qza  --i-reference-taxonomy ref-tax.qza --o-classifier silva_ver138_1_99_classifier_bac.qza
