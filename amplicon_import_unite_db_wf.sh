#awk '/^>/ {print($0)}; /^[^>]/ {print(toupper($0))}' developer/sh_refs_qiime_ver9_dynamic_29.11.2022_dev.fasta | tr -d ' ' > developer/sh_refs_qiime_ver9_dynamic_29.11.2022_dev_uppercase.fasta
#qiime tools import --type "FeatureData[Sequence]" --input-path ./sh_refs_qiime_ver9_dynamic_29.11.2022_dev_uppercase.fasta --output-path ./unite-ver9-seqs-fungi_29.11.2022.qza
#qiime tools import --type "FeatureData[Taxonomy]" --input-path ./sh_taxonomy_qiime_ver9_dynamic_29.11.2022_dev.txt --output-path ./unite-ver9-tax-fungi_29.11.2022.qza --input-format HeaderlessTSVTaxonomyFormat
#qiime feature-classifier fit-classifier-naive-bayes --i-reference-reads unite-ver9-seqs-fungi_29.11.2022.qza  --i-reference-taxonomy ./unite-ver9-tax-fungi_29.11.2022.qza --o-classifier unite-ver9-dynamic-classifier-fungi_29.11.2022.qza


#awk '/^>/ {print($0)}; /^[^>]/ {print(toupper($0))}' developer/sh_refs_qiime_ver9_dynamic_29.11.2022_dev.fasta | tr -d ' ' > developer/sh_refs_qiime_ver9_dynamic_29.11.2022_dev_uppercase.fasta
qiime tools import --type "FeatureData[Sequence]" --input-path ./fungi_unite_insdc.fasta --output-path ./fungi_unite_insdc_seq_latest.qza
qiime tools import --type "FeatureData[Taxonomy]" --input-path ./fungi_unite_insdc.txt --output-path ./fungi_unite_insdc_tax_latest.qza --input-format HeaderlessTSVTaxonomyFormat
qiime feature-classifier fit-classifier-naive-bayes --i-reference-reads ./fungi_unite_insdc_seq_latest.qza  --i-reference-taxonomy ./fungi_unite_insdc_tax_latest.qza --o-classifier fungi_unite_insdc_classifier_latest.qza
