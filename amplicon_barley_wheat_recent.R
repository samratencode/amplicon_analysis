##Load package
library(plyr)
library(dplyr)
library(conflicted)
library(phyloseq)
library(metacoder)
library(Biostrings)
library(aaronsaundeR)
library(MicEco)
library(microbiome)
library(MicrobiotaProcess)
library(microbiomeutilities)
library(ggplot2)
library(gplots)
library(ggpubr)
library(UpSetR)
library(vegan)
library(pairwiseAdonis)
library(micro4all)
library(DESeq2)
library(ComplexHeatmap)
library(tidyverse)
library(sna)
library(igraph)
library(centiserve)
library(microbiomeMarker)
library(ape)
library(ggClusterNet)
library(network)
library(phylosmith)
library(metagMisc)
library(qgraph)
library(WGCNA)
library(networkD3)
library(threejs)
library(gautils2)
library(openxlsx)
#library(xlsx)
library(RColorBrewer)
library(randomForest)
library(strucchange)
##Set-working-area
setwd("/Users/sash0009/Desktop/Poorva/ANALYSIS_barley")

##Profile-Phyloseq
feature <- read.table(file = "16S_feature-table.tsv", sep = "\t", header = T, row.names = 1, skip = 1, comment.char = "")
head(feature)
taxonomy <- read.table(file = "16S_taxonomy.tsv", sep = "\t", header = T ,row.names = 1)
head(taxonomy)
reference_seqs <- readDNAStringSet(file = "16S_dna-sequences.fasta",format = "fasta", nrec = -1L, skip = 0L, seek.first.rec = FALSE, use.names = TRUE)
head(reference_seqs)
##clean the taxonomy
#For-SILVA-UNITE-database
taxonomy$Taxon <- gsub("[a-z]__","", taxonomy$Taxon)
head(taxonomy)
tax <- taxonomy %>% separate("Taxon", c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep=";")
tax[is.na(tax)] <- ""
tax[tax=="__"] <- ""
head(tax)

#Metadata-loading
metadata <- read.delim2("metadata_bacteria.tsv")
head(metadata)
metadata1 <- read.delim2("metadata_fungi.tsv")
head(metadata1)

#OTU
otu = otu_table(as.matrix(feature), taxa_are_rows = TRUE)
rownames(otu) <- paste0("ASV", 1:nrow(otu))
head(otu)
#otu <- cbind(AsvId = rownames(otu), otu)
#rownames(otu) <- 1:nrow(otu)
#otu_df=as.data.frame(otu)
#write.xlsx(otu,'/Users/sash0009/Desktop/Poorva/ANALYSIS_barley/asv_fun_barley.xlsx',colNames = TRUE,rowNames=TRUE)
#write.xlsx(otu,'/Users/sash0009/Desktop/Poorva/ANALYSIS_barley/asv_fun_barley.xlsx',colNames = TRUE,rowNames=TRUE)

#Taxonomy
TAX = phyloseq::tax_table(as.matrix(tax[1:7]))
colnames(TAX) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
rownames(TAX) <- rownames(otu)
head(TAX)
#write.xlsx(TAX,'/Users/sash0009/Desktop/Poorva/ANALYSIS_barley/taxa_fun_barley.xlsx',colNames = TRUE,rowNames=TRUE)
#write.xlsx(TAX,'/Users/sash0009/Desktop/Poorva/ANALYSIS_barley/taxa_fun_barley.xlsx',colNames = TRUE,rowNames=TRUE)

#sequence
names(reference_seqs) <- rownames(TAX)
head(reference_seqs)
ref = as.data.frame(reference_seqs)
head(ref)
#write.xlsx(ref,'/Users/sash0009/Desktop/Poorva/ANALYSIS_barley/sequence_fun_barley.xlsx',colNames = TRUE,rowNames=TRUE)
#write.xlsx(ref,'/Users/sash0009/Desktop/Poorva/ANALYSIS_barley/sequence_fun_barley.xlsx',colNames = TRUE,rowNames=TRUE)

#Merging
ps0 <-phyloseq(otu,TAX,reference_seqs)
ps0

#Sampling
sam = sample_data(metadata)
samdf=data.frame(sample_data(metadata))
rownames(sam) <- sample_names(otu)
head(sam)

#Merging
ps1<-phyloseq(otu,TAX,reference_seqs,sam)
ps1
###Tree-generation
TREE = read.tree("16S_rooted_tree.nwk")
random_tree = rtree(ntaxa(ps1), rooted=TRUE, tip.label=taxa_names(ps1))
#Merging
ps2<-phyloseq(otu,TAX,sam,reference_seqs,random_tree)
ps2
nsamples(ps2)
sample_names(ps2)
sample_variables(ps2)

###Rare-fraction
nice_colors = c("#999999", "#E69F00", "#56B4E9","#e98756","#c08160","#5800e6", "#CDDC49", "#C475D3", "#E94B30", "#233F57", "#FEE659", "#A1CFDD", "#F4755E", "#D6F6F7","#EB6D58", "#6898BF")
rarecurve(t(otu_table(ps2)), ylab="Species (ASVs)",step=100,lwd=1.5,col = nice_colors,cex=0.5,label=F,main="Rarefaction Curve for all 32 Fungal samples")
ggsave(filename = "fun_rare_faction_barley.svg")

###Pruning
#phy_abund <- ps_prune(phy, min.samples = 5, min.reads = 10)
#tree_top_100 <- subset_taxa(ps2_relabund, Kingdom=="funteria")
#tree_top_100 <- prune_taxa(names(sort(taxa_sums(tree_top_100),TRUE)[1:100]), tree_top_100)

###Rarefied
set.seed(1) # keep result reproductive
pruned <- prune_samples(sample_sums(ps2)>0, ps2)
pruned
ps2.rarefied = rarefy_even_depth(pruned, rngseed=1, sample.size=0.9*min(sample_sums(pruned)), replace=F)
ps2.rarefied

##Venn
ps_venn(ps2.rarefied, group = "dsRNAFusarium",fraction = 0,relative = TRUE,weight = FALSE,plot = TRUE)
ggsave(filename = "bac_ven_dsrna_fusarium_barley.svg")

AsvId_common=common_taxa(ps2.rarefied, treatment = "Treatment", subset = NULL, n = "all")
AsvId_common_dsrna_fusarium=as.data.frame(AsvId_common)
AsvId_common_treatment=as.data.frame(AsvId_common)
AsvId_uniq=unique_taxa(ps2.rarefied, treatment = "Treatment", subset = NULL)

No_treatment_uniq_dsrna_fusarium=as.data.frame(AsvId_uniq$`1_No_treatment`)
Fusarium_graminearum_uniq_dsrna_fusarium=as.data.frame(AsvId_uniq$`2_Fusarium_graminearum`)
dsRNA_Fusarium_graminearum_uniq_dsrna_fusarium=as.data.frame(AsvId_uniq$`3_dsRNA_Fusarium_graminearum`)
dsRNA_uniq_dsrna_fusarium=as.data.frame(AsvId_uniq$`4_dsRNA`)


No_treatment_uniq_treatment=as.data.frame(AsvId_uniq$`1_No_treatment`)
No_treatment_Fg_uniq_treatment=as.data.frame(AsvId_uniq$`2_No_treatment_Fg`)
Non_specific_uniq_treatment=as.data.frame(AsvId_uniq$`3_Non_specific`)
Non_specific_Fg_uniq_treatment=as.data.frame(AsvId_uniq$`4_Non_specific_Fg`)
Cyp51_uniq_treatment=as.data.frame(AsvId_uniq$`5_Cyp51`)
Cyp51_Fg_uniq_treatment=as.data.frame(AsvId_uniq$`6_Cyp51_Fg`)
SdhB_uniq_treatment=as.data.frame(AsvId_uniq$`7_SdhB`)
SdhB_Fg_uniq_treatment=as.data.frame(AsvId_uniq$`8_SdhB_Fg`)

asvid_common_uniq <-rbind.fill(AsvId_common_treatment,No_treatment_uniq_treatment,No_treatment_Fg_uniq_treatment,Non_specific_uniq_treatment,Non_specific_Fg_uniq_treatment,Cyp51_uniq_treatment,Cyp51_Fg_uniq_treatment,SdhB_uniq_treatment,SdhB_Fg_uniq_treatment)
asvid_common_uniq1 <-rbind.fill(AsvId_common_dsrna_fusarium,No_treatment_uniq_dsrna_fusarium,Fusarium_graminearum_uniq_dsrna_fusarium,dsRNA_Fusarium_graminearum_uniq_dsrna_fusarium,dsRNA_uniq_dsrna_fusarium)

write.xlsx(asvid_common_uniq,file = "bac_AsvId_common_uniq_treatment_barley.xlsx", sep = "\t", quote = F, rowNames = F, colNames = T )

##upset
upsetdata <- get_upset(obj=ps2.rarefied,factorNames="dsRNAFusarium")
upset(upsetdata, empty.intersections = "on", order.by = "freq",sets.bar.color = "#56B4E9",matrix.color="blue")
upset(upsetdata, empty.intersections = "on", order.by = "freq",sets.bar.color = c("#56B4E9", "#CDDC49", "#E69F00", "#FEE659"),matrix.color="blue")
ggsave(filename = "fun_upset_TreatmentGroup_barley.svg")

co=taxa_core(ps2.rarefied, treatment = "Treatment", frequency = 0.2, abundance_threshold = 0.01)
otu_table(co)

tab_dom <- dominance(ps2.rarefied, index = "all")
tab_rar <- rarity(ps2.rarefied, index = "all")
tab_cov <- coverage(ps2.rarefied, threshold = 0.5)
tab_cor <- core_abundance(ps2.rarefied, detection = .1/100, prevalence = 50/100)
tab_inq <- inequality(ps2.rarefied)
tab_evn <- evenness(ps2.rarefied, "all")
head(tab_evn)
core.rel <- microbiome::transform(ps2.rarefied, "compositional")
core.rel.f <- format_to_besthit(core.rel)
core.rel.f <- microbiome::add_besthit(core.rel)
taxa_names(core.rel.f)[1:3]
core.taxa.standard <- core_members(core.rel.f, detection = 0.0001, prevalence = 60/100)
core.taxa.standard

#Set different detection levels and prevalence
det <- c(0, 0.1, 0.5, 2, 5, 20)/100
prevalences <- seq(.5, 1, .5) #0.5 = 95% prevalence
detections <- round(10^seq(log10(1e-3), log10(.2), length = 10),3)
#(1e-3) = 0.001% abundance; change "-3" to -2 to increase to 0.01%
p <- plot_core(core.rel.f, plot.type = "heatmap", colours = rev(brewer.pal(10, "Spectral")),min.prevalence = 0.5, prevalences = prevalences, detections = detections) + xlab("Detection Threshold (Relative Abundance (%))")
print(p)
ggsave(filename = "fun_core_prevalence_barley.svg")
p1 <-plot_core(core.rel.f, prevalences = prevalences, detections = det, plot.type = "lineplot") + xlab("Relative Abundance (%)") + theme_bw()
print(p1)

##Top-phyla

##Normalization
df_norm=ps2.rarefied %>% transform_sample_counts(function(x) {x/sum(x)}*100)  %>% psmelt() %>% arrange(OTU) %>% rename(AsvId = OTU) %>%  select(AsvId,Kingdom,Phylum, Class, Order, Family, Genus, Species,Sample,Treatment,dsRNAFusarium,Abundance)
df_norm_avg=df_norm %>% group_by (SamplingStage,Species) %>% summarise(Abundance = mean(Abundance)) %>% arrange(-Abundance) 
# sanity check: is total abundance of each sample 100%?
df_norm_avg %>% group_by(SamplingStage) %>% summarise(Abundance = sum(Abundance)) %>% pull(Abundance) %>% `==`(100) %>% all()
write.xlsx(df_norm_avg,file = "fun_relative_abundance_demo.xlsx", sep = "\t", quote = F, row.names = F, col.names = T )

write.xlsx(ps2.rarefied %>% transform_sample_counts(function(x) {x/sum(x)}*100)  %>% psmelt() %>% dplyr::arrange(OTU) %>% dplyr::rename(AsvId = OTU) %>%  dplyr::select(AsvId,Kingdom,Phylum, Class, Order, Family, Genus, Species,Sample,Treatment,dsRNAFusarium,Abundance), file = "fun_barley_relative_abundance.xlsx", sep = "\t", quote = F, rowNames = F, colNames = T)
write.table(ps2.rarefied %>% transform_sample_counts(function(x) {x/sum(x)}*100) %>% psmelt() %>% dplyr::arrange(OTU) %>% dplyr::rename(AsvId = OTU) %>%  dplyr::select(AsvId,Kingdom,Phylum, Class, Order, Family, Genus, Species,Sample,Treatment,dsRNAFusarium,Abundance), file = "fun_barley_relative_abundance.tsv", sep = "\t", quote = F, row.names = F, col.names = T)


##Specific-Normalization
ps2_glom <- tax_glom(ps2.rarefied, "Species")
ps2_relabund <- transform_sample_counts(ps2_glom, function(x) {x/sum(x)}*100)
ps2_relabund_merge <- merge_samples(ps2_relabund, "Treatment")
ps2_relabund_final <- transform_sample_counts(ps2_relabund_merge, function(x) {x/sum(x)}*100)
ps2_relabund_final
#sort(taxa_sums(ps2_relabund_final), TRUE)[1:50]/nsamples(ps2_relabund_final)
write.xlsx(ps2_relabund_final %>% psmelt() %>% dplyr::arrange(OTU) %>% dplyr::rename(AsvId = OTU) %>% select(AsvId,Phylum,Sample,Abundance), file ="fun_barley_relative_abundance_phylum_treatment.xlsx", sep = "\t", quote = F, rowNames = F, colNames = T)
write.xlsx(ps2_relabund_final %>% psmelt() %>% dplyr::arrange(OTU) %>% dplyr::rename(AsvId = OTU) %>% select(AsvId,Genus,Sample,Abundance), file ="fun_barley_relative_abundance_genus_treatment.xlsx", sep = "\t", quote = F, rowNames = F, colNames = T)
write.xlsx(ps2_relabund_final %>% psmelt() %>% dplyr::arrange(OTU) %>% dplyr::rename(AsvId = OTU) %>% select(AsvId,Species,Sample,Abundance), file ="fun_barley_relative_abundance_species_treatment.xlsx", sep = "\t", quote = F, rowNames = F, colNames = T)

comp_top_50 <- prune_taxa(names(sort(taxa_sums(ps2_relabund_final),TRUE)[1:50]), ps2_relabund_final)
comp_top_50


###Composition-plot
#Phylum
plot_bar(ps2_relabund_final,fill = "Phylum")+ labs(x = "", y="Relative Abundance (%)\n")
ggsave(filename = "fun_comp_phylum_treatment_barley.svg")
#Genus
plot_bar(ps2_relabund_final,fill = "Genus")+ labs(x = "", y="Relative Abundance (%)\n")
ggsave(filename = "fun_comp_genus_treatment_barley_top50.svg")
#Species
plot_bar(ps2_relabund_final,fill = "Species")+ labs(x = "", y="Relative Abundance (%)\n")
ggsave(filename = "fun_comp_species_treatment_barley_top50.svg")


#Diversity-plot(alpha)
richness <- estimate_richness(ps2.rarefied,measures = c ("Observed", "Shannon"))
richness_sig <- cbind(SampleID = rownames(richness), richness)
rownames(richness_sig) <- 1:nrow(richness_sig)
richness_group<-merge(richness_sig,samdf,by = 'SampleID') %>% as_tibble () %>% select (SampleID,Treatment,dsRNAFusarium,Observed,Shannon)
write.xlsx(richness_group,file = "fun_alpha_diversity_measure_barley.xlsx", sep = "\t", quote = F, rowNames = F, colNames = T )
#Normality-test
hist(richness_group$Observed, col='steelblue')
hist(richness_group$Shannon, col='steelblue')
shapiro.test(richness_group$Observed)
shapiro.test(richness_group$Shannon)
ggplot(richness_group, aes(x = Observed, fill = Treatment)) + geom_histogram(color='black', alpha=0.4, position='identity') + scale_fill_manual(values=c('red', 'blue', 'purple','green','yellow','magenta', 'grey','lightgreen'))
ggplot(richness_group, aes(x = Observed, fill = dsRNAFusarium)) + geom_histogram(color='black', alpha=0.4, position='identity') + scale_fill_manual(values=c('red', 'blue', 'purple','green'))
ggplot(richness_group, aes(x = Shannon, fill = Treatment)) + geom_histogram(color='black', alpha=0.4, position='identity') + scale_fill_manual(values=c('red', 'blue', 'purple','green','yellow','magenta', 'grey','lightgreen'))
ggplot(richness_group, aes(x = Shannon, fill = dsRNAFusarium)) + geom_histogram(color='black', alpha=0.4, position='identity') + scale_fill_manual(values=c('red', 'blue', 'purple','green'))
richness_group %>%group_by(Treatment) %>% dplyr::summarise(statistic = shapiro.test(richness_group$Observed)$statistic,p.value = shapiro.test(richness_group$Observed)$p.value)
richness_group %>%group_by(dsRNAFusarium) %>% dplyr::summarise(statistic = shapiro.test(richness_group$Observed)$statistic,p.value = shapiro.test(richness_group$Observed)$p.value)
richness_group %>%group_by(Treatment) %>% dplyr::summarise(statistic = shapiro.test(richness_group$Shannon)$statistic,p.value = shapiro.test(richness_group$Shannon)$p.value)
richness_group %>%group_by(dsRNAFusarium) %>% dplyr::summarise(statistic = shapiro.test(richness_group$Shannon)$statistic,p.value = shapiro.test(richness_group$Shannon)$p.value)
#Parametric
anova.observed = aov(richness_group$Observed ~Treatment, data=samdf)
anova.shannon = aov(richness_group$Shannon ~Treatment, data=samdf)
summary(anova.observed)
summary(anova.shannon)
#anova.sh2 = aov(richness$Shannon ~Treatment*dsRNAFusarium, data=samdf)
#summary(anova.sh2)
tukeyhsd.observed <- TukeyHSD(anova.observed)
tukeyhsd.shannnon <- TukeyHSD(anova.shannon)
tukeyhsd.observed_dff=as.data.frame(tukeyhsd.observed$Treatment)
tukeyhsd.shannon_dff=as.data.frame(tukeyhsd.shannnon$Treatment)
write.xlsx(tukeyhsd.observed_dff,file="fun_alpha_diversity_measure_observed_tukeyhsd_treatment_barley.xlsx",sheetName = "Sheet1",colNames = TRUE,rowNames = TRUE,append = FALSE,showNA = TRUE,password = NULL)
write.xlsx(tukeyhsd.shannon_dff,file="fun_alpha_diversity_measure_shannon_tukeyhsd_treatment_barley.xlsx",sheetName = "Sheet1",colNames = TRUE,rowNames = TRUE,append = FALSE,showNA = TRUE,password = NULL)
#Non-parametric
kruskal.test(richness_group$Observed ~ Treatment, data=samdf)
kruskal.test(richness_group$Shannon ~ Treatment, data=samdf)
wilcox.observed <- pairwise.wilcox.test(richness_group$Observed, samdf$Treatment,  p.adjust.method = "BH")
wilcox.shannnon <- pairwise.wilcox.test(richness_groupm$Shannon, samdf$Treatment,  p.adjust.method = "BH")
wilcox.observed_dff=as.data.frame(wilcox.observed$p.value)
wilcox.shannnon_dff=as.data.frame(wilcox.shannnon$p.value)
write.xlsx(wilcox.observed_dff,file="fun_alpha_diversity_measure_observed_wilcox_treatment_barley.xlsx", sheetName = "Sheet1",colNames = TRUE,rowNames = TRUE,append = FALSE,showNA = TRUE,password = NULL)
write.xlsx(wilcox.shannnon_dff,file="fun_alpha_diversity_measure_shannon_wilcox_treatment_barley.xlsx", sheetName = "Sheet1",colNames = TRUE,rowNames = TRUE,append = FALSE,showNA = TRUE,password = NULL)

#a_my_comparisons <- list( c("YY", "MM"), c("CA", "TT"), c("NN", "SS"))
symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))
#plot_richness(physeq1, x="Group", measures="Shannon", color = "Group") + geom_boxplot(alpha=0.6) + theme (legend.position="none", axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=12)) + stat_compare_means(method = "wilcox.test", comparisons = a_my_comparisons, label = "p.signif", symnum.args = symnum.args)
plot_richness(ps2.rarefied,x="Treatment",measures=c("Observed", "Shannon")) + geom_boxplot() + theme_classic() + theme(strip.background = element_blank(), axis.text.x.bottom = element_text(angle = 90),axis.text.x = element_text(size = 10)) + stat_compare_means(method = "anova")
plot_richness(ps2.rarefied,x="dsRNAFusarium",measures=c("Observed", "Shannon")) + geom_boxplot() + theme_classic() + theme(strip.background = element_blank(), axis.text.x.bottom = element_text(angle = 90),axis.text.x = element_text(size = 10)) + stat_compare_means(method = "anova",label = "p.signif",symnum.args = symnum.args) 
plot_richness(ps2.rarefied,x="Treatment",measures=c("Observed", "Shannon"),color = "HotSpring") + geom_boxplot() + theme_classic() + theme(strip.background = element_blank(),legend.position="none", axis.text.x.bottom = element_text(angle = 90),axis.text.x = element_text(size = 10)) + stat_compare_means(method = "wilcox.test",label = "p.signif", symnum.args = symnum.args,ref.group = ".all.")
ggsave(filename = "fun_alpha_diversity_treatment_barley.svg")

###Beta-diversity
##PCoA and PERMANOVA/ADONIS(Bray curtis)
dist = phyloseq::distance(ps2.rarefied, method="bray",binary = TRUE)
dist_adonis_al <-adonis2(dist ~dsRNAFusarium, data=samdf)
dist_adonis_al
summary(dist_adonis_al)
dist_adonis_pw <- pairwise.adonis(dist,samdf$dsRNAFusarium)
dist_adonis_pw
write.xlsx(dist_adonis_pw, file ="fun_beta_diversity_pcoa_adonis_pairwise_dsrna_fusarium_barley.xlsx", sep = "\t", quote = F, rowNames = F, colNames = T)
ordination = ordinate(ps2.rarefied, method="PCoA", distance=dist)
plot_ordination(ps2.rarefied, ordination, type = "samples",color = "Field", shape="FieldNature") + theme_classic() + theme(strip.background = element_blank()) + stat_ellipse (geom = "polygon",alpha=0.0 ) + geom_point(size=2.5)
plot_ordination(ps2.rarefied, ordination, type = "samples",color = "dsRNAFusarium") + theme_classic() + theme(strip.background = element_blank()) + stat_ellipse (geom = "polygon",alpha=0.0 ) + geom_point(size=2.5)

ggsave(filename = "fun_beta_diversity_pcoa_dsrna_fusarium_barley.svg")
#beta <- betadisper(dist, samdf$Treatment)
#permutest(beta)
#anova(beta)

##CAP-analysis
ordination_cap = ordinate(ps2.rarefied, method="CAP", distance=dist,formula = ~ dsRNAFusarium)
#anova(ordination_cap, by="terms", permu=999)
anova(ordination_cap)
RsquareAdj(ordination_cap)
plot_scree(ordination_cap)
plot_ordination(ps2.rarefied, ordination_cap, type = "samples",color = "Treatment", shape="dsRNAFusarium") + theme_classic() + theme(strip.funkground = element_blank()) + stat_ellipse (geom = "polygon",alpha=0.0 ) + geom_point(size=2.5)
plot_ordination(ps2.rarefied, ordination_cap, type = "samples",color = "dsRNAFusarium") + theme_classic() + theme(strip.funkground = element_blank()) + stat_ellipse (geom = "polygon",alpha=0.0 ) + geom_point(size=2.5)

#pcoa_phyloseq(ps2.rarefied,  'Sample', circle = TRUE, method = 'bray')
ggsave(filename = "fun_beta_diversity_cap_dsrna_fusarium_barley.svg")

##NMDS and ANOSIM(Binary Jaccard)
dist = phyloseq::distance(ps2.rarefied, method="bray", binary = TRUE)
dist_anosim_al <-anosim(dist,samdf$Treatment)
dist_anosim_al
summary(dist_anosim_al)
ordination = ordinate(ps2.rarefied, method="NMDS", distance=dist)
ordination$stress
plot_ordination(ps2.rarefied, ordination,type = "samples", color="Treatment",shape="dsRNAFusarium") + theme_classic() + theme(strip.funkground = element_blank()) + stat_ellipse (geom = "polygon",alpha=0.0 ) + geom_point(size=2.5)
plot_ordination(ps2.rarefied, ordination,type = "samples", color="dsRNAFusarium") + theme_classic() + theme(strip.funkground = element_blank()) + stat_ellipse (geom = "polygon",alpha=0.0 ) + geom_point(size=2.5)

#nmds_phyloseq(ps2.rarefied,  'Sample', circle = FALSE)
ggsave(filename = "fun_beta_diversity_nmds_treatment_barley.svg")

###Tree-plotting
phy_genus <- tax_glom(ps2_relabund, "Genus")
phy_genus
tree_top_100 <- prune_taxa(names(sort(taxa_sums(phy_genus),TRUE)[1:75]), phy_genus)

plot(TREE)
plot_tree(tree_top_100)
plot_tree(tree_top_100,nodelabf=nodeplotblank, size="abundance",color="Location", label.tips="Genus", ladderize="left", plot.margin=0.4,text.size=2.9,base.spacing=0.01)

###Heat-tree
obj_metacode <- parse_phyloseq(ps2.rarefied,class_regex = "(.*)", class_key = "taxon_name")

obj_metacode1<- parse_phyloseq(ps2.rarefied)
obj_metacode1$data$otu_table <- calc_obs_props(obj_metacode1,data = "otu_table", cols = obj_metacode1$data$sample_data$sample_id)
obj_metacode1$data$tax_table <- calc_taxon_abund(obj_metacode1,data = "otu_table", cols = obj_metacode1$data$sample_data$sample_id)
obj_metacode1$data$diff_table <- compare_groups(obj_metacode1,data = "tax_table",cols = obj_metacode1$data$sample_data$sample_id, groups = obj_metacode1$data$sample_data$HotSpring)
obj_metacode1_data_diff_table_df = as.data.frame(obj_metacode1$data$diff_table)
obj_metacode1_data_tax_data_df = as.data.frame(obj_metacode1$data$tax_data)
obj_metacode1_df <- merge(obj_metacode1_data_tax_data_df,obj_metacode1_data_diff_table_df, by = 'taxon_id')
write.xlsx(obj_metacode1_df,file = "fun_hot_spring_heattree.xlsx", sep = "\t", quote = F, rowNames = F, colNames = T )
#obj_metacode1$data$n_mean <- calc_group_mean(obj_metacode1,data = "tax_table",cols = obj_metacode1$data$sample_data$sample_id, groups = obj_metacode1$data$sample_data$FusariumTreatment)
#obj_metacode1$data$n_rsd <- calc_group_rsd(obj_metacode1,data = "tax_table",cols = obj_metacode1$data$sample_data$sample_id, groups = obj_metacode1$data$sample_data$FusariumTreatment)
obj_metacode1$data$n_samples <- calc_n_samples(obj_metacode1,data = "tax_table", cols = obj_metacode1$data$sample_data$sample_id)
#heat_tree(obj_metacode1, node_size = n_obs, node_color = log2_median_ratio, node_label = taxon_names, node_size_trans = "linear", node_size_axis_label = "OTU count", node_color_axis_label = "Number of reads", edge_color = n_samples, edge_color_axis_label = "Number of samples", layout = "davidson-harel",initial_layout = "reingold-tilford")
obj_metacode1 %>% heat_tree( node_size = n_obs, node_color = log2_median_ratio, node_label = taxon_names, node_label_size_range = c(0.006, 0.03), node_size_trans = "log10 area", node_size_axis_label = "ASV count", node_color_axis_label = "Log2 ratio of median proportions",layout = "davidson-harel",initial_layout = "reingold-tilford")
ggsave(filename = "fun_heattree.svg")


##For-network
ps2_glom <- tax_glom(ps2.rarefied, "Genus")
ps2_relabund_net <- transform_sample_counts(ps2.rarefied, function(x) {x/sum(x)})
#ps2_relabund_filter_net <- ps_prune(ps2_relabund_net, min.abundance = 0.001)
ps2_relabund_final_net <- subset_taxa(ps2_relabund_net, !is.na(Genus) & !Genus %in% c(""))
ps2_relabund_final_net


####Microbiome-Network-analysis-correlation-based-phylosmith

co_oc_table <- co_occurrence(ps2_relabund_final_net,treatment = "dsRNAFusarium",subset = '3_dsRNA_Fusarium_graminearum',rho = 0.8, p = 0.05, method = "spearman")
ig_net <- network_ps(ps2_relabund_final_net, treatment = "dsRNAFusarium", subset = '3_dsRNA_Fusarium_graminearum', co_occurrence_table = co_oc_table)
ig_net_layout <-network_layout_ps(ps2_relabund_final_net, treatment = "dsRNAFusarium", subset = '3_dsRNA_Fusarium_graminearum', co_occurrence_table = co_oc_table)
View(ig_net_layout)
write.xlsx(ig_net_layout,file="fun_genus_network_dsrna_fusarium_graminearum_dsrna_fusarium_barley.xlsx", sheetName = "Sheet1",colNames = TRUE,rowNames = TRUE,append = FALSE,showNA = TRUE,password = NULL)
co_occurrence_network(ps2_relabund_final_net,co_occurrence_table = co_oc_table, treatment = "dsRNAFusarium", subset = '3_dsRNA_Fusarium_graminearum', classification = 'Genus',cluster=FALSE)


ggsave(filename = "fun_genus_network_dsrna_fusarium_graminearum_dsrna_fusarium_barley.svg")

variable_correlation_network(ps2_relabund_final_net, treatment = "HotSpring", subset = c("Chitu_Hotspring", "Shalla_Hotspring"), classification = 'Genus')


##Topological-features
###each-node-toplogical-features
node_top<-function(ig){
  e.count<-ecount(ig_net)
  v.count<-vcount(ig_net)
  betweeness.centrality<-igraph::betweenness(ig,v=V(ig),directed = FALSE, weights = NA,nobigint = TRUE, normalized = FALSE)
  closeness.centrality<-igraph::closeness(ig,vids=V(ig),weights = NA, normalized = FALSE)
  node.transitivity<-igraph::transitivity(ig,type = c("local"), vids = NULL,weights = NULL)
  ec<-igraph::eigen_centrality(ig)
  node.degree<-igraph::degree(ig,v=V(ig),mode="all")
  hub.score <- igraph::hub_score(ig)$vector
  fgc <- cluster_fast_greedy(ig)
  wtc<-cluster_walktrap(ig)
  fgc.modularity<-modularity(fgc)
  wtc.modularity<-modularity(wtc)
  mean.distance<- mean_distance(ig)
  node.topology<-data.frame(e.count,v.count,betweeness.centrality,closeness.centrality,node.transitivity,node.degree,mean.distance,hub.score,fgc.modularity,wtc.modularity,ec)
  return(node.topology)
}
nod_top_feature <- node_top(ig_net)
write.xlsx(nod_top_feature,file="fun_genus_network_topology_dsrna_fusarium_graminearum_dsrna_fusarium_barley.xlsx", sheetName = "Sheet1",colNames = TRUE,rowNames = TRUE,append = FALSE,showNA = TRUE,password = NULL)

network_results = function(net) {
  E(net)$weight = NA
  #  modularity
  fc = cluster_fast_greedy(net,weights =NULL)
  modularity = modularity(net,membership(fc))
  # graph properties
  results = data.frame(clustering.coefficient=transitivity(net), # clustering.coefficient
                       modularity=modularity,# modularity
                       mean.degree= mean(igraph::degree(net)),# mean degree
                       size=length(E(net)),# size - number of edges
                       order=length(V(net)),# order - number of vertices
                       mean.distance=mean_distance(net),#average path length
                       norm.degree=mean(igraph::degree(net)/length(V(net))),#normalized degree
                       betweenness.centrality=centralization.betweenness(net)$centralization #betweenness centrality
  )
  return(results)
}

network_results(sc.se)

ZiPiPlot()
net_properties()
node_properties()
#Other-network
network.2(phy_genus)
#
qgraph(cor(otu),groups = sam)

#igraph-based-distance-based
ig <- make_network(ps2_relabund, dist.fun="jaccard", max.dist=0.8)
plot_network(ig, ps2_relabund, color="Treatment", line_weight=0.4, label=NULL,type='samples')
ig <- make_network(ps2_relabund, dist.fun="jaccard", max.dist=0.8,type='taxa')
plot_network(ig, ps2_relabund, color = "Phylum", line_weight=0.4, label=NULL,type='taxa')

# Use igraph to make the graph and find membership
karate <- make_graph(ig)
wc <- cluster_walktrap(ig)
members <- membership(wc)
# Convert to object suitable for networkD3
karate_d3 <- igraph_to_networkD3(ig, group = members)
# Create force directed network plot
forceNetwork(Links = karate_d3$links, Nodes = karate_d3$nodes, Source = 'source', Target = 'target', NodeID = 'name', Group = 'group')
net.js <- threejs::igraph2graphjs(ig)
graphjs(ig,layout=layout_with_fr(ig, dim=3))

#
page_rank()
##Layout
layout_randomly(ig_net)
layout_in_circle(ig_net)
layout_on_sphere(ig_net)
layout_with_fr(ig_net)
layout_with_kk(ig_net)
layout_with_lgl(ig_net)




###Lefse-to-test-for-differential-abundance-between-categories
lefse <- run_lefse (ps2.rarefied, wilcoxon_cutoff = 0.05,norm = "CPM",taxa_rank = "all", group = "Treatment", kw_cutoff = 0.05,multigrp_strat = TRUE,lda_cutoff = 4)
head(marker_table(lefse))

# bar plot
plot_ef_bar(lefse)
ggsave(filename = "fun_lda_bar_treatment_barley.svg")
# dot plot
plot_ef_dot(lefse)
ggsave(filename = "fun_lda_dot_treatment_barley.svg")

###Marker-like-lefse-edgeR
mm_edger <- run_edger(ps2.rarefied, group = "FieldNature", pvalue_cutoff = 0.05, p_adjust = "fdr")
plot_ef_bar(mm_edger)
plot_ef_dot(mm_edger)
#For-more-than-two-group
cid <- phyloseq::subset_samples(ps2.rarefied,Treatment %in% c("1_No_treatment","2_No_treatment_Fg","3_Non_specific","4_Non_specific_Fg","5_Cyp51","6_Cyp51_Fg","7_SdhB","8_SdhB_Fg"))
cid1 <- phyloseq::subset_samples(ps2.rarefied,dsRNAFusarium %in% c("1_No_treatment","2_No_treatment_Fg","3_dsRNA_Fusarium_graminearum","4_dsRNA"))
mm_edger_mg <- run_edger(cid,group = "Treatment", method  = "QLFT", pvalue_cutoff = 0.05,p_adjust = "fdr")
mm_edger_mg
write.xlsx(mm_edger_mg@marker_table,file="fun_edger_dsrna_fusarium_barley.xlsx", sheetName = "Sheet1",colNames = TRUE,rowNames = TRUE,append = FALSE,showNA = TRUE,password = NULL)
plot_ef_bar(mm_edger_mg)
ggsave(filename = "fun_edger_bar_treatment_barley.svg")
plot_ef_dot(mm_edger_mg)
ggsave(filename = "fun_edger_dot_treatment_barley.svg")

###Differential-abundance
ds = phyloseq_to_deseq2(ps2.rarefied, ~dsRNAFusarium)
ds = DESeq(ds)
resultsNames(ds)
colData(ds)
res = DESeq2::results(ds,contrast=c("dsRNAFusarium","3_dsRNA_Fusarium_graminearum","1_No_treatment"))
res = res[base::order(res$padj, na.last=NA), ]
alpha = 0.05
res_sig = res[(res$padj < alpha), ]
head(res_sig)
res_sig = cbind(as(res_sig, "data.frame"), as(tax_table(ps2.rarefied)[rownames(res_sig), ], "matrix"))
head(res_sig)
res_sig <- cbind(AsvId = rownames(res_sig), res_sig)
rownames(res_sig) <- 1:nrow(res_sig)
head(res_sig)
ggplot(res_sig, aes(x=Genus, y=-log10(pvalue), color=Phylum)) +geom_jitter(size=3, width = 0.2) + theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))
write.xlsx(res_sig,file="fun_diff_abundance_measure_deseq2_5_Cyp51_1_No_treatment_treatment_barley.xlsx", sheetName = "Sheet1",colNames = TRUE,rowNames = TRUE,append = FALSE,showNA = TRUE,password = NULL)
write.xlsx(res_sig,file="fun_diff_abundance_measure_deseq2_reference_diversified_field_nature.xlsx", sheetName = "Sheet1",colNames = TRUE,rowNames = TRUE,append = FALSE,showNA = TRUE,password = NULL)
#
res_sig_d= res_sig %>% select (AsvId,Kingdom,Phylum, Class, Order, Family, Genus, Species)
res_sig_df=as.data.frame(res_sig_d)
head(res_sig_df)

res1 = DESeq2::results(ds,tidy=TRUE,contrast=c("FieldNature","outside","reference")) %>% tbl_df() %>% dplyr::arrange(row) %>% dplyr::rename(AsvId = row)
head(res1)
res1_df=as.data.frame(res1)
head(res1_df)
goi <- res1$AsvId
stopifnot(all(goi %in% names(ds)))
head(goi)
tcounts <- t(log2((counts(ds[goi, ], normalized=TRUE, replaced=FALSE)+.5))) %>% merge(colData(ds), ., by="row.names") %>% gather(AsvId, Abundance, (ncol(.)-length(goi)+1):ncol(.)) %>% select(AsvId,FieldNature, Abundance)
head(tcounts)
tcounts_df=as.data.frame(tcounts)
head(tcounts_df)
res1_tcounts_df <- merge(res1_df,tcounts_df, by = 'AsvId')
head(res1_tcounts_df)
ggplot(res1_tcounts_df, aes(x=log2FoldChange , y=-log10(pvalue), col=FieldNature )) + geom_jitter(size=3, width = 0.2) + theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))
#
res1_sig_tcounts <- merge(res1_tcounts_df,res_sig_df,by = 'AsvId') %>% select (AsvId,Kingdom,Phylum,Class,Order,Family,Genus,Species,Abundance,baseMean,log2FoldChange,lfcSE,stat,pvalue,padj)
res1_sig_tcounts_df=as.data.frame(res1_sig_tcounts)
head(res1_sig_tcounts_df)






