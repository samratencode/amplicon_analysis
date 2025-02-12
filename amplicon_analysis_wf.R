##Load package
library(plyr)
library(dplyr)
library(tidyverse)
library(centiserve)
library(conflicted)
library(phyloseq)
library(metacoder)
library(Biostrings)
library(aaronsaundeR)
library(MicEco)
library(microeco)
library(microbiome)
library(MicrobiotaProcess)
library(microbiomeutilities)
library(micro4all)
library(microbiomeMarker)
library(ggpicrust2)
library(ggplot2)
library(gplots)
library(ggpubr)
library(ComplexHeatmap)
library(RColorBrewer)
library(UpSetR)
library(vegan)
library(pairwiseAdonis)
library(DESeq2)
library(ANCOMBC)
library(MiscMetabar)
library(ape)
library(mia)
library(sna)
library(igraph)
library(SpiecEasi)
library(SPRING)
library(NetCoMi)
library(ggClusterNet)
library(network)
library(phylosmith)
library(qgraph)
library(WGCNA)
library(networkD3)
library(threejs)
library(metagMisc)
library(gautils)
library(openxlsx)
library(randomForest)
library(strucchange)
library(upstartr)
library(ExploreMetabar)

##Set-working-area
setwd("/Users/xxxx/xxxx/xxxx")

##Overlapping

conflicts_prefer(phyloseq::tax_glom)
conflicts_prefer(phyloseq::transform_sample_counts)
conflicts_prefer(phyloseq::psmelt)
conflicts_prefer(phyloseq::plot_bar)

##Profile-Phyloseq
feature <- read.table(file = "16S_feature_table.tsv", sep = "\t", header = T, row.names = 1, skip = 1, comment.char = "")
head(feature)
taxonomy <- read.table(file = "16S_taxonomy.tsv", sep = "\t", header = T ,row.names = 1)
head(taxonomy)
reference_seqs <- readDNAStringSet(file = "16S_dna_sequences.fasta",format = "fasta", nrec = -1L, skip = 0L, seek.first.rec = FALSE, use.names = TRUE)
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
write.xlsx(otu,'Users/xxxx/xxxx/xxxx/asv_bac_barley.xlsx',colNames = TRUE,rowNames=TRUE)
#write.xlsx(otu,'/Users/xxxx/xxxx/xxxx/asv_fun_barley.xlsx',colNames = TRUE,rowNames=TRUE)

#Taxonomy
TAX = phyloseq::tax_table(as.matrix(tax[1:7]))
colnames(TAX) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
rownames(TAX) <- rownames(otu)
head(TAX)
write.xlsx(TAX,'/Users/xxxx/xxxx/xxxx/taxa_bac_barley.xlsx',colNames = TRUE,rowNames=TRUE)
#write.xlsx(TAX,'/Users/xxxx/xxxx/xxxx/taxa_fun_barley.xlsx',colNames = TRUE,rowNames=TRUE)

#sequence
names(reference_seqs) <- rownames(TAX)
head(reference_seqs)
ref = as.data.frame(reference_seqs)
head(ref)
write.xlsx(ref,'/Users/xxxx/xxxx/xxxx/sequence_bac_barley.xlsx',colNames = TRUE,rowNames=TRUE)
#write.xlsx(ref,'/Users/xxxx/xxxx/xxxx/sequence_fun_barley.xlsx',colNames = TRUE,rowNames=TRUE)

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
otu.rare = otu_table(ps2)
otu.rare_df = as.data.frame(t(otu.rare))
rarecurve(otu.rare_df, ylab="Species (ASVs)",step=1000,lwd=1.5,col = nice_colors,cex=0.5,label=F,main="Rarefaction Curve for all Bacterial samples")

ggsave(filename = "bac_rare_faction_pb.svg")

###Pruning
#phy_abund <- ps_prune(phy, min.samples = 5, min.reads = 10)
#tree_top_100 <- subset_taxa(ps2_relabund, Kingdom=="Bacteria")
#tree_top_100 <- prune_taxa(names(sort(taxa_sums(tree_top_100),TRUE)[1:100]), tree_top_100)

###Rarefied
set.seed(1) # keep result reproductive
pruned <- prune_samples(sample_sums(ps2)>0, ps2)
pruned
ps2.rarefied = rarefy_even_depth(pruned, rngseed=1, sample.size=0.9*min(sample_sums(pruned)), replace=F)
ps2.rarefied

##upset-plot
#c("#56B4E9", "#CDDC49", "#E69F00","#e98756","#c08160","#F4755E", "#D6F6F7","#EB6D58", "#6898BF")
upsetdata <- get_upset(obj=ps2.rarefied,factorNames="TreatmentGroup")
upsetdata <- cbind(AsvId = rownames(upsetdata), upsetdata)
upsetdata_df=as.data.frame(upsetdata)
write.xlsx(upsetdata_df,file = "bac_upset_treatmentgroup_pb.xlsx", sep = "\t", quote = F, rowNames = F, colNames = T )
upset(upsetdata, sets=unique(as.vector(sample_data(ps2.rarefied)$TreatmentGroup)),empty.intersections = "on", order.by = "freq",sets.bar.color = "#56B4E9",matrix.color="blue")
upset(upsetdata, sets=unique(as.vector(sample_data(ps2.rarefied)$TreatmentGroup)),empty.intersections = "on", order.by = "freq",sets.bar.color = c("#56B4E9", "#CDDC49", "#E69F00","#e98756"),matrix.color="blue")

ggsave(filename = "bac_upset_treatmentgroup_pb.svg")



##Core-microbiome-analysis
sample_variables(ps2.rarefied)
phy_1w <- subset_samples(ps2.rarefied, TreatmentGroup == "4_dsRNA")

core.rel <- microbiome::transform(phy_1w, "compositional")
core.rel.f <- format_to_besthit(core.rel)
taxa_names(core.rel.f)[1:5]

#Set different detection levels and prevalence
det <- c(0, 0.1, 0.5, 2, 5, 20)/100
prevalences <- seq(.5, 1, .5) #0.5 = 95% prevalence
detections <- round(10^seq(log10(1e-3), log10(.2), length = 10),3)
#(1e-3) = 0.001% abundance; change "-3" to -2 to increase to 0.01%
p <- plot_core(core.rel.f, plot.type = "heatmap", colours = rev(brewer.pal(10, "Spectral")),min.prevalence = 0.5, prevalences = prevalences, detections = detections) + xlab("Detection Threshold (Relative Abundance (%))")
print(p)
ggsave(filename = "bac_core_prevalence_treatmentgroup_dsrna_pb.svg")

p1 <-plot_core(core.rel.f, prevalences = prevalences, detections = det, plot.type = "lineplot") + xlab("Relative Abundance (%)") + theme_bw()
print(p1)



##Normalization
df_norm=ps2.rarefied %>% transform_sample_counts(function(x) {x/sum(x)}*100)  %>% psmelt() %>% arrange(OTU) %>% rename(AsvId = OTU) %>%  select(AsvId,Kingdom,Phylum, Class, Order, Family, Genus, Species,Sample,Genotype,Abundance)
df_norm_avg=df_norm %>% group_by (SamplingStage,Species) %>% summarise(Abundance = mean(Abundance)) %>% arrange(-Abundance) 
# sanity check: is total abundance of each sample 100%?
df_norm_avg %>% group_by(SamplingStage) %>% summarise(Abundance = sum(Abundance)) %>% pull(Abundance) %>% `==`(100) %>% all()
write.xlsx(df_norm_avg,file = "bac_relative_abundance_demo.xlsx", sep = "\t", quote = F, row.names = F, col.names = T )
write.xlsx(ps2.rarefied %>% transform_sample_counts(function(x) {x/sum(x)}*100)  %>% psmelt() %>% dplyr::arrange(OTU) %>% dplyr::rename(AsvId = OTU) %>%  dplyr::select(AsvId,Kingdom,Phylum, Class, Order, Family, Genus, Species,Sample,Treatment,TreatmentGroup,Abundance), file = "bac_pb_relative_abundancee.xlsx", sep = "\t", quote = F, rowNames = F, colNames = T)
write.table(ps2.rarefied %>% transform_sample_counts(function(x) {x/sum(x)}*100) %>% psmelt() %>% dplyr::arrange(OTU) %>% dplyr::rename(AsvId = OTU) %>%  dplyr::select(AsvId,Kingdom,Phylum, Class, Order, Family, Genus, Species,Sample,Treatment,TreatmentGroup,Abundance), file = "bac_pb_relative_abundancee.tsv", sep = "\t", quote = F, row.names = F, col.names = T)


##Specific-Normalization
ps2_glom <- tax_glom(ps2.rarefied, "Species")
ps2_relabund <- transform_sample_counts(ps2_glom, function(x) {x/sum(x)}*100)
ps2_relabund_merge <- merge_samples(ps2_relabund, "TreatmentGroup")
ps2_relabund_final <- transform_sample_counts(ps2_relabund_merge, function(x) {x/sum(x)}*100)
ps2_relabund_final
#sort(taxa_sums(ps2_relabund_final), TRUE)[1:50]/nsamples(ps2_relabund_final)
write.xlsx(ps2_relabund_final %>% psmelt() %>% dplyr::arrange(OTU) %>% dplyr::rename(AsvId = OTU) %>% select(AsvId,Phylum,Sample,Abundance), file ="bac_pb_relative_abundance_phylum_treatmentgroup.xlsx", sep = "\t", quote = F, rowNames = F, colNames = T)
write.xlsx(ps2_relabund_final %>% psmelt() %>% dplyr::arrange(OTU) %>% dplyr::rename(AsvId = OTU) %>% select(AsvId,Genus,Sample,Abundance), file ="bac_pb_relative_abundance_genus_treatmentgroup.xlsx", sep = "\t", quote = F, rowNames = F, colNames = T)
write.xlsx(ps2_relabund_final %>% psmelt() %>% dplyr::arrange(OTU) %>% dplyr::rename(AsvId = OTU) %>% select(AsvId,Species,Sample,Abundance), file ="bac_pb_relative_abundance_species_treatmentgroup.xlsx", sep = "\t", quote = F, rowNames = F, colNames = T)

comp_top_20 <- prune_taxa(names(sort(taxa_sums(ps2_relabund_final),TRUE)[1:20]), ps2_relabund_final)
comp_top_20


###Composition-plot
#Phylum
plot_bar(ps2_relabund_final,fill = "Phylum")+ labs(x = "", y="Relative Abundance (%)\n")
ggsave(filename = "bac_comp_phylum_treatmentgroup_pb.svg")
#Genus
plot_bar(comp_top_20,fill = "Genus")+ labs(x = "", y="Relative Abundance (%)\n")
ggsave(filename = "bac_comp_genus_treatmentgroup_pb_top20.svg")
#Species
plot_bar(comp_top_20,fill = "Species")+ labs(x = "", y="Relative Abundance (%)\n")
ggsave(filename = "bac_comp_species_treatmentgroup_pb_top20.svg")

###Alternative-composition-plot
#Phylum
plot_taxa_composition(ps2_relabund_final, taxonomic.level = "Phylum", plot.type = "barplot", transform = "compositional", palette = brewer.pal(12, "Paired")) + scale_y_percent() + labs(x = "", y="Relative Abundance (%)\n")
ggsave(filename = "fun_comp_phylum_timepoint_pp.svg")

#Genus
plot_taxa_composition(comp_top_15, taxonomic.level = "Genus", plot.type = "barplot", transform = "compositional", palette = brewer.pal(12, "Paired")) + scale_y_percent() + labs(x = "", y="Relative Abundance (%)\n")
ggsave(filename = "fun_comp_genus_timepoint_pp_top15.svg")

#Species
plot_taxa_composition(comp_top_15, taxonomic.level = "Species", plot.type = "barplot", transform = "compositional", palette = brewer.pal(12, "Paired")) + scale_y_percent() + labs(x = "", y="Relative Abundance (%)\n")
ggsave(filename = "fun_comp_species_timepoint_pp_top15.svg")


###alpha-diversity
richness <- estimate_richness(ps2.rarefied,measures = c ("Shannon"))
richness_sig <- cbind(SampleId = rownames(richness), richness)
rownames(richness_sig) <- 1:nrow(richness_sig)
richness_group<-merge(richness_sig,samdf,by = 'SampleId') %>% as_tibble () %>% select (SampleId,Treatment,TreatmentGroup,Shannon)
write.xlsx(richness_group,file = "bac_alpha_diversity_measure_pb.xlsx", sep = "\t", quote = F, rowNames = F, colNames = T )

#Normality-test
hist(richness_group$Shannon, col='steelblue')
shapiro.test(richness_group$Shannon)
ggplot(richness_group, aes(x = Shannon, fill = Treatment)) + geom_histogram(color='black', alpha=0.4, position='identity') + scale_fill_manual(values=c('red', 'blue', 'purple'))
richness_group %>%group_by(Treatment) %>% dplyr::summarise(statistic = shapiro.test(richness_group$Shannon)$statistic,p.value = shapiro.test(richness_group$Shannon)$p.value)


#Parametric (significant pvalue)
anova.shannon = aov(richness_group$Shannon ~TreatmentGroup, data=samdf)
summary(anova.shannon)
tukeyhsd.shannnon <- TukeyHSD(anova.shannon)
tukeyhsd.shannon_dff=as.data.frame(tukeyhsd.shannnon$TreatmentGroup)
write.xlsx(tukeyhsd.shannon_dff,file="bac_alpha_diversity_measure_shannon_tukeyhsd_treatmentgroup_pb.xlsx",sheetName = "Sheet1",colNames = TRUE,rowNames = TRUE,append = FALSE,showNA = TRUE,password = NULL)


#Non-parametric (non significant pvalue)
kruskal.test(richness_group$Shannon ~ TreatmentGroup, data=samdf)
wilcox.shannnon <- pairwise.wilcox.test(richness_group$Shannon, samdf$TreatmentGroup,  p.adjust.method = "BH")
wilcox.shannnon_dff=as.data.frame(wilcox.shannnon$p.value)

write.xlsx(wilcox.shannnon_dff,file="bac_alpha_diversity_measure_shannon_wilcox_treatmentgroup_pw.xlsx", sheetName = "Sheet1",colNames = TRUE,rowNames = TRUE,append = FALSE,showNA = TRUE,password = NULL)

#Diversity-plot(alpha)
plot_richness(ps2.rarefied,x="TreatmentGroup",measures=c("Shannon")) + geom_boxplot() + theme_classic() + theme(strip.background = element_blank(), axis.text.x.bottom = element_text(angle = 90),axis.text.x = element_text(size = 10)) + stat_compare_means(method = "anova",label.x = 1.5) 
plot_richness(ps2.rarefied,x="TreatmentGroup",measures=c("Shannon")) + geom_boxplot() + theme_classic() + theme(strip.background = element_blank(), axis.text.x.bottom = element_text(angle = 90),axis.text.x = element_text(size = 10)) + stat_compare_means(label.x = 1.5) 

ggsave(filename = "bac_alpha_diversity_treatmentgroup_pb.svg")


###Beta-diversity
##PCoA and PERMANOVA/ADONIS(Bray curtis)
dist = phyloseq::distance(ps2.rarefied, method="bray",binary = TRUE)
dist_adonis_al <-adonis2(dist ~TreatmentGroup, data=samdf)
dist_adonis_al
summary(dist_adonis_al)
dist_adonis_pw <- pairwise.adonis(dist,samdf$TreatmentGroup)
dist_adonis_pw
write.xlsx(dist_adonis_pw, file ="bac_beta_diversity_pcoa_adonis_pairwise_treatmentgroup_pb.xlsx", sep = "\t", quote = F, rowNames = F, colNames = T)
ordination = ordinate(ps2.rarefied, method="PCoA", distance=dist)
plot_ordination(ps2.rarefied, ordination, type = "samples",color = "Treatment", shape="TreatmentGroup") + theme_classic() + theme(strip.background = element_blank()) + stat_ellipse (geom = "polygon",alpha=0.0 ) + geom_point(size=2.5)
plot_ordination(ps2.rarefied, ordination, type = "samples",color = "TreatmentGroup") + theme_classic() + theme(strip.background = element_blank()) + stat_ellipse (geom = "polygon",alpha=0.0 ) + geom_point(size=2.5)

ggsave(filename = "bac_beta_diversity_pcoa_treatmentgroup_pb.svg")


##NMDS and ANOSIM(Binary Jaccard)
dist = phyloseq::distance(ps2.rarefied, method="bray", binary = TRUE)
dist_anosim_al <-anosim(dist,samdf$TreatmentGroup)
dist_anosim_al
summary(dist_anosim_al)
ordination = ordinate(ps2.rarefied, method="NMDS", distance=dist)
ordination$stress
dist_adonis_pw <- pairwise.adonis(dist,samdf$TreatmentGroup)
dist_adonis_pw
write.xlsx(dist_adonis_pw, file ="bac_beta_diversity_nmds_adonis_pairwise_treatmentgroup_pw.xlsx", sep = "\t", quote = F, rowNames = F, colNames = T)
plot_ordination(ps2.rarefied, ordination,type = "samples", color="Treatment",shape="TreatmentGroup") + theme_classic() + theme(strip.background = element_blank()) + stat_ellipse (geom = "polygon",alpha=0.0 ) + geom_point(size=2.5)
plot_ordination(ps2.rarefied, ordination,type = "samples", color="TreatmentGroup") + theme_classic() + theme(strip.background = element_blank()) + stat_ellipse (geom = "polygon",alpha=0.0 ) + geom_point(size=2.5)


ggsave(filename = "bac_beta_diversity_nmds_treatmentgroup_pw.svg")




##For-network
ps2_glom <- tax_glom(ps2.rarefied, "Phylum")
ps2_relabund_net <- transform_sample_counts(ps2.rarefied, function(x) {x/sum(x)})
ps2_relabund_filter_net <- ps_prune(ps2_relabund_net, min.abundance = 0.001)
ps2_relabund_final_net <- subset_taxa(ps2_relabund_filter_net, !is.na(Phylum) & !Phylum %in% c(""))
ps2_relabund_final_net


####Microbiome-Network-analysis-correlation-based-phylosmith

co_oc_table <- co_occurrence(ps2_relabund_final_net,treatment = "TreatmentGroup",subset = '4_dsRNA',rho = 0.8, p = 0.05, method = "spearman")
ig_net <- network_ps(ps2_relabund_final_net, treatment = "TreatmentGroup", subset = '4_dsRNA', co_occurrence_table = co_oc_table)
ig_net_layout <-network_layout_ps(ps2_relabund_final_net, treatment = "TreatmentGroup", subset = '4_dsRNA', co_occurrence_table = co_oc_table)
View(ig_net_layout)
write.xlsx(ig_net_layout,file="bac_phylum_network_treatmentgroup_dsrna.xlsx", sheetName = "Sheet1",colNames = TRUE,rowNames = TRUE,append = FALSE,showNA = TRUE,password = NULL)
co_occurrence_network(ps2_relabund_final_net,co_occurrence_table = co_oc_table, treatment = "TreatmentGroup", subset = '4_dsRNA', classification = 'Phylum',cluster=FALSE)


ggsave(filename = "bac_phylum_network_treatmentgroup_dsrna.svg")


##Topological-features
###each-node-toplogical-features
node_top<-function(ig){
  E(ig)$weight = NA
  e.count<-ecount(ig)
  size=length(E(ig)) #size- number of edges
  v.count<-vcount(ig)
  order=length(V(ig)) #order - number of vertices
  betweeness.centrality<-igraph::betweenness(ig,v=V(ig),directed = FALSE, weights = NA,nobigint = TRUE, normalized = FALSE)
  closeness.centrality<-igraph::closeness(ig,vids=V(ig),weights = NA, normalized = FALSE)
  node.transitivity<-igraph::transitivity(ig,type = c("local"), vids = NULL,weights = NULL) # clustering.coefficient
  ec<-igraph::eigen_centrality(ig)
  node.degree<-igraph::degree(ig,v=V(ig),mode="all")
  node.degree.mean <-mean(igraph::degree(ig,v=V(ig),mode="all")) # mean degree
  node.degree.norm=mean(igraph::degree(ig,v=V(ig),mode="all")/length(V(ig))) #normalized degree
  hub.score <- igraph::hub_score(ig)$vector
  fgc <- cluster_fast_greedy(ig,weights =NULL)
  wtc<-cluster_walktrap(ig,weights =NULL)
  fgc.modularity<-modularity(fgc,membership(fgc))
  wtc.modularity<-modularity(wtc,membership(wtc))
  mean.distance<- mean_distance(ig) #average path length
  node.topology<-data.frame(e.count,size,v.count,order,betweeness.centrality,closeness.centrality,node.transitivity,node.degree,node.degree.mean,node.degree.norm,mean.distance,hub.score,fgc.modularity,wtc.modularity,ec)
  return(node.topology)
}
nod_top_feature <- node_top(ig_net)
write.xlsx(nod_top_feature,file="bac_phylum_network_topology_treatmentgroup_dsrna_pb.xlsx", sheetName = "Sheet1",colNames = TRUE,rowNames = TRUE,append = FALSE,showNA = TRUE,password = NULL)


#ZiPi-plot
fgc <- cluster_fast_greedy(ig_net,weights =NULL)
#fgc.modularity<-modularity(fgc,membership(fgc))
ig_deg=igraph::degree(ig_net)
N_nodes=length(ig_deg)

Z=ig_deg
Z[]=0
P=Z

Membership=membership(fgc)
Seq=seq(1:N_nodes)
for(i in 1:N_nodes){
  L=Membership==Membership[i]         
  neighbs=neighbors(ig_net,i)               
  Kis=sum(L[neighbs])
  SUM=0
  SUMsq=0	
  SUMP=0
  Miv=Seq[L]
  for(j in 1:sum(L)){
    neighbsj=neighbors(ig_net,Miv[j])
    Kjs=sum(L[neighbsj])
    SUM=SUM+Kjs
    SUMsq=SUMsq+Kjs^2
  }
  Z[i]=(Kis-SUM/sum(L))/sqrt(SUMsq/sum(L)-(SUM/sum(L))^2)
  if(Kis-SUM/sum(L)==0){Z[i]=0}
  for(k in 1:max(Membership)){
    Lp=Membership==k
    Kisp=sum(Lp[neighbs])
    SUMP=SUMP+(Kisp/ig_deg[i])^2}
  P[i]=1-SUMP
}
attribute_node.group1=cbind(degree=ig_deg,module=Membership,Pi=P,Zi=Z)
attribute_node.group2=cbind(degree=ig_deg,module=Membership,Pi=P,Zi=Z)
attribute_node.group3=cbind(degree=ig_deg,module=Membership,Pi=P,Zi=Z)
attribute_node.group4=cbind(degree=ig_deg,module=Membership,Pi=P,Zi=Z)

attribute_node.group1_df = data.frame(attribute_node.group1)
attribute_node.group2_df = data.frame(attribute_node.group2)
attribute_node.group3_df = data.frame(attribute_node.group3)
attribute_node.group4_df = data.frame(attribute_node.group4)

write.xlsx(attribute_node.group1_df,file="bac_phylum_network_zipi_treatmentgroup_nt_pb.xlsx", sheetName = "Sheet1",colNames = TRUE,rowNames = TRUE,append = FALSE,showNA = TRUE,password = NULL)
write.xlsx(attribute_node.group2_df,file="bac_phylum_network_zipi_treatmentgroup_ntfg_pb.xlsx", sheetName = "Sheet1",colNames = TRUE,rowNames = TRUE,append = FALSE,showNA = TRUE,password = NULL)
write.xlsx(attribute_node.group3_df,file="bac_phylum_network_zipi_treatmentgroup_dsrnafg_pb.xlsx", sheetName = "Sheet1",colNames = TRUE,rowNames = TRUE,append = FALSE,showNA = TRUE,password = NULL)
write.xlsx(attribute_node.group4_df,file="bac_phylum_network_zipi_treatmentgroup_dsrna_pb.xlsx", sheetName = "Sheet1",colNames = TRUE,rowNames = TRUE,append = FALSE,showNA = TRUE,password = NULL)


#zipi-graph-on-one-plot
par(mfrow=c(1,1),mar=c(4,4,2,8))
plot(attribute_node.group1[,3],attribute_node.group1[,4],xlim=c(0,1),ylim=c(-4,4),xlab="Among-module connectivity (Pi)",ylab=("Within-module connectivity (Zi)"),col=2,pch=1,cex=0.8)
abline(v=0.62,h=2.5,col=8)
points(attribute_node.group2[,3],attribute_node.group2[,4],col=3,pch=6,cex=0.8)
points(attribute_node.group3[,3],attribute_node.group3[,4],col=4,pch=6,cex=0.8)
points(attribute_node.group4[,3],attribute_node.group4[,4],col="grey27",pch=6,cex=0.8)
text(0.15,4,"Module hubs")
text(0.8,4,"Network hubs")
text(0.15,-4,"Peripherals")
text(0.8,-4,"Connectors")
legend(1.05,4,legend=c("1_NT","2_NT_Fg","3_dsRNA_Fg","4_dsRNA"),pch=c(1,6),col=c(2,3,4,"grey27"),xpd=T,bty="n",pt.lwd = 2)
legend(1.05,4,legend=c("3_dsRNA_Fg","4_dsRNA"),pch=c(1,6),col=c(2,"grey27"),xpd=T,bty="n",pt.lwd = 2)

ggsave(filename = "bac_phylum_network_zipi_plot.svg")



###Lefse-to-test-for-differential-abundance-between-categories
lefse <- run_lefse (ps2.rarefied, wilcoxon_cutoff = 0.05,norm = "CPM",taxa_rank = "all", group = "Treatment", kw_cutoff = 0.05,multigrp_strat = TRUE,lda_cutoff = 3)
head(marker_table(lefse))
write.xlsx(lefse@marker_table,file="bac_lefse_treatmentgroup_pb.xlsx", sheetName = "Sheet1",colNames = TRUE,rowNames = TRUE,append = FALSE,showNA = TRUE,password = NULL)

# bar plot
plot_ef_bar(lefse)
ggsave(filename = "bac_lda_bar_treatmentgroup_pb.svg")
# dot plot
plot_ef_dot(lefse)
ggsave(filename = "bac_lda_dot_treatmentgroup_pb.svg")

###ANCOMBC-Differential-abundance-analysis

ps_ancom = subset_samples(ps2.rarefied, TreatmentGroup %in% c("0_ND","4_dsRNA"))

ancom_out = ancombc(ps_ancom, formula = "TreatmentGroup", p_adj_method = "fdr",group = "TreatmentGroup")

results_ancom_bc = data.frame(ASVId = ancom_out$res$lfc$taxon,lfc = ancom_out$res$lfc$TreatmentGroup4_dsRNA, se = ancom_out$res$se$TreatmentGroup4_dsRNA, W = ancom_out$res$W$TreatmentGroup4_dsRNA, p_val = ancom_out$res$p_val$TreatmentGroup4_dsRNA, q_value = ancom_out$res$q_val$TreatmentGroup4_dsRNA, Diff_ab = ancom_out$res$diff_abn$TreatmentGroup4_dsRNA)
results_ancom_bc$lfc = results_ancom_bc$lfc * -1
results_ancom_bc
results_ancom_bc$group = ifelse(results_ancom_bc$q_value < 0.05 & results_ancom_bc$lfc > 0, "4_dsRNA", "Not singificant")
results_ancom_bc$group = ifelse(results_ancom_bc$q_value < 0.05 & results_ancom_bc$lfc < 0, "0_ND",results_ancom_bc$group)
results_ancom_bc
##Volcano-plot
ggplot(results_ancom_bc, aes(x = as.numeric(lfc), y = -log10(as.numeric(q_value)), color = group)) + geom_point() + labs(x = "Log2 Fold Change", y = "-log10(p-value)", title = "Volcano Plot") + geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +theme_bw()

TAXtab = as.data.frame(TAX)
TAXtab
name_tmp = sprintf("ASV%s",seq(1:nrow(TAXtab)))
TAXtab$name <- NA
TAXtab$name = ifelse(is.na(TAXtab$name), TAXtab$Genus, TAXtab$name)
TAXtab$name = ifelse(is.na(TAXtab$name), TAXtab$Family, TAXtab$name)
TAXtab$name = ifelse(is.na(TAXtab$name), TAXtab$Order, TAXtab$name)
TAXtab$name = ifelse(is.na(TAXtab$name), TAXtab$Class, TAXtab$name)
TAXtab$name = ifelse(is.na(TAXtab$name), TAXtab$Phylum, TAXtab$name)
unique(is.na(TAXtab$name)) 
TAXtab$name = paste(name_tmp, TAXtab$name, sep = "-")
TAXtab$name[1:10]

results_ancom_bc$name = TAXtab$name[match(results_ancom_bc$ASVId,rownames(TAXtab))]
results_ancom_bc
results_ancom_bc_df = data.frame(results_ancom_bc)

write.xlsx(results_ancom_bc_df,file="bac_ancombc_fc_0_ND_4_dsRNA_treatmentgroup_pp.xlsx", sheetName = "Sheet1",colNames = TRUE,rowNames = TRUE,append = TRUE,showNA = TRUE,password = NULL)

results_sig = results_ancom_bc[results_ancom_bc$q_value < 0.05,]
results_sig
table(results_sig$group)

Sum_abundance= colSums(decostand(ps2.rarefied@otu_table, method = "total"))
Sum_abundance
results_sig$sum_abundance = Sum_abundance[match(results_sig$ASVId, names(Sum_abundance))]
results_sig = results_sig[base::order(results_sig$sum_abundance, decreasing = TRUE),]
results_sig = results_sig[1:9,]
results_sig
results_sig$name = factor(results_sig$name, levels = rev(results_sig$name)) 

##ancombc-plot
a = ggplot(results_sig, aes(lfc, name, color = group ))+ geom_point() +theme_bw() + geom_segment( aes(x=0, xend=lfc, y=name, yend=name, color = group)) + geom_vline(xintercept = 0, size=0.3) +xlab("log2FoldChange") + ylab(NULL)
print(a)

ggsave(filename = "bac_ancombc_fc_0_ND_4_dsRNA_treatmentgroup_pp.svg")





