options(java.parameters = "-Xmx64g", stringsAsFactors = F)
setwd("C:\Users\user\Desktop\Github\Sky0903\DWTPs\DWTPs")
set.seed(20230320)

library("phyloseq")
library("tidyverse")
library("ggplot2")
library("xlsx")
library("ggtree")
library("treeio")
library("ggstance")
library("vegan")
library("RColorBrewer")
library("dplyr")
library("ggforce")
library("reshape2")
library("knitr")
library("extrafont")
library("readxl")
library("patchwork")
library("MASS")

#Identity threshold : 97%
tap97 <- qiime2R::qza_to_phyloseq(features = "table_tap5.qza",
                                  tree = "rooted_tree_tap5.qza",
                                  metadata = "metadata.tsv",
                                  taxonomy = "tapwater_97_vsearch.qza")

#Checking the sequencing depth per sample 
tab <- otu_table(tap97)
class(tab) <- "matrix"
nice_colors = c("#999999", "#E69F00", "#56B4E9","#e98756","#c08160","#5970e6", "#CDDC49", "#C475D3","#E94B30", "#233F57", "#FEE659", "#A1CFDD", "#F4755E", "#D6F6F7","#EB6D58", "#6898BF", "#3300FF", "#0000CC", "#0099CC", "#66FF99")
tab <- t(tab)
out = rarecurve(tab, col=nice_colors, step=100 , lwd=2, main="Rarefaction Curve for all samples", labe=F)
names(out) = rownames(t(otu_table(tap97)))
protox <- mapply(FUN = function(x, y) {
  mydf <- as.data.frame(x)
  colnames(mydf) <- "value"
  mydf$SampleID <- y
  mydf$subsample <- attr(x, "Subsample")
  mydf
}, x = out, y = as.list(names(out)), SIMPLIFY = FALSE)
xy <- do.call(rbind, protox)
rownames(xy) <- NULL
xy = data.frame(xy, sample_data(tap97)[match(xy$SampleID, rownames(sample_data(tap97))), ])
ggplot(xy, aes(x = subsample, y = value, color = SampleID))+ theme_bw() 
+ theme(axis.title = element_text(face = "bold", family = "Tahoma", size = 18), axis.title.x = element_text(vjust = -0.25), axis.title.y = element_text(vjust = 2) ,axis.text = element_text(family = "Tahoma", size= 11.5) ,legend.title = element_text(family = "Tahoma", face = "bold", size = 18), legend.text = element_text(family = "Tahoma", size = 11.5)) +geom_line()
+ xlab("Sequenced Reads")
+ ylab('ASVs Detected')
+ guides(fill=guide_legend(title="Sample"))
+ geom_vline(xintercept =35549 , color = "red", size= 1, linetype = "dashed")

#After checking the relocation curve, remove samples that do not have sufficient sequencing depth and proceed with analysis
otu_mat<- read_excel("tap97.xlsx", 6)
tax_mat<- read_excel("tap97.xlsx", 8)
samples_df <- read_excel("tap97.xlsx", 7)
otu_mat <- otu_mat %>% tibble::column_to_rownames("ASV")
tax_mat <- tax_mat %>% tibble::column_to_rownames("ASV")
samples_df <- samples_df %>% tibble::column_to_rownames("Sample")
otu_mat <- as.matrix(otu_mat)
tax_mat <- as.matrix(tax_mat)
OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
TAX = tax_table(tax_mat)
samples = sample_data(samples_df)
tap97_new<- phyloseq(OTU, TAX, samples) 

tap97_edit_g2 <- tax_glom(tap97_new, taxrank = "Species") # agglomerate by species level
tap97_edit_g2.melt=psmelt(tap97_edit_g2)

otu_mat<- read_excel("tap97_edit.species.glom2.xlsx", 3)
tax_mat<- read_excel("tap97_edit.species.glom2.xlsx", 5)
samples_df <- read_excel("tap97_edit.species.glom2.xlsx", 4)
otu_mat <- otu_mat %>% tibble::column_to_rownames("ASV")
tax_mat <- tax_mat %>% tibble::column_to_rownames("ASV")
samples_df <- samples_df %>% tibble::column_to_rownames("Sample")
otu_mat <- as.matrix(otu_mat)
tax_mat <- as.matrix(tax_mat)
OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
TAX = tax_table(tax_mat)
samples = sample_data(samples_df)
tap97_edit.glom2<- phyloseq(OTU, TAX, samples)


#Decontamination by Djurhuus et al.
#Outlier removal, for PCR replicates that are very different from their sister replicates

allOTUs<-as.data.frame(otu_table(tap97_edit.glom2, taxa_are_rows=TRUE))
samples_all <- allOTUs
(replicates<-substr(colnames(samples_all),1,nchar(colnames(samples_all))-1))
distmat<-vegdist(t(samples_all))
distList=list(NA) ; distList.tri=list(NA); index=1
for (i in unique(replicates)){
  rowMatch<-which(replicates%in%i)
  distList[[index]]<-as.matrix(distmat)[rowMatch, rowMatch]
  distList.tri[[index]]<-	distList[[index]][upper.tri(	distList[[index]])]
  index=index+1
}

#visualize, if desired 
#hist(unlist(distList.tri))
#distList.tri 
normparams=fitdistr(unlist(distList.tri), "normal")$estimate #fit distribution to bray-curtis dissimilarities. lognormal, beta, normal, etc, may have less-good fits
probs=pnorm(unlist(distList.tri), normparams[1], normparams[2])
outliers =which(probs>0.95)
# minOutlier<-min(unlist(distList.tri)[outliers]) #minimum outlier value
minOutlier<-0.5  #here, if desired, instead of fitting a distribution you can a hard cutoff because of the distribution visualized
namesOutliers=list(NA)
for (i in 1:length(distList)){
  namesOutliers[[i]]<-intersect(
    names(which(colSums(distList[[i]]>=minOutlier)>0)),
    names(which.max(rowMeans(distList[[i]])))
  )
}
samples_all2<-samples_all[,-match(unlist(namesOutliers),colnames(samples_all))]
samples_all3 <- samples_all2[rowSums(samples_all2)>0,]


tap97_edit_p2<- prune_taxa(taxa_sums(tap97_edit.glom2) > 9, tap97_edit.glom2) # Remove ASVs represented low abundance

tap97.rare2 <- rarefy_even_depth(tap97_edit_p2, rngseed=20230320, sample.size=0.9*min(sample_sums(tap97_edit_p2)), replace=F)