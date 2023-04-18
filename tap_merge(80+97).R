#merge the annotations to maximize taxonomic resolution (identity threshold : 80% + 97%)
options(java.parameters = "-Xmx64g", stringsAsFactors = F)
setwd("C:/Users/user/Desktop/Github/Sky0903/DWTPs/DWTPs")
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



otu_mat<- read_excel("tap97.rare2.xlsx", 3)
tax_mat<- read_excel("tap97.rare2.xlsx", 6)
samples_df <- read_excel("tap97.rare2.xlsx", 4)
otu_mat <- otu_mat %>% tibble::column_to_rownames("OTU")
tax_mat <- tax_mat %>% tibble::column_to_rownames("OTU")
samples_df <- samples_df %>% tibble::column_to_rownames("Sample")
otu_mat <- as.matrix(otu_mat)
tax_mat <- as.matrix(tax_mat)
OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
TAX = tax_table(tax_mat)
samples = sample_data(samples_df)
tap_merge<- phyloseq(OTU, TAX, samples)

# Select only invertebrate 
tap_merge_inv <- subset_taxa(tap_merge, Phylum%in% c("Annelida", "Arthropoda", "Cnidaria", "Gastrotricha", "Mollusca", "Nematoda", "Platyhelminthes", "Rotifera"))

# Number of ASVs box plot
plot_alpha_invertebrate <- plot_richness(tap_merge_inv, x = "Stage", measures = c("Observed"))
+ geom_point( aes(fill = Stage), shape = 21, size = 3, color = "black")
+ theme_bw()
+ theme(strip.text.x = element_text(size = 18, face = "bold", family = "Tahoma") ,axis.title = element_text(face = "bold", family = "Tahoma", size = 18), axis.title.x = element_text(vjust = -0.25), axis.title.y = element_text(vjust = 2),axis.text.x = element_text(family = "Tahoma", size= 10.5) ,axis.text.y = element_text(family = "Tahoma", size= 11.5), legend.position = "none" , plot.margin = unit(c(1,6,1,6), "cm"))
+ geom_boxplot(aes(fill = Stage, group = Stage), alpha = 0.5, outlier.alpha = 0)
+ scale_fill_discrete(name="Stage")
+ ylab("Observed")

# Major class
Branchiopoda <- subset_taxa(tap_merge_inv, Class%in% c("Branchiopoda"))
Catenulida <- subset_taxa(tap_merge_inv, Class%in% c("Catenulida"))
Eurotatoria <- subset_taxa(tap_merge_inv, Class%in% c("Eurotatoria"))

class_alpha <- read.xlsx(file = "figure.xlsx", 1)
p_Ca <- ggplot(class_alpha, aes(x = Stage , y = Observed, fill = Stage))
+ theme_bw()
+ theme(text = element_text(family = "Tahoma"), axis.text = element_text(size= 13), axis.text.x = element_blank(), axis.title = element_text(size = 22, face = "bold"),axis.title.x = element_blank(), axis.title.y = element_text(vjust = 2),legend.text = element_text(size = 14), strip.text.x = element_text(size = 15, face = "bold", family = "Tahoma"), legend.title = element_text(family = "Tahoma", face = "bold" , size = 17) ,strip.background = element_rect(colour = "black", fill = "white"),  panel.border = element_rect(color = "black", fill = NA, linewidth= 0.5) ) + geom_boxplot(aes(group = Station), width = 7, alpha = 1, outlier.alpha = 0 )
+ facet_grid(~Phylum+Class, scales = "free", space = "free_y") 
+ labs(fill= "Stage")

class_Zs <- read.xlsx(file = "figure.xlsx", 2
                      class_Zs$Class <- factor(class_Zs$Class, levels = c("Branchiopoda", "Catenulida", "Eurotatoria"))
                     
p_Cz <- ggplot(class_Zs, aes(x=Stage, y=St.Zs, fill = Z_Score))
+ geom_tile(colour = "black",size = 0.58)+ theme(axis.text.x = element_text(family = "Tahoma", size = 13),axis.text.y = element_blank(),axis.title.x = element_text(family = "Tahoma", face = "bold", vjust=-0.25, size = 22), axis.title.y = element_text(face = "bold", family = "Tahoma", vjust = 6, size = 22 ), legend.title = element_text(family = "Tahoma", face = "bold", size = 18), legend.text = element_text(family = "Tahoma", size = 14 ) ,panel.background = element_blank(), panel.grid = element_blank(), strip.background = element_blank(), strip.text = element_blank(), axis.ticks.y = element_blank()) +scale_fill_gradientn(colours = c("snow2", "red"), breaks = c(0,1,2))
+ facet_grid(cols = vars(Class))
+ ylab("Abundance")
+ labs(fill = "Z-score")
+ coord_fixed(ratio = 4)

p_Ca/p_Cz

# Relative mean abundance of major taxa
fig_2b_sp <- read.xlsx(file ="figure.xlsx",4)
fig_2b_inv <- read.xlsx(file ="figure.xlsx",5)
fig_2b_phylum <- read.xlsx(file ="figure.xlsx",3)

fig_2b_sp$Species <-factor(fig_2b_sp$Species, levels = c("Daphnia galeata", "Stenostomum sp.", "Brachionus sp.", "Keratella cochlearis", "Synchaeta grandis"))
fig_2b_phylum$Phylum <- factor(fig_2b_phylum$Phylum, labels = c("Arthropoda","Platyhelminthes", "Rotifera"))

p_Sarea_abu  <- ggplot() 
+ geom_area(data = fig_2b_phylum, aes(x=Stage, y = Ratio, fill = Phylum, group = Phylum),alpha = 0.4,  position = position_dodge(width = 0.001))
+ scale_fill_discrete(labels = c("Arthropoda" ,"Platyhelminthes" , "Rotifera"))
+ geom_line(data = fig_2b_sp ,aes(x=Stage, y=Ratio, color =Species, group = Species))
+ geom_errorbar(data=fig_2b_sp, aes(x=Stage, ymin=Ratio-Ratio_STDEV, ymax=Ratio+Ratio_STDEV ,color =Species),width =0.15)
+ geom_point(data = fig_2b_sp, aes(x=Stage, y=Ratio, color = Species), size = 3,shape = 21, color = "black")
+ geom_point(data = fig_2b_sp, aes(x=Stage, y=Ratio, color = Species), size= 2.8,shape = 19)
+ scale_color_manual(values = c("Daphnia galeata" = "red", "Stenostomum sp." = "light green" , "Brachionus sp." = "deepskyblue" ,"Keratella cochlearis" = "dodgerblue3" ,"Synchaeta grandis" = "blue3"))
+ theme_bw()
+ theme(axis.title = element_text(face = "bold", family = "Tahoma", size = 18), axis.title.x = element_text(vjust = -0.25), axis.title.y = element_text(vjust = 2) ,axis.text = element_text(family = "Tahoma", size= 11.5) ,legend.title = element_text(family = "Tahoma", face = "bold", size = 18), legend.text = element_text(family = "Tahoma", size = 11.5))
+ scale_x_discrete(expand = c(0.02,0.02))
+ ylab("Read counts (%)")
+ guides( fill = guide_legend(order = 1), color = guide_legend(order = 0))

p_Sarea_abu 
+ geom_line(data = fig_2b_inv ,aes(x=Stage, y=Ratio, group = Invertebrate), size = 1.5)
+ geom_point(data = fig_2b_inv, aes(x=Stage, y=Ratio), size = 3, shape = 19)
+ geom_errorbar(data =fig_2b_inv, aes(x=Stage, ymin= Ratio-Ratio_STDEV, ymax = Ratio+Ratio_STDEV), width = 0.15)


#Comparison between eDNA index and Relative abundance over stages
fig_index <- read.xlsx(file = "figure.xlsx", 8)
fig_index$Species <- factor(fig_index$Species, c("Daphnia galeata", "Stenostomum sp.", "Brachionus sp.", "Keratella cochlearis", "Synchaeta grandis"))
ky_fig <- ggplot()
+ geom_line(data = fig_index, aes(x=Stage, y=Abundance, group = Type, color = Species, linetype= Type), size=1)
+ geom_point(data = fig_index, aes(x=Stage, y= Abundance, color = Species, shape = Type), size=2, show.legend = FALSE)
+ geom_line(data = fig_index, aes(x=Stage, y=Abundance, group = Type, color = Species, linetype= Type) ,size=1)
+ geom_point(data = fig_index, aes(x=Stage, y= Abundance, color = Species, shape = Type), size=2, show.legend = FALSE)
+ geom_errorbar(data=fig_index, aes(x=Stage, ymin=Abundance-Stdev, ymax=Abundance+Stdev ,color =Species),width =0.15)
+ geom_errorbar(data=fig_index, aes(x=Stage, ymin=Abundance-Stdev, ymax=Abundance+Stdev  ,color =Species),width =0.15)
+ facet_grid(Species~.)
+ scale_color_manual(values = c("Daphnia galeata" = "red", "Stenostomum sp." = "light green" , "Brachionus sp." = "deepskyblue" ,"Keratella cochlearis" = "dodgerblue3" ,"Synchaeta grandis" = "blue3"))
+ theme_bw()
+ theme(axis.title = element_text(face = "bold", family = "Tahoma", size = 18), axis.title.x = element_text(vjust = -0.25), axis.title.y = element_text(vjust = 2.47) ,axis.text = element_text(family = "Tahoma", size= 11.5) ,legend.title = element_text(family = "Tahoma", face = "bold", size = 18), legend.text = element_text(family = "Tahoma", size = 11.5),strip.background = element_rect(colour = "black", fill = "white"), strip.text.y =  element_blank(),  panel.border = element_rect(color = "black", fill = NA, linewidth= 0.5))
+ scale_x_discrete(expand = c(0.02,0.02))
+ ylab("Relative abundance")


ggsave("Index - Relative abundance.png", dpi=300, dev="png", height=6, width=10, units="in")


# Beta diversity plot
tap_merge_inv.ord <- ordinate(tap_merge_inv, "NMDS","bray")

plot_ord_invertebrate <- plot_ordination(tap_merge_inv, tap_merge_inv.ord,  type = "samples", color = "Stage")
+ theme_bw()
+ theme(axis.title = element_text(face = "bold", family = "Tahoma", size = 18), axis.title.x = element_text(vjust = -0.25), axis.title.y = element_text(vjust = 2) ,axis.text = element_text(family = "Tahoma", size= 11.5), legend.title = element_text(family = "Tahoma", face = "bold", size = 18), legend.text = element_text(family = "Tahoma", size = 13))
+ geom_point(size = 5, shape = 19)
+ geom_mark_ellipse(show.legend = FALSE)
+ scale_color_discrete(name = "Stage")

# Bray-curtis distance
otu_mat<- read_excel("tap_merge_invertebrate.xlsx", 4)
tax_mat<- read_excel("tap_merge_invertebrate.xlsx", 6)
samples_df <- read_excel("tap_merge_invertebrate.xlsx", 5)
otu_mat <- otu_mat %>% tibble::column_to_rownames("ASV")
tax_mat <- tax_mat %>% tibble::column_to_rownames("ASV")
samples_df <- samples_df %>% tibble::column_to_rownames("Sample")
otu_mat <- as.matrix(otu_mat)
tax_mat <- as.matrix(tax_mat)
OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
TAX = tax_table(tax_mat)
samples = sample_data(samples_df)
tap_merge_inv_fig3b<- phyloseq(OTU, TAX, samples)

hc2 <- phyloseq::distance(tap_merge_inv_fug3b, method = "bray")
hc.c2 <- hclust(hc2, method = "average")
par(cex=2)
plot(hc.c2, hang = -1, main = "Bray-Curtis distance", ylab = "", xlab= "", sub = "")

