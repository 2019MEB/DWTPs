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

# Limit of detection (LOD) and limit of quantification (LOQ) for the universal primers 
lm_Cq <- lm(Cq~SQ, data = Cq_Standard)

Cq <- ggplot(lm_Cq, aes(x= log10(SQ), y = Cq))
+ geom_point(fill = "grey64", shape = 21, size = 5, color = "black")
+ theme_bw() + theme(text = element_text(family = "Tahoma"), axis.text = element_text(size= 13), axis.title = element_text(size = 18, face = "bold", family = "Tahoma"), axis.title.y = element_text(vjust = 2), axis.title.x = element_text(vjust = -0.25))
+ xlab(expression(bold(paste('Standard concentrations (log'['10']~ 'copies/reaction)'))))
+ ylab("Cq-value")
+ stat_smooth(method = "lm", color = "black", size = 1.25)



Cq_final <- Cq
+ geom_vline(xintercept = log10(6.89), color = "red", linetype = "dashed" ,size= 1)
+ geom_vline(xintercept = log10(14), color = "black", size= 1, linetype = "dashed")
+ geom_text(x = log10(8.3), y=38.75, label = "LOD", size = 5, angle = 90, color = "red", family = "Tahoma")
+ geom_text(x = log10(16.75), y=38, label = "LOQ", size = 5, angle = 90, color = "black", family = "Tahoma",)
+ geom_text(x= log10(100), y= 38.5, label = "R-Squared : 0.9808", size =5, family = "Tahoma")
+ geom_text(x= log10(100), y= 37.75, label = "y = -3.253 x + 40.66", size =5, family = "Tahoma")

ggsave("Cq_standard.png", dpi=300, dev="png", height=6, width=9, units="in")



# Read counts of C. riparius COI sequence by Illumina
m <- nls(Abundance~a*exp(b*Sample), NGS_standard ,start = list(a = 1, b=1) )

summary(m) # Check the equation

CR <- ggplot(NGS_standard, aes(x= Sample, y = Abundance))
+ geom_point(fill = "grey64", size = 5, shape = 21, color = "black")
+ geom_text(x=1, y=23450, label = expression(paste('y = 3.456 e'^'1.574x')), size = 5, family = "Tahoma")
+ geom_text(x=1, y=25000, label="R-Squared : 0.9867 ", size = 5, family = "Tahoma")
+ stat_smooth(size = 1.25, method = "nls",  formula = (y ~ a*exp(b*x)) ,se=F ,color = "black", method.args = list(start   = list(a = 1, b = 1), control = list(minFactor = 1/ 1024, maxiter = 100)))
+ theme_bw()
+ theme(text = element_text(family = "Tahoma"), axis.text = element_text(size= 13), axis.title = element_text(size = 19, face = "bold", family = "Tahoma"), axis.title.y = element_text(vjust = 2), axis.title.x = element_text(vjust = -0.25))
+ scale_x_continuous(labels = c("0","1","2","3","4","5"), breaks = c(0,1,2,3,4,5))
+ xlab(expression(bold(paste('Standard concentrations (log'['10']~ 'copies/reaction)'))))
+ ylab("Read counts")

ggsave("NGS_Read counts", dpi=300, dev="png", height=6, width=10, units="in")


# A standard curve was generated based on the copy numbers of synthetic dsDNA fragments of C. riparius COI.
EX <-read.xlsx(file = "standard.xlsx", 1)

lm.mc <- lm(Ct~Lgcopy, data = EX)

summary(lm.mc) 

mc.st <- ggplot(lm.mc, aes(x= Lgcopy, y = Ct))
+ geom_point(shape = 21, size = 5 ,color = "black", fill = "grey64")
+ theme_bw()
+ theme(text = element_text(family = "Tahoma"), axis.text = element_text(size= 13), axis.title = element_text(size = 18, face = "bold", family = "Tahoma"), axis.title.y = element_text(vjust = 2), axis.title.x = element_text(vjust = -0.25))
+ xlab(expression(bold(paste('Standard concentrations (log'['10']~ 'copies/reaction)'))))
+ ylab("Cq-value")
+ stat_smooth(method = "lm", color = "black", size = 1.25)
+ scale_y_continuous(limits = c(20,40), breaks = seq(20,40,5))

mc_final <- mc.st
+ geom_text(x = 3.25, y = 36.75, label = "R-Squared : 0.9976", size =6, family ="Tahoma")
+ geom_text(x = 3.25, y = 35.75, label = "y = -3.19133 x + 41.47285", size = 6, family ="Tahoma")

# Microcosms
RM <- read.xlsx(file = "riparius.xlsx", 1)
RM_new <- RM
RM_new$Time <- as.character(RM$Time)
my_y_title <- expression(paste(italic("C. riparius"), " DNA Copies / L water")) 

RMP <- ggplot(RM_new, aes(x= Time, y = Average , group = Sample))
+ geom_point(aes(shape = Sample, color = Sample), size = 4)
+ scale_color_manual(values = c("black", "purple", "darkorange3"), label = c("9 individuals","6 individuals", "3 individuals"))
+ scale_shape_manual(values=c(18, 15, 4), label = c("9 individuals","6 individuals", "3 individuals"))
+ theme_classic()
+ theme(text = element_text(family = "Tahoma"), axis.text = element_text(size= 13) , axis.title = element_text(size = 18, face = "bold"),axis.title.x = element_text(vjust = -0.25), axis.title.y = element_text(vjust = 2) ,legend.text = element_text(size = 13) ,legend.title = element_blank())
+ geom_errorbar(data = RM_new, aes(x = Time, ymin = Average-SEM, ymax=Average+SEM,width = 0.15, color = Sample), size =1)
+ scale_x_discrete(limits= c("24", "48","96"))
+ ylab(expression(bold(paste(bolditalic("C. riparius"), " DNA Copies / L water")))) + xlab("Time (hours)")

RMP2 <- RMP
+ scale_y_log10(limits=c(10^0, 10^6), breaks=10^(0:7), label=trans_format("log10",math_format(10^.x)))
+ guides(colour = guide_legend(reverse = TRUE), shape = guide_legend(reverse = TRUE))
+ annotation_logticks(sides = "l")

RMP3 <- RMP2
+ annotate("segment",0.5,xend = 3.5,6.89,yend = 6.89, linetype = "dashed", size = 1, color = "red")
+ annotate("segment",0.5,xend = 3.5,14,yend = 14, linetype = "dashed", size = 1)
+ annotate("text", x = 2.67, y = 5, label = "LOD for PCR assay (6.89 copies/reaction)", size= 4, color = "red", family = "Tahoma")
+ annotate("text", x = 2.71, y = 22, label = "LOQ for PCR assay (14 copies/reaction)", size= 4, family = "Tahoma")




