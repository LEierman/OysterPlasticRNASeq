## Open edgeR package and vegan package

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("edgeR")

library(edgeR)
library(vegan)

#### Two files - gonad samples in one and gill samples in the other
#### with all genes present.
### First is the gonad samples

x <- read.delim("all_genes_gonad_samples.txt", row.names="RefName")

head(x)

## Create the DGEList data class for analysis in edgeR

y <- DGEList(counts=x)
y$samples


# reading in data to identify groups for the sample data
targets <- read.table("sample_info_gonad.txt", header=TRUE)

# creating a model matrix for the GLM
design <- model.matrix(~Sex*Substrate, data=targets)
colnames(design)

design

##Filter low read count genes out of analysis
keep <- filterByExpr(y, design)
table(keep)

## Reset DGElist object to only have genes with sufficient read counts

y <- y[keep, , keep.lib.sizes=FALSE]
y$samples


AveLogCPM <- aveLogCPM(y)
hist(AveLogCPM)

## Normalize read counts

y <- calcNormFactors(y)
y$samples
y.counts <- y$counts
head(y.counts)
write.table(y.counts, "TMM_Gonad_Samples_Counts.txt",quote=FALSE)


#estimating dispersions for the GLM model
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTrendedDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)

# GLM fit to the fit the negative binomial GLM for each tag

fit <- glmFit(y,design)
colnames(fit)

# glmLRT to carry out likelihood ratio test; 
# select one or more coefficients to drop from full model to serve as a null model
# against which the full model is compared in a likelihood ratio test


lrtSE <- glmLRT(fit, coef=2)
lrtSE.all <- topTags(lrtSE, n=66625)
head(lrtSE.all)
write.table(lrtSE.all, "lrtSE_Se*SuModel_gonad_samples_only_TMM.txt", sep="\t")
is.de.lrtSE <- decideTestsDGE(lrtSE)
summary(is.de.lrtSE)
plotMD(lrtSE, status=is.de.lrtSE)

lrtSU <- glmLRT(fit, coef=3)
lrtSU.all <- topTags(lrtSU, n=66625)
head(lrtSU.all)
write.table(lrtSU.all, "lrtSU_Se*SuModel_gonad_samples_only_TMM.txt", sep="\t")
is.de.lrtSU <- decideTestsDGE(lrtSU)
summary(is.de.lrtSU)
plotMD(lrtSU, status=is.de.lrtSU)

lrtSExSU <- glmLRT(fit, coef=4)
lrtSExSU.all <- topTags(lrtSExSU, n=66625)
head(lrtSExSU.all)
write.table(lrtSExSU.all, "lrtSExSU_Se*SuModel_gonad_samples_only_TMM.txt", sep="\t")
is.de.lrtSExSU <- decideTestsDGE(lrtSExSU)
summary(is.de.lrtSExSU)
plotMD(lrtSExSU, status=is.de.lrtSExSU)

## interaction term analysis using contrasts
# first creating design object and groups
Group <- factor(paste(targets$Sex, targets$Substrate, sep="."))
cbind(targets, Group=Group)
design <- model.matrix(~0+Group)
colnames(design) <- levels(Group)
table(Group)

design

##Filter low read count genes out of analysis
keep <- filterByExpr(y, design)
table(keep)

## Reset DGElist object to only have genes with sufficient read counts

y <- y[keep, , keep.lib.sizes=FALSE]
y$samples


AveLogCPM <- aveLogCPM(y)
hist(AveLogCPM)

## Normalize read counts

y <- calcNormFactors(y)
y$samples

# next generation dispersions
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTrendedDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)

plotMDS(y)
pch <- c(0,15,1,16)
colors <- rep(c("dodgerblue2","dodgerblue2","red","red"))
parameter<-plotMDS(y, col=colors[Group], pch=pch[Group],cex=2.1,cex.axis=1.4,cex.lab=1.4)

legend("topright", legend=levels(Group), pch=pch, col=colors, cex=1.2, ncol=1)


#PERMANOVA using adonis function in vegan package
y_2<-t(y.counts)
adonis2(y_2~Sex*Substrate, data=targets, permutations=999)

# GLM fit to the fit the negative binomial GLM for each tag
fit <- glmFit(y,design)
colnames(fit)

# creating contrasts for comparisons within the interaction term

my.contrasts <-makeContrasts(
  SeM.SuSvsP=Male.Shell-Male.Plastic,
  SeF.SuSvsP=Female.Shell-Female.Plastic,
  SuS.SeMvsF=Male.Shell-Female.Shell,
  SuP.SeMvsF=Male.Plastic-Female.Plastic,
  Full=(Male.Shell-Female.Shell)-(Male.Plastic-Female.Plastic),
  MSvsFP = Male.Shell-Female.Plastic,
  MPvsFS = Male.Plastic-Female.Shell,
  Sex=(Male.Shell-Female.Shell)-(Female.Plastic-Male.Plastic),
  Substrate=(Male.Shell-Male.Plastic)-(Female.Plastic-Female.Shell),
  levels=design)

lrt1 <- glmLRT(fit, contrast=my.contrasts[,"SeM.SuSvsP"])
lrt1.all <- topTags(lrt1, n=73220)
head(lrt1.all)
write.table(lrt1.all, "lrtSeM.SuSvsP_gonad_samples_TMM.txt", sep="\t")
is.de.lrt1 <- decideTestsDGE(lrt1)
summary(is.de.lrt1)
plotMD(lrt1, status=is.de.lrt1)

lrt2 <- glmLRT(fit, contrast=my.contrasts[,"SeF.SuSvsP"])
lrt2.all <- topTags(lrt2, n=73220)
head(lrt2.all)
write.table(lrt2.all, "lrtSeF.SuSvsP_gonad_samples_TMM.txt", sep="\t")
is.de.lrt2 <- decideTestsDGE(lrt2)
summary(is.de.lrt2)
plotMD(lrt2, status=is.de.lrt2)

lrt3 <- glmLRT(fit, contrast=my.contrasts[,"SuS.SeMvsF"])
lrt3.all <- topTags(lrt3, n=73220)
head(lrt3.all)
write.table(lrt3.all, "lrtSuS.SeMvsF_gonad_samples_TMM.txt", sep="\t")
is.de.lrt3 <- decideTestsDGE(lrt3)
summary(is.de.lrt3)
plotMD(lrt3, status=is.de.lrt3)

lrt4 <- glmLRT(fit, contrast=my.contrasts[,"SuP.SeMvsF"])
lrt4.all <- topTags(lrt4, n=73220)
head(lrt4.all)
write.table(lrt4.all, "lrtSuP.SeMvsFL_gonad_samples_TMM.txt", sep="\t")
is.de.lrt4 <- decideTestsDGE(lrt4)
summary(is.de.lrt4)
plotMD(lrt4, status=is.de.lrt4)

lrt5 <- glmLRT(fit, contrast=my.contrasts[,"Full"])
lrt5.all <- topTags(lrt5, n=73220)
head(lrt5.all)
write.table(lrt1.all, "lrtFull_gonad_samples_TMM.txt", sep="\t")
is.de.lrt5 <- decideTestsDGE(lrt5)
summary(is.de.lrt5)
plotMD(lrt5, status=is.de.lrt5)

lrt6 <- glmLRT(fit, contrast=my.contrasts[,"MSvsFP"])
lrt6.gonad <- topTags(lrt6, n=73220)
head(lrt6.gonad)
write.table(lrt6.gonad, "lrtMSvsFP_gonad_samples_TMM.txt", sep="\t")
is.de.lrt6 <- decideTestsDGE(lrt6)
summary(is.de.lrt6)
plotMD(lrt6, status=is.de.lrt6)

lrt7 <- glmLRT(fit, contrast=my.contrasts[,"MPvsFS"])
lrt7.gonad <- topTags(lrt7, n=73220)
head(lrt7.gonad)
write.table(lrt7.gonad, "lrtMPvsFS_gonad_samples_TMM.txt", sep="\t")
is.de.lrt7 <- decideTestsDGE(lrt7)
summary(is.de.lrt7)
plotMD(lrt7, status=is.de.lrt7)

lrt8 <- glmLRT(fit, contrast=my.contrasts[,"Sex"])
lrt8.gonad <- topTags(lrt8, n=73220)
head(lrt8.gonad)
write.table(lrt8.gonad, "lrtSexContrast_gonad_samples_TMM.txt", sep="\t")
is.de.lrt8 <- decideTestsDGE(lrt8)
summary(is.de.lrt8)
plotMD(lrt8, status=is.de.lrt8)

lrt9 <- glmLRT(fit, contrast=my.contrasts[,"Substrate"])
lrt9.gonad <- topTags(lrt9, n=73220)
head(lrt9.gonad)
write.table(lrt9.gonad, "lrtSubstrateContrast_gonad_samples_TMM.txt", sep="\t")
is.de.lrt9 <- decideTestsDGE(lrt9)
summary(is.de.lrt9)
plotMD(lrt9, status=is.de.lrt9)



#### Now all of the gill samples

x <- read.delim("all_genes_gill_samples.txt", row.names="RefName")
head(x)

## Create the DGEList data class for analysis in edgeR
## Keeping it simple with just counts and adding groups later

y <- DGEList(counts=x)
plotMDS(y)


# reading in data to identify groups for the sample data
targets <- read.table("sample_info_gill.txt", header=TRUE)

# creating a model matrix for the GLM
design <- model.matrix(~Sex*Substrate, data=targets)
colnames(design)

design

##Filter low read count genes out of analysis
keep <- filterByExpr(y, design)
table(keep)

## Reset DGElist object to only have genes with sufficient read counts

y <- y[keep, , keep.lib.sizes=FALSE]
y$samples


AveLogCPM <- aveLogCPM(y)
hist(AveLogCPM)

## Normalize read counts

y <- calcNormFactors(y)
y$samples

y.counts <- y$counts

write.table(y.counts, "TMM_Gill_Samples_Counts.txt",quote=FALSE)

#estimating dispersions for the GLM model
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTrendedDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)

# GLM fit to the fit the negative binomial GLM for each tag
fit <- glmFit(y,design)
colnames(fit)

# glmLRT to carry out likelihood ratio test; 
# select one or more coefficients to drop from full model to serve as a null model
# against which the full model is compared in a likelihood ratio test


lrtSE <- glmLRT(fit, coef=2)
lrtSE.all <- topTags(lrtSE, n=66625)
head(lrtSE.all)
write.table(lrtSE.all, "lrtSE_gill_samples_only_TMM.txt", sep="\t")
is.de.lrtSE <- decideTestsDGE(lrtSE)
summary(is.de.lrtSE)
plotMD(lrtSE, status=is.de.lrtSE)

lrtSU <- glmLRT(fit, coef=3)
lrtSU.all <- topTags(lrtSU, n=66625)
head(lrtSU.all)
write.table(lrtSU.all, "lrtSU_gill_samples_only_TMM.txt", sep="\t")
is.de.lrtSU <- decideTestsDGE(lrtSU)
summary(is.de.lrtSU)
plotMD(lrtSU, status=is.de.lrtSU)


lrtSExSU <- glmLRT(fit, coef=4)
lrtSExSU.all <- topTags(lrtSExSU, n=66625)
head(lrtSExSU.all)
write.table(lrtSExSU.all, "lrtSExSU_gill_samples_only_TMM.txt", sep="\t")
is.de.lrtSExSU <- decideTestsDGE(lrtSExSU)
summary(is.de.lrtSExSU)
plotMD(lrtSExSU, status=is.de.lrtSExSU)



## interaction term analysis using contrasts
# first creating design object and groups
Group <- factor(paste(targets$Sex, targets$Substrate, sep="."))
cbind(targets, Group=Group)
design <- model.matrix(~0+Group)
colnames(design) <- levels(Group)

##Filter low read count genes out of analysis
keep <- filterByExpr(y, design)
table(keep)
## Reset DGElist object to only have genes with sufficient read counts

y <- y[keep, , keep.lib.sizes=FALSE]
y$samples


AveLogCPM <- aveLogCPM(y)
hist(AveLogCPM)

## Normalize read counts

y <- calcNormFactors(y)
y$samples

# next generation dispersions
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTrendedDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)

pch <- c(0,15,1,16)
colors <- rep(c("dodgerblue2","dodgerblue2","red","red"))
parameter<-plotMDS(y, col=colors[Group], pch=pch[Group],cex=2.1,cex.axis=1.4,cex.lab=1.4)
legend(0.1,0.2, legend=levels(Group), 
       xpd=TRUE,inset=c(0,0),
       pch=pch, col=colors, cex=1, ncol=1)

#PERMANOVA using adonis function in vegan package
y_2<-t(y.counts)
adonis2(y_2~Sex*Substrate, data=targets, permutations=999)

# GLM fit to the fit the negative binomial GLM for each tag
fit <- glmFit(y,design)
colnames(fit)

# creating contrasts for comparisons within the interaction term
my.contrasts <-makeContrasts(
  SeM.SuSvsP=Male.Shell-Male.Plastic,
  SeF.SuSvsP=Female.Shell-Female.Plastic,
  SuS.SeMvsF=Male.Shell-Female.Shell,
  SuP.SeMvsF=Male.Plastic-Female.Plastic,
  Full=(Male.Shell-Female.Shell)-(Male.Plastic-Female.Plastic),
  MSvsFP = Male.Shell-Female.Plastic,
  MPvsFS = Male.Plastic-Female.Shell,
  Sex=(Male.Shell-Female.Shell)-(Female.Plastic-Male.Plastic),
  Substrate=(Male.Shell-Male.Plastic)-(Female.Plastic-Female.Shell),
  levels=design)

lrt1 <- glmLRT(fit, contrast=my.contrasts[,"SeM.SuSvsP"])
lrt1.all <- topTags(lrt1, n=73220)
head(lrt1.all)
write.table(lrt1.all, "lrtSeM.SuSvsP_gill_samples_TMM.txt", sep="\t")
is.de.lrt1 <- decideTestsDGE(lrt1)
summary(is.de.lrt1)
plotMD(lrt1, status=is.de.lrt1)


lrt2 <- glmLRT(fit, contrast=my.contrasts[,"SeF.SuSvsP"])
lrt2.all <- topTags(lrt2, n=73220)
head(lrt2.all)
write.table(lrt2.all, "lrtSeF.SuSvsP_gill_samples_TMM.txt", sep="\t")
is.de.lrt2 <- decideTestsDGE(lrt2)
summary(is.de.lrt2)
plotMD(lrt2, status=is.de.lrt2)

lrt3 <- glmLRT(fit, contrast=my.contrasts[,"SuS.SeMvsF"])
lrt3.all <- topTags(lrt3, n=73220)
head(lrt3.all)
write.table(lrt3.all, "lrtSuS.SeMvsF_gill_samples_TMM.txt", sep="\t")
is.de.lrt3 <- decideTestsDGE(lrt3)
summary(is.de.lrt3)
plotMD(lrt3, status=is.de.lrt3)

lrt4 <- glmLRT(fit, contrast=my.contrasts[,"SuP.SeMvsF"])
lrt4.all <- topTags(lrt4, n=73220)
head(lrt4.all)
write.table(lrt4.all, "lrtSuP.SeMvsFL_gill_samples_TMM.txt", sep="\t")
is.de.lrt4 <- decideTestsDGE(lrt4)
summary(is.de.lrt4)
plotMD(lrt4, status=is.de.lrt4)

lrt5 <- glmLRT(fit, contrast=my.contrasts[,"Full"])
lrt5.all <- topTags(lrt5, n=73220)
head(lrt5.all)
write.table(lrt5.all, "lrtFull_gill_samples_TMM.txt", sep="\t")
is.de.lrt5 <- decideTestsDGE(lrt5)
summary(is.de.lrt5)
plotMD(lrt5, status=is.de.lrt5)

lrt6 <- glmLRT(fit, contrast=my.contrasts[,"MSvsFP"])
lrt6.gill <- topTags(lrt6, n=73220)
head(lrt6.gill)
write.table(lrt6.gill, "lrtMSvsFP_gill_samples_TMM.txt", sep="\t")
is.de.lrt6 <- decideTestsDGE(lrt6)
summary(is.de.lrt6)
plotMD(lrt6, status=is.de.lrt6)

lrt7 <- glmLRT(fit, contrast=my.contrasts[,"MPvsFS"])
lrt7.gill <- topTags(lrt7, n=73220)
head(lrt7.gill)
write.table(lrt7.gill, "lrtMPvsFS_gill_samples_TMM.txt", sep="\t")
is.de.lrt7 <- decideTestsDGE(lrt7)
summary(is.de.lrt7)
plotMD(lrt7, status=is.de.lrt7)

lrt8 <- glmLRT(fit, contrast=my.contrasts[,"Sex"])
lrt8.gill <- topTags(lrt8, n=73220)
head(lrt8.gill)
write.table(lrt8.gill, "lrtSexContrast_gill_samples_TMM.txt", sep="\t")
is.de.lrt8 <- decideTestsDGE(lrt8)
summary(is.de.lrt8)
plotMD(lrt8, status=is.de.lrt8)

lrt9 <- glmLRT(fit, contrast=my.contrasts[,"Substrate"])
lrt9.gill <- topTags(lrt9, n=73220)
head(lrt9.gill)
write.table(lrt9.gill, "lrtSubstrateContrast_gill_samples_TMM.txt", sep="\t")
is.de.lrt9 <- decideTestsDGE(lrt9)
summary(is.de.lrt9)
plotMD(lrt9, status=is.de.lrt9)


### Gonad minus male P8

x <- read.delim("all_genes_gonad_samples-P8.txt", row.names="RefName")

head(x)

## Create the DGEList data class for analysis in edgeR
## Keeping it simple with just counts and adding groups later

y <- DGEList(counts=x)
y$samples


# reading in data to identify groups for the sample data
targets <- read.table("sample_info_gonad-P8.txt", header=TRUE)

# creating a model matrix for the GLM
design <- model.matrix(~Sex*Substrate, data=targets)
colnames(design)

design

##Filter low read count genes out of analysis
keep <- filterByExpr(y, design)
table(keep)

## Reset DGElist object to only have genes with sufficient read counts

y <- y[keep, , keep.lib.sizes=FALSE]
y$samples


AveLogCPM <- aveLogCPM(y)
hist(AveLogCPM)

## Normalize read counts

y <- calcNormFactors(y)
y$samples


#estimating dispersions for the GLM model
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTrendedDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)

y.counts<-(y$counts)

# GLM fit to the fit the negative binomial GLM for each tag

fit <- glmFit(y,design)
colnames(fit)

# glmLRT to carry out likelihood ratio test; 
# select one or more coefficients to drop from full model to serve as a null model
# against which the full model is compared in a likelihood ratio test


lrtSE <- glmLRT(fit, coef=2)
lrtSE.all <- topTags(lrtSE, n=66625)
head(lrtSE.all)
write.table(lrtSE.all, "lrtSE_Se*SuModel_gonad_samples_only_TMM-P8.txt", sep="\t")
is.de.lrtSE <- decideTestsDGE(lrtSE)
summary(is.de.lrtSE)
plotMD(lrtSE, status=is.de.lrtSE)

lrtSU <- glmLRT(fit, coef=3)
lrtSU.all <- topTags(lrtSU, n=66625)
head(lrtSU.all)
write.table(lrtSU.all, "lrtSU_Se*SuModel_gonad_samples_only_TMM-P8.txt", sep="\t")
is.de.lrtSU <- decideTestsDGE(lrtSU)
summary(is.de.lrtSU)
plotMD(lrtSU, status=is.de.lrtSU)

lrtSExSU <- glmLRT(fit, coef=4)
lrtSExSU.all <- topTags(lrtSExSU, n=66625)
head(lrtSExSU.all)
write.table(lrtSExSU.all, "lrtSExSU_Se*SuModel_gonad_samples_only_TMM-P8.txt", sep="\t")
is.de.lrtSExSU <- decideTestsDGE(lrtSExSU)
summary(is.de.lrtSExSU)
plotMD(lrtSExSU, status=is.de.lrtSExSU)

## interaction term analysis using contrasts
# first creating design object and groups
Group <- factor(paste(targets$Sex, targets$Substrate, sep="."))
cbind(targets, Group=Group)
design <- model.matrix(~0+Group)
colnames(design) <- levels(Group)
table(Group)



##Filter low read count genes out of analysis

keep <- filterByExpr(y, design)
table(keep)

## Reset DGElist object to only have genes with sufficient read counts

y <- y[keep, , keep.lib.sizes=FALSE]
y$samples


AveLogCPM <- aveLogCPM(y)
hist(AveLogCPM)

## Normalize read counts

y <- calcNormFactors(y)
y$samples

# next generation dispersions
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTrendedDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)

plotMDS(y)
pch <- c(0,15,1,16)
colors <- rep(c("dodgerblue2","dodgerblue2","red","red"))
parameter<-plotMDS(y, col=colors[Group], pch=pch[Group],cex=2.1,cex.axis=1.4,cex.lab=1.4)
legend("topright", legend=levels(Group), pch=pch, col=colors, cex=1.2, ncol=1)

y_2<-t(y.counts)
adonis2(y_2~Sex*Substrate, data=targets, permutations=999)

# GLM fit to the fit the negative binomial GLM for each tag
fit <- glmFit(y,design)
colnames(fit)

# creating contrasts for comparisons within the interaction term

my.contrasts <-makeContrasts(
  SeM.SuSvsP=Male.Shell-Male.Plastic,
  SeF.SuSvsP=Female.Shell-Female.Plastic,
  SuS.SeMvsF=Male.Shell-Female.Shell,
  SuP.SeMvsF=Male.Plastic-Female.Plastic,
  Full=(Male.Shell-Female.Shell)-(Male.Plastic-Female.Plastic),
  MSvsFP = Male.Shell-Female.Plastic,
  MPvsFS = Male.Plastic-Female.Shell,
  Sex=(Male.Shell-Female.Shell)-(Female.Plastic-Male.Plastic),
  Substrate=(Male.Shell-Male.Plastic)-(Female.Plastic-Female.Shell),
  levels=design)

lrt1 <- glmLRT(fit, contrast=my.contrasts[,"SeM.SuSvsP"])
lrt1.all <- topTags(lrt1, n=73220)
head(lrt1.all)
write.table(lrt1.all, "lrtSeM.SuSvsP_gonad_samples_TMM-P8.txt", sep="\t")
is.de.lrt1 <- decideTestsDGE(lrt1)
summary(is.de.lrt1)
plotMD(lrt1, status=is.de.lrt1)

lrt2 <- glmLRT(fit, contrast=my.contrasts[,"SeF.SuSvsP"])
lrt2.all <- topTags(lrt2, n=73220)
head(lrt2.all)
write.table(lrt2.all, "lrtSeF.SuSvsP_gonad_samples_TMM-P8.txt", sep="\t")
is.de.lrt2 <- decideTestsDGE(lrt2)
summary(is.de.lrt2)
plotMD(lrt2, status=is.de.lrt2)

lrt3 <- glmLRT(fit, contrast=my.contrasts[,"SuS.SeMvsF"])
lrt3.all <- topTags(lrt3, n=73220)
head(lrt3.all)
write.table(lrt3.all, "lrtSuS.SeMvsF_gonad_samples_TMM-P8.txt", sep="\t")
is.de.lrt3 <- decideTestsDGE(lrt3)
summary(is.de.lrt3)
plotMD(lrt3, status=is.de.lrt3)

lrt4 <- glmLRT(fit, contrast=my.contrasts[,"SuP.SeMvsF"])
lrt4.all <- topTags(lrt4, n=73220)
head(lrt4.all)
write.table(lrt4.all, "lrtSuP.SeMvsFL_gonad_samples_TMM-P8.txt", sep="\t")
is.de.lrt4 <- decideTestsDGE(lrt4)
summary(is.de.lrt4)
plotMD(lrt4, status=is.de.lrt4)

lrt5 <- glmLRT(fit, contrast=my.contrasts[,"Full"])
lrt5.all <- topTags(lrt5, n=73220)
head(lrt5.all)
write.table(lrt1.all, "lrtFull_gonad_samples_TMM-P8.txt", sep="\t")
is.de.lrt5 <- decideTestsDGE(lrt5)
summary(is.de.lrt5)
plotMD(lrt5, status=is.de.lrt5)

lrt6 <- glmLRT(fit, contrast=my.contrasts[,"MSvsFP"])
lrt6.gonad <- topTags(lrt6, n=73220)
head(lrt6.gonad)
write.table(lrt6.gonad, "lrtMSvsFP_gonad_samples_TMM-P8.txt", sep="\t")
is.de.lrt6 <- decideTestsDGE(lrt6)
summary(is.de.lrt6)
plotMD(lrt6, status=is.de.lrt6)

lrt7 <- glmLRT(fit, contrast=my.contrasts[,"MPvsFS"])
lrt7.gonad <- topTags(lrt7, n=73220)
head(lrt7.gonad)
write.table(lrt7.gonad, "lrtMPvsFS_gonad_samples_TMM-P8.txt", sep="\t")
is.de.lrt7 <- decideTestsDGE(lrt7)
summary(is.de.lrt7)
plotMD(lrt7, status=is.de.lrt7)

lrt8 <- glmLRT(fit, contrast=my.contrasts[,"Sex"])
lrt8.gonad <- topTags(lrt8, n=73220)
head(lrt8.gonad)
write.table(lrt8.gonad, "lrtSexContrast_gonad_samples_TMM-P8.txt", sep="\t")
is.de.lrt8 <- decideTestsDGE(lrt8)
summary(is.de.lrt8)
plotMD(lrt8, status=is.de.lrt8)

lrt9 <- glmLRT(fit, contrast=my.contrasts[,"Substrate"])
lrt9.gonad <- topTags(lrt9, n=73220)
head(lrt9.gonad)
write.table(lrt9.gonad, "lrtSubstrateContrast_gonad_samples_TMM-P8.txt", sep="\t")
is.de.lrt9 <- decideTestsDGE(lrt9)
summary(is.de.lrt9)
plotMD(lrt9, status=is.de.lrt9)
