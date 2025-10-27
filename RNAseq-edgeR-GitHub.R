# Open libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(DESeq2)
  library(extrafont)
  library(GEOquery)
  library(EnsDb.Hsapiens.v86)
  library(edgeR)
  library(circlize)
  library(grid)
})

# Set directory
workingDir={workingDir}
setwd(workingDir)
source('RNAseq-Functions.R')

# Get count file
countMatrix <- read.table(paste0(workingDir, 'readcounts.txt'), fill=TRUE, header=TRUE) 
countMatrix <- countMatrix %>% column_to_rownames('Gene')

# Read meta data file
metaData <- read.csv('metaData.csv', row.names='X')

# RAW COUNTS ###################################################################
# Get gene symbols
geneSymbol <- select(EnsDb.Hsapiens.v86, 
                     keys=as.character(row.names(countMatrix)), 
                     keytype="GENEID", 
                     columns=c("GENEID", "GENENAME", 'TXBIOTYPE', 'UNIPROTID', 'UNIPROTDB'))

# Colnames
colnames(geneSymbol) <- c("Gene", "geneName", "geneType", "UniprotID", 'UniprotType')

# Filter
geneSymbol <- geneSymbol %>% 
  dplyr::filter(!geneType%in%c('misc_RNA','TEC','unprocessed_pseudogene','processed_pseudogene', 'lincRNA', 'IG_C_pseudogene', 'antisense', 'IG_V_gene')) %>%
  dplyr::filter(UniprotType=='SWISSPROT') %>%
  dplyr::select(c('Gene', 'geneName')) %>% 
  unique()

# Add to raw counts
rawCountsGeneSym <- left_join(countMatrix %>% rownames_to_column('Gene'), geneSymbol, by="Gene") %>% 
  dplyr::filter(!is.na(geneName)) %>% 
  dplyr::select(-c('Gene')) %>%
  as.data.frame()

# Remove duplicate genes
rawCountsGeneSym <- rawCountsGeneSym[!duplicated(rawCountsGeneSym[c('geneName')]), ] %>% as.data.frame()
rawCountsGeneSym <- rawCountsGeneSym %>% remove_rownames %>% column_to_rownames('geneName')

# EDGER ########################################################################
# DESEQ filter
# Filter low read counts
keep <- rowSums(rawCountsGeneSym>5) >= 10
rawCountsGeneSym <- rawCountsGeneSym[keep,]

# Filter data 
cpm_log <- cpm(rawCountsGeneSym, log = TRUE)
median_log2_cpm <- apply(cpm_log, 1, median)
hist(median_log2_cpm)
expr_cutoff <- -1
abline(v = expr_cutoff, col = "red", lwd = 3)
sum(median_log2_cpm > expr_cutoff)
rawCountsGeneSymClean <- rawCountsGeneSym[median_log2_cpm > expr_cutoff, ]
heatmap(cor(cpm_log))

# PCA
pca <- prcomp(t(cpm_log), scale. = TRUE)
plot(pca$x[, 1], pca$x[, 2], pch = ".", xlab = "PC1", ylab = "PC2")
text(pca$x[, 1], pca$x[, 2], labels = colnames(cpm_log))
summary(pca)

# Groups
group <- metaData$Group
group

# Control vs affected ###########################################################
# Filter for groups of interest
groupsofI <- c() # populate with patient IDs
group <- (metaData[groupsofI,])$Group
rawCountsGeneSymCleanFilt <- rawCountsGeneSymClean %>% dplyr::select(groupsofI)

# DGE list
y <- DGEList(counts = rawCountsGeneSymCleanFilt, group = group)
y

# Normalise
y <- calcNormFactors(y)
y$samples

# model variance
y <- estimateDisp(y)
sqrt(y$common.dispersion) 
plotBCV(y)

# Test
et <- exactTest(y)
results_edgeR <- topTags(et, n = nrow(rawCountsGeneSymCleanFilt), sort.by = "none")
head(results_edgeR$table)

# Sig genes
sum(results_edgeR$table$FDR < .1)
plotSmear(et, de.tags = rownames(results_edgeR)[results_edgeR$table$FDR < .1])
abline(h = c(-2, 2), col = "blue")
paste0(rownames(results_edgeR)[results_edgeR$table$FDR < .1], collapse=', ')

# Save
df <- results_edgeR[results_edgeR$table$FDR < .1,] %>% as.data.frame()

# Covariates
age <- (metaData[groupsofI,])$Age
sex <- (metaData[groupsofI,])$Sex

# Design matrix
y <- DGEList(rawCountsGeneSymCleanFilt)
y <- calcNormFactors(y)
design <- model.matrix(~group + age + sex)
design

# Test
y <- estimateDisp(y, design)
fit <- glmFit(y, design)
lrt <- glmLRT(fit, coef = 2)
results_edgeR_lrt <- topTags(lrt, n = nrow(lrt$table))

# One gene
boxplot(as.numeric(rawCountsGeneSymCleanFilt["GNRHR", ]) ~ group)

# Get results
affDF <- results_edgeR_lrt %>% as.data.frame() %>% dplyr::filter(FDR<0.05)

# Column
results_edgeR_lrt_plot <- results_edgeR_lrt %>% as.data.frame() %>% rownames_to_column('geneName')

# Edit colours and labels
results_edgeR_lrt_plot$Label <- ifelse((results_edgeR_lrt_plot$FDR<0.05)&(abs(results_edgeR_lrt_plot$logFC)>0.5), results_edgeR_lrt_plot$geneName, NA)
results_edgeR_lrt_plot$Colour <- ifelse((results_edgeR_lrt_plot$FDR<0.05)&((results_edgeR_lrt_plot$logFC)>0.5), 'Increased', 'None')
results_edgeR_lrt_plot$Colour <- ifelse((results_edgeR_lrt_plot$FDR<0.05)&((results_edgeR_lrt_plot$logFC)<(-0.5)), 'Decreased', results_edgeR_lrt_plot$Colour)
results_edgeR_lrt_plot <- results_edgeR_lrt_plot %>% dplyr::filter(!is.na(FDR))
nrow(results_edgeR_lrt_plot %>% dplyr::filter(FDR<0.05))

# Plot
ggplot(results_edgeR_lrt_plot, aes(x=logFC, y=-log10(FDR), label=Label, colour=Colour, size=abs(logFC)*(-log10(FDR))))+
  geom_hline(yintercept=-log10(0.05), colour='grey')+
  geom_vline(xintercept=c(0.5,-0.5), colour='grey')+
  geom_point(alpha=0.75)+
  ggrepel::geom_text_repel(family='Roboto', colour='black', size=4)+
  theme_classic()+
  theme(text=element_text(family='Roboto'))+
  scale_colour_manual(values=c('#90bdcf', '#75a450', '#ebebeb'))+
  guides(colour='none', size='none')+
  scale_x_continuous(limits=c(-5,5))+
  scale_size(range=c(1,10))

# GSEA #########################################################################
# Column
results_edgeR_lrt <- results_edgeR_lrt %>% as.data.frame() %>% rownames_to_column('geneName')

# Get significant genes
paste0((results_edgeR_lrt %>% dplyr::filter(FDR<0.05))$geneName, collapse=', ')

# Create gene list
resLFCList <- results_edgeR_lrt$logFC
names(resLFCList) <- results_edgeR_lrt$geneName

# Clean
resLFCList <- na.omit(resLFCList)
resLFCList <- sort(resLFCList, decreasing=T)
resLFCList <- resLFCList[!duplicated(names(resLFCList))]

# Get GO BP gmt file
GO_file <- "Human_GO_bp_with_GO_iea_symbol_2023.gmt"
GO=fgsea::gmtPathways(GO_file)

# Run FGSEA
set.seed(1)
resGSEA <- fgsea::fgsea(pathways=GO,
                        stats=resLFCList,
                        minSize=10, 
                        maxSize=200
)

# Get significant pathways
resGSEAFiltPADJ <- resGSEA %>%
  dplyr::filter(padj<=0.05) %>% 
  arrange(pval)

# Leading edge
resGSEAFiltPADJ$leadingEdge[[15]]

# Collapse
resGSEACollapse <- fgsea::collapsePathways(
  resGSEA,
  pathways=GO,
  stats=resLFCList,
  pval.threshold = 0.05,
  nperm = 10/0.05,
  gseaParam = 1
)

# Filter
mainPathways <- resGSEA[pathway %in% resGSEACollapse$mainPathways][order(-NES), pathway]
resGSEAFilt <- resGSEA %>% dplyr::filter(pathway %in% mainPathways)

# Enrichment labels
resGSEAFilt$Enrichment=ifelse((resGSEAFilt$NES > 0)&(resGSEAFilt$padj<0.05), "Upregulated", 'None')
resGSEAFilt$Enrichment=ifelse((resGSEAFilt$NES < 0)&(resGSEAFilt$padj<0.05), "Downregulated", resGSEAFilt$Enrichment)

# GSEA volcano plot
ggplot(resGSEAFilt, aes(x=NES, y=-log10(padj), size=size)) +
  geom_hline(yintercept=-log10(0.05), color="#ebebeb", linewidth=0.5)+
  geom_vline(xintercept=1, color="#ebebeb", linewidth=0.5)+
  geom_vline(xintercept=-1, color="#ebebeb", linewidth=0.5) +
  geom_point(aes(colour=Enrichment), alpha=0.6) +
  scale_size_continuous(range=c(1,12))+
  theme_bw() +
  ggrepel::geom_text_repel(aes(label=stringr::str_wrap(pathway, 50)), family='Roboto', 
                           data=subset(resGSEAFilt, NES>1&padj<0.05),
                           nudge_y=-0.1,min.segment.length = 5,
                           lineheight=0.8,
                           size=4) +
  ggrepel::geom_text_repel(aes(label=stringr::str_wrap(pathway, 50)), family='Roboto', 
                           data=subset(resGSEAFilt, NES<(-1.5)&padj<0.05),
                           nudge_y=-0.075, min.segment.length = 5,
                           lineheight=0.8,
                           size=4) +
  xlab("NES") + ylab("-Log10 adjusted p-value") +
  theme(#legend.title=element_blank(),
    text=element_text(size=16, family='Roboto'),
    strip.background=element_rect(colour=NA, fill='white'),
    panel.grid=element_blank()) +
  scale_color_manual(values=c("#9c759b",'#ebebeb', "#75a450"))+
  scale_x_continuous(limits=c(-3,3))+
  guides(size='none',
         colour='none')

# Bar plot
ggplot(rbind(head(resGSEAFilt %>% dplyr::filter(NES>0, padj<0.05) %>% arrange(padj), n=10), head(resGSEAFilt %>% dplyr::filter(NES<0,padj<0.05) %>% arrange(padj), n=10)),
       aes(x=NES, y=fct_reorder(stringr::str_wrap(pathway, 40), NES), fill=NES))+
  geom_col(colour='black')+
  theme_classic()+
  theme(text=element_text(family='Roboto'))+
  labs(fill='-log10(Padj)', y='', x='NES')+
  scale_fill_gradientn(colors=c("#90bdcf",'#ebebeb', "#75a450"),
                       guide = guide_colorbar(frame.colour = "black", ticks.colour = "black"))+
  scale_x_continuous(expand=c(0,0))+
  guides(fill='none')
