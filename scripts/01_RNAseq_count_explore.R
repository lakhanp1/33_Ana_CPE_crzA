# library(ballgown)
library(here)
library(tidyverse)
library(FactoMineR)
library(factoextra)
library(corrplot)
library(ggrepel)
library(GGally)
library(data.table)
library(DESeq2)
library(tximport)
library(matrixStats)

## This script uses ballgown to extract the FPKM matrix from the stringTie output
## It also perform PCA analysis using FPKM values
rm(list = ls())

analysisName <- "count_data_Af293"

readLength <- 100
file_sampleInfo <- here::here("data", "sample_info.txt")

outDir <- here("analysis", "01_RNAseq_Af293", analysisName)
outPrefix <- paste(outDir, "/", analysisName, sep = "")
# orgDb <- org.AFumigatus.Af293.eg.db

if(!dir.exists(outDir)){
  dir.create(path = outDir, recursive = T)
}


###########################################################################
exptInfo <- suppressMessages(readr::read_tsv(file = file_sampleInfo)) %>% 
  as.data.frame()

rownames(exptInfo) <- exptInfo$sampleId

## set the factor levels.
## control levels should be first
# exptInfo$gt <- factor(exptInfo$gt, levels = c("WT", "tet90"))
# exptInfo$treatment <- factor(exptInfo$treatment, levels = c("YPD", "YPD_DOX", "YPD_GdA"))
exptInfo$condition <- factor(exptInfo$condition, levels = unique(exptInfo$condition))

# ## Add gene symbol for each Ensembl ID
# geneInfo <- AnnotationDbi::select(x = orgDb,
#                                   keys = keys(x = orgDb, keytype = "GID"),
#                                   columns = c("GENE_NAME", "DESCRIPTION"),
#                                   keytype = "GID") %>% 
#   dplyr::rename(geneId = GID)

design <- ~ condition

###########################################################################
## import counts data: either by tximport or as raw count matrix

## import the counts data using tximport and run DESeq2
path_stringtie <- here::here("data", "stringTie_Af293")
filesStringtie <- paste(path_stringtie, "/stringTie_", exptInfo$sampleId, "/t_data.ctab", sep = "")
names(filesStringtie) <- exptInfo$sampleId

tmp <- data.table::fread(file = filesStringtie[1], sep = "\t", header = T, stringsAsFactors = F)
tx2gene <- tmp[, c("t_name", "gene_id")]

txi <- tximport(files = filesStringtie, type = "stringtie",
                tx2gene = tx2gene, readLength = readLength)

ddsTxi <- DESeqDataSetFromTximport(txi = txi, colData = exptInfo, design = design)
# assay(ddsTxi)
# colData(ddsTxi)
# rowData(ddsTxi)

## Run DESeq2
dds <- DESeq(ddsTxi)

# ## import raw counts data and run DESeq2
# file_rawCounts <- here::here("data", "MatrixCountsPerGeneBySample.Reneto.tab")
# 
# countsDf <- suppressMessages(readr::read_tsv(file = file_rawCounts, col_names = T)) %>%
#   as.data.frame()
# rownames(countsDf) <- countsDf$geneId
# countsDf$geneId <- NULL
# 
# 
# if(all(rownames(exptInfo) %in% colnames(countsDf))){
#   countsDf <- countsDf[, rownames(exptInfo)]
# } else{
#   stop("Column names in count matrix does not match with row names in experiment data")
# }
# 
# ## select only those sample rows which are part of current comparison
# # exptInfo <- droplevels(subset(exptInfo, condition %in% compare))
# # countsDf <- countsDf[, rownames(exptInfo)]
# 
# ## run DESeq2 and extract the processed data
# ddsCount <- DESeqDataSetFromMatrix(countData = countsDf, colData = exptInfo, design = design)
# 
# ## Run DESeq2
# dds <- DESeq(ddsCount)

###########################################################################
## raw counts
rawCounts <- tibble::rownames_to_column(as.data.frame(counts(dds, normalized = FALSE)), var = "geneId")

readr::write_tsv(x = rawCounts, path = paste(outPrefix, ".rawCounts.tab", sep = ""))


## FPKM
fpkmCounts <- tibble::rownames_to_column(as.data.frame(fpkm(dds)), var = "geneId")

# if(all(fpkmCounts$geneId %in% geneInfo$geneId)){
  fpkmCounts <- dplyr::left_join(x = fpkmCounts, y = geneInfo, by = c("geneId" = "geneId")) %>%
    dplyr::select(geneId, exptInfo$sampleId, everything()) %>%
    dplyr::filter(!is.na(geneId))
# }

fwrite(x = fpkmCounts, file = paste(outPrefix, ".FPKM.tab", sep = ""),
       sep = "\t", row.names = F, col.names = T, quote = F)

## normalized counts matrix
normCounts <- tibble::rownames_to_column(as.data.frame(counts(dds, normalized = TRUE)), var = "geneId")
readr::write_tsv(x = normCounts, path = paste0(c(outPrefix,".normCounts.tab"), collapse = ""))

## r-log normalized counts
rld <- rlog(dds, blind = FALSE)
rldCount <- rownames_to_column(as.data.frame(assay(rld)), var = "geneId")

readr::write_tsv(x = rldCount, path = paste(outPrefix, ".rlogCounts.tab", sep = ""))


pcaVars <- c("genotype", "treatment")
plotPCA(rld, intgroup = pcaVars, ntop = 4000)

pcaData <- plotPCA(rld, intgroup = pcaVars, returnData = TRUE, ntop = 4000)
percentVar <- sprintf("%.2f", 100 * attr(pcaData, "percentVar"))

pcaData <- dplyr::mutate_at(.tbl = pcaData, .vars = pcaVars, .funs = as.factor)

pltTitle <- "Principal Component Analysis"
pointCol <- base::structure(RColorBrewer::brewer.pal(n = length(unique(pcaData$genotype)), name = "Set1"),
                            names = levels(pcaData$genotype))


pcaPlot <- ggplot(pcaData, aes(x = PC1, y = PC2)) +
  geom_point(mapping = aes(color = genotype, shape = treatment), size=4) +
  geom_text_repel(mapping = aes(label = name), size = 3, point.padding = 0.5) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  scale_color_manual(values = pointCol) +
  ggtitle(pltTitle) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(face = "bold", size = 15),
        axis.title.y = element_text(face = "bold", size = 15),
        plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm"),
        legend.text = element_text(size = 13),
        legend.title = element_text(face = "bold", size = 15)
  )


png(filename = paste(outPrefix, ".PCA.png", sep = ""), width = 4000, height = 3000, res = 350)
pcaPlot
dev.off()


###########################################################################
## read the FPKM data from stringTie output using ballgown package
## no clearity on how ballgown calculate gene level FPKM scores using transcript score.
## so use stringTie -> tximport way
## tximport calculate gene level FPKM/TPM scores by using weighted average over all transcripts
## for a gene. See tximport code for the details
# filesStringtie <- paste(path_stringtie, "/stringTie_", exptInfo$sampleId, sep = "")
# 
# bg <- ballgown(samples = filesStringtie, meas='all')
# 
# geneExpression <- data.frame(gexpr(bg)) %>%
#   rownames_to_column(var = "geneId")
# 
# 
# ## rename the sample names
# renameCols <- base::structure(gsub(pattern = "FPKM.stringTie_", replacement = "", x = names(geneExpression)),
#                               names = names(geneExpression))
# 
# 
# renameCols[colnames(geneExpression)]
# 
# colnames(geneExpression) <- renameCols[colnames(geneExpression)]
# 
# 
# fwrite(x = geneExpression, file = paste(outPrefix, "FPKM_ballgown.tab", sep = "_"),
#        sep = "\t", col.names = T, row.names = F, quote = F)
# 
# 



###########################################################################
## PCA based on rld counts

normCountMat <- as.matrix(rldCount[, c(exptInfo$sampleId), drop = FALSE])
rownames(normCountMat) <- rldCount$geneId

## remove low count rows
keep <- rowSums(normCountMat > 1) >= 2
normCountMat <- normCountMat[keep, ]

## transform the data such that the genes are columns and each sample is a row
## also append the additional information for each sample using left_join()
df2 <- as.data.frame(t(normCountMat)) %>% 
  tibble::rownames_to_column(var = "sampleId") %>% 
  dplyr::left_join(y = exptInfo, by = c("sampleId" = "sampleId")) %>% 
  dplyr::select(!!!colnames(exptInfo), dplyr::everything())

# df2 <- df2 %>%
#   dplyr::filter(!treatment %in% c("C")) %>% 
#   dplyr::filter(!genotype %in% c("CEA17_dzipD"))

row.names(df2) <- df2$sampleId

res.pca <- PCA(df2, graph = FALSE, scale.unit = TRUE,
               quali.sup = 1:ncol(exptInfo), ncp = 10)

eig.val <- get_eigenvalue(res.pca)

## scree plot: variance by PC
fviz_eig(res.pca, addlabels = TRUE, ylim = c(0, 50))

## Graph of individuals
ind <- get_pca_ind(res.pca)

fviz_pca_ind(res.pca,
             # col.ind = exprData$treatment,
             fill.ind = df2$condition,
             pointshape = 21,
             repel = TRUE,
             mean.point = FALSE,
             legend.title = "Study",
             pointsize = 3
)


## prepare the plot dataframe for ggplot
plotData <- as.data.frame(ind$coord) %>%
  tibble::rownames_to_column(var = "sampleId") %>%
  dplyr::left_join(y = exptInfo, by = c("sampleId" = "sampleId"))

pairs(x = plotData[, 2:6],
      pch = 19,  cex = 1,
      col = plotData$condition,
      lower.panel=NULL)


## set the factor levels
plotData <- dplyr::mutate_at(.tbl = plotData, .vars = pcaVars, .funs = as.factor)

pointCol <- base::structure(
  RColorBrewer::brewer.pal(n = length(unique(exptInfo$genotype)), name = "Set1"),
  names = unique(exptInfo$genotype))

pltTitle <- "Principal Component Analysis: all samples"


## decide which PCs to use for plotting
pcToPlot <- c(1, 2)
pcCols <- grep(pattern = "Dim.", x = colnames(plotData), value = T)[pcToPlot]
fillColumn <- "genotype"
shapeColumn <- "treatment"

pcaPlot <- ggplot(data = plotData,
                  mapping = aes(x = !!sym(pcCols[1]), y = !!sym(pcCols[2]), label = sampleId)) +
  geom_point(mapping = aes(color = !!sym(fillColumn), shape = !!sym(shapeColumn)),
             size = 4, stroke = 2) +
  scale_shape_manual(values = c("C" = 1, "02CAS" = 15, "2CAS" = 17)) +
  # guides(fill=guide_legend(override.aes=list(shape=21))) +
  scale_color_manual(values = pointCol) +
  geom_text_repel(size = 3, point.padding = 0.5) +
  geom_hline(yintercept = 0, linetype = 2, alpha = 0.5) + 
  geom_vline(xintercept = 0, linetype = 2, alpha = 0.5) +
  xlab( paste("PC",pcToPlot[1]," (", sprintf("%.2f", eig.val[pcToPlot[1], "variance.percent"]), "%)", sep = "") ) +
  ylab( paste("PC",pcToPlot[2]," (", sprintf("%.2f", eig.val[pcToPlot[2], "variance.percent"]), "%)", sep = "") ) +
  ggtitle(pltTitle) + 
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(face = "bold"),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 13, face = "bold"))

print(pcaPlot)

# pdf(file = paste(outPrefix, ".rld_PCA.pdf", sep = ""), width = 10, height = 10)
png(filename = paste(outPrefix, ".rld_PCA.png", sep = ""), width = 4000, height = 3000, res = 350)
print(pcaPlot)
dev.off()

#############################################################################
## correlation scatter plot
pt <- ggpairs(data = as.data.frame(normCountMat),
              upper = list(continuous = wrap("points", size = 0.1)),
              lower = list(continuous = wrap("cor", size = 8)),
              diag = list(continuous = "densityDiag")) +
  theme_bw() +
  theme(
    strip.text.y = element_text(size = 18, angle = 0, hjust = 0, face = "bold"),
    strip.text.x = element_text(size = 18, angle = 90, hjust = 0, face = "bold")
  )

png(filename = paste(outPrefix, ".scatter_matrix.png", sep = ""),
    width = 10000, height = 10000, res = 250)

pt
dev.off()

# ###########################################################################
# ## calculate mean FPKM
# 
# exptInfo <- read.table(file = sampleInfoFile, header = T, sep = "\t", row.names = "sampleId")
# 
# ## IF NEEDED: modify the row names and column names
# # colnames(rawData) <- sub("_WT", "_", colnames(rawData))
# # rownames(exptInfo) <- sub("_WT", "_", rownames(exptInfo))
# 
# ## select only those sample rows which are part of current comparison
# designInfo <- exptInfo
# # designInfo <- droplevels(subset(exptInfo, condition %in% compare))
# 
# if(! all( rownames(designInfo) %in% colnames(geneExpression) ) ){
#   stop("Column names in FPKM matrix does not match with row names in experiment data")
# }
# 
# 
# ## get the sample names for each condition under comparison
# grp1 <- sapply(rownames(designInfo[designInfo$condition %in% compare[1], ]), FUN = as.name, USE.NAMES = F, simplify = T)
# name1 <- paste(compare[1], "_meanFPKM", sep = "")
# grp1Len <- length(grp1)
# 
# grp2 <- sapply(rownames(designInfo[designInfo$condition %in% compare[2], ]), FUN = as.name, USE.NAMES = F, simplify = T)
# name2 <- paste(compare[2], "_meanFPKM", sep = "")
# grp2Len <- length(grp2)
# 
# ## can add multiple groups
# 
# ## calculate the average FPKM
# fpkmData <- geneExpression %>% 
#   dplyr::select(geneId, !!!c(grp1, grp2)) %>%
#   rowwise() %>%
#   mutate(!!name1 := sum(!!!grp1) / !!grp1Len,
#          !!name2 := sum(!!!grp2) / !!grp2Len)
# 
# 
# 
# ## write the data to file
# fwrite(x = fpkmData, file = "FPKM_matrix.tab", sep = "\t", quote = F, col.names = T)
# 
# 


