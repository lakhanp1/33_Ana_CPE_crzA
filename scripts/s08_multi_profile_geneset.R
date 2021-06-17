suppressPackageStartupMessages(library(chipmine))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(org.AFumigatus.A1163.eg.db))
suppressPackageStartupMessages(library(TxDb.Afumigatus.A1163.AspGD.GFF))


## This script plots the profile heatmaps for multiple samples. If the expression values are present, it
## also plots a simple heatmap using these expression values

rm(list = ls())

##################################################################################

## IMP: the first sampleID will be treated primary and clustering will be done/used for/of this sample
comparisonName <- "figure8"
outDir <- here::here("analysis", "06_crzA_binding_profile")
outPrefix <- paste(outDir, "/", comparisonName, ".geneset_profile", sep = "")

file_plotSamples <- paste(outDir, "/", "samples.txt", sep = "")
file_geneset <- paste(outDir, "/", "geneset1.txt", sep = "")

matrixType <- "2kb_TSS_1kb"
up <- 2000
body <- 0
down <- 1000
binSize <- 10

matrixDim = c(c(up, body, down)/binSize, binSize)

showExpressionHeatmap <- FALSE

## genes to read
file_exptInfo <- here::here("data", "reference_data", "ChIPseq_sample_info.tab")

TF_dataPath <- here::here("data", "ChIPseq_CEA17")
polII_dataPath <- here::here("data", "polII_data")
hist_dataPath <- here::here("data", "histone_data")
other_dataPath <- here::here("data", "other_data")

orgDb <- org.AFumigatus.A1163.eg.db
txDb <- TxDb.Afumigatus.A1163.AspGD.GFF

## colors
colList <- list()


##################################################################################
sampleList <- suppressMessages(readr::read_tsv(file = file_plotSamples, comment = "#"))

tempSInfo <- get_sample_information(
  exptInfoFile = file_exptInfo, samples = sampleList$sampleId,
  dataPath = TF_dataPath, profileMatrixSuffix = matrixType
)

tfIds <- tempSInfo$sampleId[which(tempSInfo$IP_tag %in% c("HA", "MYC", "TAP", "TF"))]
inputIds <- tempSInfo$sampleId[which(tempSInfo$IP_tag %in% c("input"))]
polII_ids <- tempSInfo$sampleId[which(tempSInfo$IP_tag == "polII")]
histIds <- tempSInfo$sampleId[which(tempSInfo$IP_tag == "HIST")]


## read the experiment sample details and select only those which are to be plotted
tfData <- get_sample_information(
  exptInfoFile = file_exptInfo, samples = tfIds,
  dataPath = TF_dataPath, profileMatrixSuffix = matrixType
)


inputData <- get_sample_information(
  exptInfoFile = file_exptInfo, samples = inputIds,
  dataPath = TF_dataPath, profileMatrixSuffix = matrixType
)

polIIData <- get_sample_information(
  exptInfoFile = file_exptInfo, samples = polII_ids,
  dataPath = polII_dataPath, profileMatrixSuffix = matrixType
)

histData <- get_sample_information(
  exptInfoFile = file_exptInfo, samples = histIds,
  dataPath = hist_dataPath, profileMatrixSuffix = matrixType
)


exptData <- dplyr::bind_rows(tfData, inputData, histData, polIIData)

exptDataList <- purrr::transpose(exptData)  %>% 
  purrr::set_names(nm = purrr::map(., "sampleId"))

polIICols <- list(
  exp = structure(polII_ids, names = polII_ids),
  is_expressed = structure(paste("is_expressed", ".", polII_ids, sep = ""), names = polII_ids)
)


tfCols <- sapply(
  X = c("hasPeak", "peakId", "peakEnrichment", "peakPval", "peakQval", "peakSummit", "peakDist", "summitDist",
        "peakType", "bidirectional", "featureCovFrac", "relativeSummitPos", "peakRegion", "peakPosition",
        "peakCoverage", "pvalFiltered", "summitSeq"),
  FUN = function(x){ structure(paste(x, ".", tfIds, sep = ""), names = tfIds) },
  simplify = F, USE.NAMES = T)

##################################################################################

# ## genes to read
# geneSet <- data.table::fread(file = file_genes, header = F,
#                              col.names = c("chr", "start", "end", "geneId", "score", "strand"))
# 
# kmClust <- dplyr::left_join(
#   x = suppressMessages(readr::read_tsv(file = tfData$clusterFile[1])),
#   y = geneSet, by = c("geneId" = "geneId")
# )
# 
# geneDesc <- AnnotationDbi::select(x = orgDb, keys = geneSet$gene, keytype = "GID",
#                                   columns = c("GENE_NAME", "DESCRIPTION")) %>% 
#   dplyr::rename(geneId = GID)
# 
# # geneInfo <- dplyr::left_join(x = kmClust, y = geneDesc, by = c("geneId" = "GID"))
# geneInfo <- geneDesc
# 
# head(geneInfo)
# 
# 
# expressionData <- get_TF_binding_data(exptInfo = tfData,
#                                       genesDf = geneInfo)
# 
# # expressionData <- get_polII_expressions(exptInfo = polIIData,
# #                                         genesDf = expressionData)
# 
# # view(dfSummary(expressionData))
# 
# peakTargetMat <- peak_target_matrix(sampleInfo = tfData, position = "best")


anLables <- list()
# anLables[[tssPeakTypeCol]] = gsub("peakType", "TSS peak type\n", tssPeakTypeCol) %>% gsub("\\(|\\)", "", .)
# anLables[[tesPeakTypeCol]] = gsub("tesPeakType", "TES peak type\n", tesPeakTypeCol) %>% gsub("\\(|\\)", "", .)
# anLables[[isExpCol]] = txt = gsub("is_expressed", "is expressed\n", isExpCol) %>% gsub("\\(|\\)", "", .)
anLables[["is_SM_gene"]] = "SM gene"
anLables[["is_TF"]] = "Transcription Factor"
anLables[["gene_length"]] = "Gene Length"

##################################################################################


## generate profile matrix for the first time
# genesGr <- rtracklayer::import(con = file_genes, format = "bed")
# 
# geneStartGr <- GenomicRanges::resize(x = genesGr, width = 1, fix = "start")
# 
# i <- 1
# for(i in 1:nrow(exptData)){
#   bwMat <- bigwig_profile_matrix(bwFile = exptData$bwFile[i],
#                                  regions = geneStartGr,
#                                  signalName = exptData$sampleId[i],
#                                  extend = c(up, down),
#                                  targetName = "ATG",
#                                  storeLocal = TRUE,
#                                  localPath = exptData$matFile[i])
# }

## color list
matList <- import_profiles(
  exptInfo = dplyr::bind_rows(tfData, inputData, histData, polIIData),
  # geneList = geneInfo$geneId,
  targetType = "point", targetName = "TSS",
  up = matrixDim[1], target = matrixDim[2], down = matrixDim[3]
)


## tf colors
tfMeanProfile <- NULL
if(length(c(tfIds)) == 1){
  tfMeanProfile <- matList[[tfIds]]
} else{
  tfMeanProfile <- getSignalsFromList(lt = matList[tfIds])
}

quantile(tfMeanProfile, c(seq(0, 0.9, by = 0.1), 0.95, 0.99, 0.992, 0.995, 0.999, 0.9999, 1), na.rm = T)
# tfMeanColor <- colorRamp2(quantile(tfMeanProfile, c(0.50, 0.995), na.rm = T), c("white", "red"))
tfColorList <- sapply(
  X = c(tfIds, inputIds),
  FUN = function(x){
    return(
      colorRamp2(
        # breaks = quantile(tfMeanProfile, c(0.50, 0.95), na.rm = T),
        # colors = unlist(strsplit(x = exptDataList[[x]]$color, split = ","))
        breaks = quantile(tfMeanProfile, c(0.6, 0.7, 0.995, 0.999), na.rm = T),
        colors = c("white", "#ffffcc", "#e31a1c", "#bd0026")
      )
    )
  }
)

## polII colors
polIIMeanProfile <- NULL
polIIColorList <- NULL
# if(nrow(polIIData) == 1){
#   polIIMeanProfile <- matList[[polIIData$sampleId]]
# } else{
#   polIIMeanProfile <- getSignalsFromList(lt = matList[polIIData$sampleId])
# }
# quantile(polIIMeanProfile, c(seq(0, 0.9, by = 0.1), 0.95, 0.99, 0.992, 0.995, 0.999, 0.9999, 1), na.rm = T)
# polIIMeanColor <- colorRamp2(quantile(polIIMeanProfile, c(0.01, 0.5, 0.995), na.rm = T), c("blue", "white", "red"))
# polIIColorList <- sapply(X = polIIData$sampleId, FUN = function(x){return(polIIMeanColor)})

## histone colors
histMeanProfile <- NULL
histColorList <- NULL
# if(nrow(histData) == 1){
#   histMeanProfile <- matList[[histData$sampleId]]
# } else{
#   histMeanProfile <- getSignalsFromList(lt = matList[histData$sampleId])
# }
# quantile(histMeanProfile, c(seq(0, 0.9, by = 0.1), 0.95, 0.99, 0.992, 0.995, 0.999, 0.9999, 1), na.rm = T)
# # histMeanColor <- colorRamp2(quantile(histMeanProfile, c(0.30, 0.995), na.rm = T), c("black", "yellow"))
# histColorList <- sapply(
#   X = histData$sampleId,
#   FUN = function(x){
#     return(colorRamp2(breaks = quantile(histMeanProfile, c(0.30, 0.995), na.rm = T),
#                       colors =  unlist(strsplit(x = exptDataList[[x]]$color, split = ","))))
#   }
# )

colorList <- unlist(list(tfColorList, polIIColorList, histColorList))
ylimList <- list()
# ylimList <- sapply(c(polII_ids, histIds), function(x){return(0.996)}, simplify = FALSE)
ylimList <- append(x = ylimList,
                   values = sapply(c(tfIds, inputIds), function(x){return(c(0, 25))}, simplify = FALSE))



##################################################################################
# plot profiles for genes of interest

geneSubset <- suppressMessages(
  readr::read_tsv(file = file_geneset)
) %>% 
  dplyr::mutate(
    category = forcats::as_factor(category)
  )

# geneSubset <- dplyr::left_join(x = geneSubset, y = peakTargetMat, by = "geneId") %>% 
#   dplyr::left_join(y = geneInfo, by = "geneId")


multiProfiles_geneset <- multi_profile_plots(
  exptInfo = exptData,
  genesToPlot = geneSubset$geneId,
  targetType = "point",
  targetName = "ATG",
  matBins = matrixDim,
  clusters = dplyr::select(geneSubset, geneId, cluster = category),
  showAnnotation = FALSE,
  drawClusterAn = FALSE,
  profileColors = colorList,
  column_title_gp = gpar(fontsize = 12),
  # row_order = geneSubset$geneId,
  show_row_names = TRUE,
  row_labels = geneSubset$geneName,
  posLineGpar = gpar(col = "black", lty = 2, alpha = 1, lwd = 2),
  ylimFraction = ylimList
)


pdfWd <- 2 + 
  (length(multiProfiles_geneset$heatmapList@ht_list) * 2) +
  (length(polII_ids) * 0.25 * showExpressionHeatmap)


# ## gene length annotation
# anGl_geneset <- gene_length_heatmap_annotation(
#   bedFile = file_genes,
#   genes = geneSubset$geneId,
#   axis_param = list(at = c(2000, 4000), labels = c("2kb", "> 4kb")),
#   pointSize = unit(4, "mm"))


geneset_htlist <- multiProfiles_geneset$heatmapList



# wd <- 500 + (nrow(exptData) * 700) + (length(polII_ids) * 500 * showExpressionHeatmap)
title_geneset = paste(comparisonName, ": genes of interest", collapse = "")

# draw Heatmap and add the annotation name decoration
pdf(file = paste(outPrefix, ".pdf", sep = ""), width = 10, height = 10)

geneset_htlist <- draw(
  geneset_htlist,
  main_heatmap = exptData$profileName[1],
  # annotation_legend_list = list(profile1$legend),
  column_title = title_geneset,
  column_title_gp = gpar(fontsize = 14, fontface = "bold"),
  row_sub_title_side = "left",
  heatmap_legend_side = "right",
  gap = unit(7, "mm"),
  padding = unit(rep(0.5, times = 4), "cm")
)

dev.off()


# rowOrderDf <- row_order(geneset_htlist) %>% 
#   purrr::map_dfr(.f = ~ tibble(rank = ., geneId = geneSubset$geneId[.]))
# 
# ordered_data <- dplyr::left_join(x = rowOrderDf, y = geneSubset, by = "geneId")
# 
# readr::write_tsv(x = ordered_data, file = paste(outPrefix_geneset, ".data.tab", sep = ""))

##################################################################################











