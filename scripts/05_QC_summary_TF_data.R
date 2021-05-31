suppressPackageStartupMessages(library(chipmine))
suppressPackageStartupMessages(library(org.AFumigatus.A1163.eg.db))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(ggbeeswarm))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(ggrepel))

## 1) peak enrichment distribution
## 2) peak p-value distribution
## 3) peak annotation pie chart
## 4) combined matrix of peak enrichment distribution plots
## 5) combined matrix of peak p-value distribution plots

rm(list = ls())

source(file = "E:/Chris_UM/GitHub/omics_util/04_GO_enrichment/s01_enrichment_functions.R")

##################################################################################
analysisName <- "TF_ChIP_summary"
outDir <- here::here("analysis", "03_QC_ChIPseq")
outPrefix <- paste(outDir, "/", analysisName, sep = "")

file_exptInfo <- here::here("data", "reference_data", "ChIPseq_sample_info.tab")

orgDb <- org.AFumigatus.A1163.eg.db
txDb <- get_TxDb_sqlite(org = "A_fumigatus_A1163", yaml = "E:/Chris_UM/Database/TxDb.config.yaml")

TF_dataPath <- here::here("data", "ChIPseq_CEA17")
file_tf_macs2 <- paste(TF_dataPath, "/", "sample_tf_macs2.list", sep = "")

matrixType <- "2kb_summit"
up <- 2000
down <- 2000
body <- 0
bin <- 10
matrixDim = c(c(up, body, down)/bin, bin)



##################################################################################

tfSampleList <- suppressMessages(readr::read_tsv(file = file_tf_macs2))

tfInfo <- get_sample_information(
  exptInfoFile = file_exptInfo,
  samples = tfSampleList$sampleId,
  dataPath = TF_dataPath,
  profileMatrixSuffix = matrixType)

tfInfoList <- purrr::transpose(tfInfo)  %>% 
  purrr::set_names(nm = purrr::map(., "sampleId"))


allPlotData <- NULL

peakCountsDf <- tibble::tibble(
  sampleId = character(), peaks_total = numeric(), peaks_pval20 = numeric(),
  peaks_fe3 = numeric(), peaks_pval20_fe3 = numeric()
)


pdf(file = paste(outPrefix, ".macs2.pdf", sep = ""), width = 15, height = 10,
    onefile = TRUE, pointsize = 10)

rowId <- 1

for (rowId in 1:nrow(tfInfo)) {
  
  cat(rowId, ":", tfInfo$sampleId[rowId], "\n")
  
  peakType <- tfInfo$peakType[rowId]
  
  ## extract the peak counts to decide best TF replicate
  peaksGr <- rtracklayer::import(con = tfInfo$peakFile[rowId], format = peakType)
  peakCountsDf <- dplyr::bind_rows(
    peakCountsDf,
    tibble(
      sampleId = tfInfo$sampleId[rowId],
      peaks_total = length(peaksGr),
      peaks_pval20 = length(which(peaksGr$pValue >= 20)),
      peaks_fe3 = length(which(peaksGr$signalValue > 3)),
      peaks_pval20_fe3 = length(which(peaksGr$pValue >= 20 | peaksGr$signalValue > 3))
    )
  )

  chipSummary <- chip_summary(
    sampleId = tfInfo$sampleId[rowId],
    peakAnnotation = tfInfo$peakAnno[rowId],
    peakFile = tfInfo$peakFile[rowId],
    peakFormat = peakType
  )


  plot(chipSummary$figure)

  allPlotData <- dplyr::bind_rows(allPlotData, chipSummary$data)
  
}

dev.off()

bestReplicate <- dplyr::left_join(x = tfInfo, y = peakCountsDf, by = "sampleId") %>% 
  dplyr::select(!!!colnames(peakCountsDf), sampleName, condition, rep) %>% 
  dplyr::group_by(condition) %>% 
  dplyr::arrange(desc(peaks_pval20), .by_group = TRUE) %>% 
  dplyr::mutate(
    totalReps = n(),
    bestRep = row_number()
  ) %>% 
  dplyr::ungroup() %>% 
  dplyr::select(sampleId, starts_with("peaks_"), bestRep, totalReps, everything())

readr::write_tsv(x = bestReplicate, file = paste(outPrefix, ".best_replicates.tab", sep = ""))

##################################################################################
## combined summary plot matrix

theme_plot_matrix <- theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.title = element_text(size = 14),
    axis.text.y = element_text(size = 16),
    axis.text.x = element_blank(),
    panel.grid = element_blank(),
    title = element_text(hjust = 0.5, size = 20),
    strip.text = element_text(size = 8, hjust = 0),
    strip.background = element_rect(fill = "white")
  )

facetCols <- 2

## combined summary plot matrix: peak enrichment
gg_all_enrichment <- ggplot(
  data = allPlotData,
  mapping = aes(x = sampleId, y = peakEnrichment)) +
  geom_quasirandom(color = "#bf9a2d") +
  geom_boxplot(width=0.1, fill = NA, outlier.colour = NA, color = alpha("black", 1)) +
  geom_hline(yintercept = 3, color = "blue") +
  labs(title = "masc2 fold enrichment distribution") +
  facet_wrap(facets = ~ sampleId, ncol = facetCols, scales = "free") +
  theme_plot_matrix

# pdf(file = paste(outPrefix, ".macs2_enrichment.pdf", sep = ""), width = 16, height = 16, onefile = TRUE)
png(filename = paste(outPrefix, ".macs2_enrichment.png", sep = ""), width = 3000, height = 3000, res = 350)
print(gg_all_enrichment)
dev.off()

## combined summary plot matrix: peak p-value
gg_all_pval <- ggplot(
  data = allPlotData,
  mapping = aes(x = sampleId, y = peakPval)) +
  geom_quasirandom(color = "#567a0f") +
  geom_boxplot(width=0.1, fill = NA, outlier.colour = NA, color = alpha("black", 1)) +
  geom_hline(yintercept = 20, color = "red") +
  labs(title = "masc2 p-value distribution") +
  facet_wrap(facets = ~ sampleId, ncol = facetCols, scales = "free") +
  theme_plot_matrix

# pdf(file = paste(outPrefix, ".macs2_enrichment.pdf", sep = ""), width = 16, height = 16, onefile = TRUE)
png(filename = paste(outPrefix, ".macs2_pval.png", sep = ""), width = 3000, height = 3000, res = 350)
print(gg_all_pval)
dev.off()


## combined summary plot matrix: peak width
gg_all_width <- ggplot(
  data = allPlotData,
  mapping = aes(x = sampleId, y = peakWidth)) +
  geom_quasirandom(color = "#8472c0") +
  geom_boxplot(width=0.1, fill = NA, outlier.colour = NA, color = alpha("black", 1)) +
  geom_hline(yintercept = 300, color = "red") +
  labs(title = "macs2 peak width distribution") +
  facet_wrap(facets = ~ sampleId, ncol = facetCols, scales = "free") +
  theme_plot_matrix

# pdf(file = paste(outPrefix, ".peak_width.png", sep = ""), width = 16, height = 16, onefile = TRUE)
png(filename = paste(outPrefix, ".peak_width.png", sep = ""), width = 3000, height = 3000, res = 350)
print(gg_all_width)
dev.off()


##################################################################################




##################################################################################
# ## profile plot
# 
# outDir <- dirname(tfInfo$matFile[rowId])
# outPrefix <- paste(outDir, "/", tfInfo$sampleId[rowId], ".raw_peaks_summary", sep = "")
# 
# ## control samples to plot alongside TF ChIP sample
# tfControls <- c("AN10300_sCopy_OE_16h_input_FLAG_ChIPMix55_1",
#                 "AN10300_sCopy_OE_16h_input_FLAG_ChIPMix55_2",
#                 "AN0148_sCopy_OE_16h_HA_ChIPMix46_1",
#                 "AN2025_sCopy_OE_16h_HA_ChIPMix36_1")
# 
# 
# ## create profile matrix of 2kb region around peak summit for control samples
# peaksGr <- rtracklayer::import(con = tfInfo$peakFile[rowId], format = peakType)
# 
# if(length(peaksGr) > 0){
#   if(peakType == "broadPeak"){
#     mcols(peaksGr)$peak <- round(width(peaksGr) / 2)
#   }
#   
#   peakSummitGr <- GenomicRanges::narrow(x = peaksGr,
#                                         start = pmax(peaksGr$peak, 1),
#                                         width = 1)
#   
#   ctrlSampleInfo <- get_sample_information(
#     exptInfoFile = file_exptInfo,
#     samples = tfControls,
#     dataPath = TF_dataPath,
#     profileMatrixSuffix = matrixType)
#   
#   for (ctrlIdx in 1:nrow(ctrlSampleInfo)) {
#     profileMat <- bigwig_profile_matrix(bwFile = ctrlSampleInfo$bwFile[ctrlIdx],
#                                         regions = peakSummitGr,
#                                         signalName = ctrlSampleInfo$profileName[ctrlIdx],
#                                         extend = c(up, down),
#                                         targetName = "summit",
#                                         storeLocal = T,
#                                         localPath = ctrlSampleInfo$matFile[ctrlIdx])
#   }
#   
#   
#   
#   tempSampleInfo <- dplyr::bind_rows(tfInfo[rowId, ], ctrlSampleInfo) %>% 
#     dplyr::distinct()
#   
#   tfProfileMat <- import_profile_from_file(
#     file = tfInfo$matFile[rowId],
#     signalName = tfInfo$profileName[rowId],
#     selectGenes = peakSummitGr$name,
#     up = matrixDim[1], target = matrixDim[2], down = matrixDim[3],
#     targetType = "point", targetName = "summit" 
#   )
#   
#   quantile(tfProfileMat, c(seq(0, 0.9, by = 0.1), 0.95, 0.99, 0.992, 0.995, 0.999, 0.9999, 1), na.rm = T)
#   # tfMeanColor <- colorRamp2(quantile(tfProfileMat, c(0.50, 0.995), na.rm = T), c("white", "red"))
#   tfColorList <- sapply(
#     X = tempSampleInfo$sampleId,
#     FUN = function(x){
#       return(colorRamp2(breaks = quantile(tfProfileMat, c(0.50, 0.995), na.rm = T),
#                         colors = c("white", "red")))
#     }
#   )
#   
#   
#   
#   profilePlots <- multi_profile_plots(exptInfo = tempSampleInfo,
#                                       genesToPlot = peakSummitGr$name,
#                                       profileColors = tfColorList,
#                                       targetType = "point", targetName = "summit",
#                                       clusters = NULL,
#                                       showAnnotation = FALSE,
#                                       matBins = matrixDim,
#                                       column_title_gp = gpar(fontsize = 12))
#   
#   
#   rowOrd_peaks <- order(peakSummitGr$signalValue, decreasing = TRUE)
#   
#   # pdf(file = paste(outPrefix, "_profiles.pdf", sep = ""), width = 18, height = 13)
#   png(file = paste(outPrefix, ".profiles.png", sep = ""), width = 3500, height = 2500, res = 250)
#   
#   ht <- draw(
#     profilePlots$heatmapList,
#     main_heatmap = tempSampleInfo$profileName[1],
#     row_order = rowOrd_peaks,
#     column_title = paste(tfInfo$sampleId[rowId], "binding signal around 2kb region of macs2 peak summit"),
#     column_title_gp = gpar(fontsize = 12, fontface = "bold"),
#     row_sub_title_side = "left",
#     heatmap_legend_side = "bottom",
#     gap = unit(7, "mm"),
#     padding = unit(rep(0.5, times = 4), "cm")
#   )
#   
#   dev.off()
#   
# }
# 
# 







