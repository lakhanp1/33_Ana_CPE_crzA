suppressPackageStartupMessages(library(chipmine))
suppressPackageStartupMessages(library(markPeaks))
suppressPackageStartupMessages(library(org.AFumigatus.A1163.eg.db))
suppressPackageStartupMessages(library(TxDb.Afumigatus.A1163.AspGD.GFF))
suppressPackageStartupMessages(library(BSgenome.Afumigatus.A1163.AspGD))
suppressPackageStartupMessages(library(here))


## 1) annotate peaks
## 2) create gene level peak annotation data

rm(list = ls())

source(file = "E:/Chris_UM/GitHub/omics_util/04_GO_enrichment/s01_enrichment_functions.R")

##################################################################################

file_exptInfo <- here::here("data", "reference_data", "ChIPseq_sample_info.tab")

orgDb <- org.AFumigatus.A1163.eg.db
# txDb <- get_TxDb_sqlite(org = "A_fumigatus_A1163", yaml = "E:/Chris_UM/Database/TxDb.config.yaml")
txDb <- TxDb.Afumigatus.A1163.AspGD.GFF
bsGenome <- BSgenome.Afumigatus.A1163.AspGD

TF_dataPath <- here::here("data", "ChIPseq_CEA17")

file_tf_macs2 <- paste(TF_dataPath, "/", "sample_tf_macs2.list", sep = "")

matrixType <- "2kb_summit"
up <- 2000
down <- 2000
body <- 0
bin <- 10
matrixDim = c(c(up, body, down)/bin, bin)

##################################################################################
geneSet <- AnnotationDbi::select(
  x = orgDb, keys = keys(x = orgDb, keytype = "GID"),
  columns = "DESCRIPTION", keytype = "GID"
) %>% 
  dplyr::rename(geneId = GID)

txInfo <- suppressMessages(
  AnnotationDbi::select(
    x = txDb, keys = AnnotationDbi::keys(x = txDb, keytype = "TXID"),
    columns = c("GENEID", "TXNAME", "TXTYPE"), keytype = "TXID")) %>%
  dplyr::mutate(TXID = as.character(TXID)) %>%
  dplyr::rename(geneId = GENEID, txName = TXNAME, txType = TXTYPE)

txInfo <- dplyr::filter(txInfo, !txType %in% c("tRNA", "rRNA", "snRNA", "snoRNA")) %>% 
  dplyr::filter(!grepl(pattern = "uORF", x = geneId))

tfSampleList <- suppressMessages(readr::read_tsv(file = file_tf_macs2))

tfInfo <- get_sample_information(
  exptInfoFile = file_exptInfo,
  samples = tfSampleList$sampleId,
  dataPath = TF_dataPath,
  profileMatrixSuffix = matrixType)

##################################################################################

i <- 1


for(i in 1:nrow(tfInfo)){
  
  cat(i, ":", tfInfo$sampleId[i], "\n")
  
  ## annotate peaks and prepare gene level annotation file
  peakType <- tfInfo$peakType[i]
  
  peakAn <- markPeaks::annotate_peaks(
    peakFile = tfInfo$peakFile[i],
    txdb = txDb,
    txIds = txInfo$TXID,
    fileFormat = peakType,
    promoterLength = 500,
    upstreamLimit = 1000,
    bidirectionalDistance = 1000,
    includeFractionCut = 0.7,
    bindingInGene = FALSE,
    insideSkewToEndCut = 0.7,
    removePseudo = TRUE,
    output = tfInfo$peakAnno[i])
  
  if( !is.null(peakAn) ){
    tfDf <- gene_level_peak_annotation(
      sampleId = tfInfo$sampleId[i],
      peakAnnotation = tfInfo$peakAnno[i],
      genesDf = geneSet,
      peakFile = tfInfo$peakFile[i],
      bwFile = tfInfo$bwFile[i],
      outFile = tfInfo$peakTargetFile[i])
  }
  
  
  # ## create profile matrix of 2kb region around peak summit
  # peaksGr <- rtracklayer::import(con = tfInfo$peakFile[i], format = peakType)
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
  #   profileMat <- bigwig_profile_matrix(bwFile = tfInfo$bwFile[i],
  #                                       regions = peakSummitGr,
  #                                       signalName = tfInfo$profileName[i],
  #                                       extend = c(up, down),
  #                                       targetName = "summit",
  #                                       storeLocal = T,
  #                                       localPath = tfInfo$matFile[i])
  # }
  
  ## extract 500bp sequence around peak summit
  peakData <- import_peaks_as_df(file = tfInfo$peakFile[i])
  summitSeq <- get_peak_summit_seq(
    file = tfInfo$peakFile[i], peakFormat = peakType,
    genome = bsGenome, length = 500
  )
  
  summitSeq <- dplyr::left_join(x = summitSeq, y = peakData, by = "peakId")
  
  readr::write_tsv(
    x = summitSeq,
    file = paste(dirname(tfInfo$peakAnno[i]),"/", tfInfo$sampleId[i], ".summitSeq500.tab", sep = "")
  )
  
}



GenomicFeatures::cds(x = txDb, filter = list(tx_id = txInfo$TXID))

cdsGr <- unlist(range(GenomicFeatures::cdsBy(x = txDb, by = "gene")))

mcols(cdsGr)$name <- names(cdsGr)

atgGr <- GenomicRanges::resize(x = cdsGr, width = 1, fix = "start")
promoterGr <- trim(promoters(x = atgGr, upstream = 400, downstream = 100))

promoterGr$seq <- BSgenome::getSeq(x = bsGenome, names = promoterGr, as.character = TRUE)

promoterBs <- BSgenome::getSeq(x = bsGenome, names = promoterGr)

Biostrings::writeXStringSet(x = promoterBs, filepath = "A_fumigatus_A1163.CDS.upstream500.seq.fasta")





