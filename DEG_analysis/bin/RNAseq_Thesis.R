#######################################################
########Differential Gene Expression with DESeq2#######
#######################################################


############################
########Envirornment########
############################
#Packages needed

packages <- c("DESeq2", "pheatmap", "genefilter", "dplyr", "ggplot2", "reshape2", "ggrepel", "pathview", "magrittr", "biomaRt",
              "gridExtra", "grid", "lattice", "ggplotify", "ape")
lapply(packages, require, character.only = TRUE)
remove(packages)

#Functions
#Antilog2 function
applyAntiLogFC <- function(logFC){
  
  antilog.2 <- logFC
  
  if (is.na(logFC)) {
    antilog.2
  } else {
    if (logFC >= 0){
      antilog.2 <- 2^logFC
    } else {
      antilog.2 <- -(1/2^logFC)
    }
  }
  
  antilog.2
}

#Volcano plot
volcanoPlot <- function(x) {
  df <- as.data.frame(x)
  df$padj[is.na(df$padj)] = 1 #change na values to 1's in padj
  mutateddf <- mutate(df, Significant = ifelse(df$padj < 0.01, "padj < 0.01", "Not Sig")) #Will have different colors depending on significance
  #Gene lables to plot
  bpGenes <- c("CYFIP1", "CYFIP2", "NIPA1", "NIPA2", "TUBGCP5") 
  annotateGenes <- rbind(head(mutateddf, 10),
                         mutateddf[mutateddf$geneName %in% bpGenes,])
  annotateGenes <- annotateGenes[!duplicated(annotateGenes$entrezID),]
  #Volcano plot
  mutateddf %>% ggplot() + geom_jitter(aes(log2FoldChange, -log10(padj), colour = Significant)) +
    theme_bw() +
    theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
          axis.text.x = element_text(angle = 90, hjust = 1, size = 20, vjust = 0.5), 
          axis.title.x = element_blank(), axis.title.y = element_text(size = 20), 
          axis.text.y = element_text(size = 20), axis.text.y.right = element_text(size = 20), 
          legend.text = element_text(size = 18), legend.title = element_text(size = 18)) + 
    guides(colour = guide_legend(override.aes = list(size = 4))) + 
    geom_text_repel(data = annotateGenes, aes(log2FoldChange, -log10(padj), label=geneName), size = 5, arrow = arrow(length = unit(0.02, "npc"), type = "open", ends = "first")) 
  
  #ggsave(paste("../output/", ".jpeg", sep = ""), device = "jpeg")
  
}


ensembl = useMart("ensembl",dataset = "hsapiens_gene_ensembl") #biomaRt access

countfilter <- 10 #minimum cut-off for gene expression levels

resources <- read.table("../resources/targets.txt", header = T, check.names = F, sep = "\t", quote = "") #data frame containing all the relevant information for each comparison 
#############################################
########Differential Gene Expression#########
#############################################

stats <- data.frame() #data frame to collect general statistics of each analysis
for (file in list.files(path = "../input")) {
  print(file)
  rawData <- read.table(paste("../input/", file, sep = ""), header = T, check.names = F, sep = "\t", quote = "", stringsAsFactors = F) #read in the count data
  
  targets <- resources[resources$sampleName %in% colnames(rawData[-1]),]
  targets <- targets[match(colnames(rawData[-1]), targets[,1]),]
  
  conditions <- paste(targets$phenotype, targets$time, sep = ".") #Define the elements for the regression model
  
  
  rownames(rawData) <- rawData$tracking_id
  rawData[,1] <- NULL
  rawData = rawData[apply(rawData, 1, function(row) {max(row) > countfilter}),] #Remove genes with less than 10 counts in at least one of the samples
  rawData <-  rawData[complete.cases(rawData), ]
  rawData <- cbind(rownames(rawData), rawData)
  colnames(rawData)[1] <- "tracking_id"
  
  exptDesign = data.frame(row.names = colnames(rawData[-1]),
                          condition = conditions) #Define the experimental design for the comparison
  
  exptObject <- DESeqDataSetFromMatrix(countData = rawData,
                                       tidy = TRUE,
                                       colData = exptDesign,
                                       design = ~ condition)
  
  analysisObject = DESeq(exptObject)
  
  rawCounts <- as.data.frame(counts(analysisObject, normalized = FALSE))
  colnames(rawCounts) <- gsub("^", "raw.counts.", colnames(rawCounts)) #add string indicating that the values are non-normalised (raw counts)
  rawCounts <- cbind(rownames(rawCounts), rawCounts)
  colnames(rawCounts)[1] <- "tracking_id"
  normalisedCounts <- as.data.frame(counts(analysisObject, normalized = TRUE))
  colnames(normalisedCounts) <- gsub("^", "norm.counts.", colnames(normalisedCounts)) #add string indicating that the values are normalised (norm counts)
  normalisedCounts <- cbind(rownames(normalisedCounts), normalisedCounts)
  colnames(normalisedCounts)[1] <- "tracking_id"
  
  mut <-  unique(conditions[grep("mut", conditions)])
  ctrl <- unique(conditions[grep("control", conditions)])
  
  result = as.data.frame(results(analysisObject, contrast = c("condition", mut, ctrl), independentFiltering = TRUE, pAdjustMethod = "BH")) #get the results of the differential gene expression by comparing mutant vs control
  result <- cbind(rownames(result), result)
  colnames(result)[1] <- "tracking_id"
  
  geneCounts <-  left_join(rawCounts, normalisedCounts, by = c("tracking_id" = "tracking_id"))
  fullData <- left_join(geneCounts, result, by = c("tracking_id" = "tracking_id"))
  fullData$FC <- unlist(lapply(fullData$log2FoldChange, applyAntiLogFC))
  
  
  #biomaRt query
  IDs <- getBM(filters = "ensembl_gene_id",
               attributes = c("ensembl_gene_id", "gene_biotype", "entrezgene", "external_gene_name", "description"),
               values = fullData$tracking_id,
               mart = ensembl)
  colnames(IDs) <- c("ensemblID", "gene_biotype", "entrezID", "geneName", "description")
  
  #remove duplicated values from biomaRt query and merge with full data
  deduped.data  <- IDs[match(fullData$tracking_id,IDs[,1]),]
  colnames(deduped.data)
  
  finalResults <- left_join(deduped.data, fullData, by = c("ensemblID" = "tracking_id"))
  finalResults = finalResults[order(finalResults$pvalue),]
  
  cell_1 <- as.character(unique(targets$cellType[grep("mut", targets$phenotype)]))
  stage_1 <- as.character(unique(targets$time[grep("mut", targets$phenotype)]))
  cell_2 <- as.character(unique(targets$cellType[grep("control", targets$phenotype)]))
  stage_2 <- as.character(unique(targets$time[grep("control", targets$phenotype)]))

  fileName <- paste(cell_1, ".", stage_1, "_vs_", cell_2, ".", stage_2, sep = "")
  outDir <- paste("../output/", fileName, sep = "")
  
  if (!dir.exists(outDir)) {
    dir.create(outDir, recursive = T)
  }
  
  write.table(finalResults, paste(outDir, "/", fileName, ".txt", sep = ""), col.names = T, row.names = F, sep = "\t", quote = F)
  
  data <- data.frame(cellLine = cell_1,
                     stage = stage_1,
                     comparison = fileName,
                     `DEGs` = nrow(finalResults),
                     `DEGs padj < 0.05` = nrow(subset(finalResults, padj < 0.05)),
                     `DEGs padj < 0.01` = nrow(subset(finalResults, padj < 0.01)),
                     `Protein coding padj < 0.05` = nrow(subset(finalResults, padj < 0.05 & gene_biotype == "protein_coding")),
                     `Protein coding padj < 0.01` = nrow(subset(finalResults, padj < 0.01 & gene_biotype == "protein_coding")))
  
  stats <- rbind(stats, data)
  
  volcano <- volcanoPlot(finalResults)
  
  topSig <- as.character(subset(finalResults, padj < 0.01)[,1]) #character vector with the top significant genes (1% FDR)
  
  resLFC <- lfcShrink(analysisObject, contrast = c("condition", ctrl, mut))
 
  
  
  #Logtransform
  vsdObject <- varianceStabilizingTransformation(analysisObject)
 
  #1Generate the matrix of expression
  mat  <- assay(vsdObject)[ topSig, ]
  #2Transform to log scale and remove Inf/NA's and normalise agains the mean row expression for each gene
  #logmat <- as.data.frame(log2(mat))
  #logmat <- data.frame(row.names = rownames(logmat), do.call(data.frame,lapply(logmat, function(f) replace(f, is.infinite(f),0))))
  normlogMat  <- mat - rowMeans(mat)
  #3 generate annotation data frame
  annDf <-  droplevels.data.frame(targets[targets$sampleName %in% gsub("raw.counts.", "", rownames(colData(vsdObject))),])
  rownames(annDf) <- annDf$sampleName
  
  annCol <- data.frame(row.names = rownames(annDf),
                       `Cell Type` = annDf$cellType)
  
  title <- paste(unique(levels(annDf$cellType))[1], "vs", unique(levels(annDf$cellType))[2], unique(levels(annDf$time)), sep = " ")
  
  
  ####PCA
  annDf$names <- rownames(annDf)
  pcaData <- plotPCA(vsdObject[ topSig, ], returnData = T)
  scores <- left_join(pcaData, annDf, by = c("name" = "sampleName"))
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  
  #png(filename = paste("~/Desktop/", names(x),".png", sep = ""))
  pcaPlot <- scores %>% ggplot() + geom_point(aes(PC1, PC2, color = cellType, shape = time), size = 4) + 
    theme(axis.text.x = element_text(hjust = 1, size = 15, vjust = 0.5),
          axis.text.y = element_text(hjust = 1, size = 15, vjust = 0.5),
          axis.title.x = element_text(size = 20),
          axis.title.y = element_text(size = 20),
          legend.text = element_text(size = 18), legend.title = element_text(size = 18)) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) + 
    ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
    geom_text_repel(data = scores, aes(PC1, PC2, label = name), size = 5, arrow = arrow(length = unit(0.02, "npc"), type = "open", ends = "first")) 
  
  #Boxplots
  colnames(rawCounts) <- gsub("raw.counts.", "", colnames(rawCounts))
  log2BoxRaw <- log2(rawCounts[-1])[ topSig, ]
  log2BoxRaw = as.data.frame(t(log2BoxRaw))
  log2BoxRaw$Sample <- rownames(log2BoxRaw)
  log2BoxRaw.m <- melt(log2BoxRaw)
  log2BoxRaw.m$Sample <- gsub("raw.counts.", "", log2BoxRaw.m$Sample)
  log2BoxRaw.m <- left_join(log2BoxRaw.m, annDf, by = c("Sample" = "sampleName"))
  rawBoxPlot <- log2BoxRaw.m %>% ggplot() + 
    geom_boxplot(aes(Sample, value, fill = phenotype), outlier.color = "red", outlier.shape = 16, outlier.size = 2, notch = T) + 
    ggtitle(label = "Raw counts") +
    theme(axis.text.x = element_text(hjust = 0.5, size = 15, vjust = 0.5),
          axis.text.y = element_text(hjust = 0.5, size = 15, vjust = 0.5),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.text = element_text(size = 18), legend.title = element_text(size = 18))
  
  colnames(normalisedCounts) <- gsub("norm.counts.", "", colnames(normalisedCounts))
  log2BoxNorm <- log2(normalisedCounts[-1])[ topSig, ]
  log2BoxNorm = as.data.frame(t(log2BoxNorm))
  log2BoxNorm$Sample <- rownames(log2BoxNorm)
  log2BoxNorm.m <- melt(log2BoxNorm)
  log2BoxNorm.m$Sample <- gsub("norm.counts.", "", log2BoxNorm.m$Sample)
  log2BoxNorm.m <- left_join(log2BoxNorm.m, annDf, by = c("Sample" = "sampleName"))
  normBoxPlot <- log2BoxNorm.m %>% ggplot() + 
    geom_boxplot(aes(Sample, value, fill = phenotype), outlier.color = "red", outlier.shape = 16, outlier.size = 2, notch = T) + 
    ggtitle(label = "Normalised counts") +
    theme(axis.text.x = element_text(hjust = 0.5, size = 15, vjust = 0.5),
          axis.text.y = element_text(hjust = 0.5, size = 15, vjust = 0.5),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.text = element_text(size = 18), legend.title = element_text(size = 18))
  
  d <- cor(rawCounts[-1][ topSig, ], method = "pearson")
  hcRaw <- hclust(dist(1 - d))
  d <- cor(normalisedCounts[-1][ topSig, ], method = "pearson")
  hcNorm <- hclust(dist(1 - d))
  

  #Report

  #pheatmap(as.matrix(normlogMat), scale = "row", annotation_col = annCol, show_rownames = FALSE, show_colnames = FALSE, labels_col = annDf$sampleName, main = title)
  par(mfrow = c(2,2))
  plotMA(resLFC,  ylim = c(-2,2)) #ma plot
  plotDispEsts(analysisObject) #dispersion estimates and fitting model

  
  pdf(paste(outDir, "/", fileName, ".pdf", sep = ""), width = 16, height = 14)
  grid.arrange(rawBoxPlot, normBoxPlot, volcano, pcaPlot,
               ncol = 2, nrow = 2)
  par(mfrow = c(2,2))
  plotMA(resLFC,  ylim = c(-2,2)) #ma plot
  plotDispEsts(analysisObject) #dispersion estimates and fitting model
  #Philogenetic trees
  plot.phylo(as.phylo(hcRaw), type = "p", show.node.label = TRUE, main = "Dendrogram: raw expression")
  plot.phylo(as.phylo(hcNorm), type = "p", show.node.label = TRUE, main = "Dendrogram: normalised expression")
  #Heatmap
  pheatmap(as.matrix(normlogMat), scale = "row", annotation_col = annCol, show_rownames = FALSE, show_colnames = FALSE, labels_col = annDf$sampleName, main = title)
  
  graphics.off()
  gc()
}

write.table(stats, file = paste("../output/", "dge_summary.txt", sep = ""), col.names = T, row.names = F, sep = "\t", quote = F)

#save.image(file = "dge.RData")
