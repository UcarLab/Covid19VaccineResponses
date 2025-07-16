library(chromVAR)
library(motifmatchr)
library(Matrix)
library(SummarizedExperiment)
library(BiocParallel)
library(BSgenome.Hsapiens.UCSC.hg38)
library(TFBSTools)
library(JASPAR2024)

register(MulticoreParam(24))
file_name = 'peak_matrix.txt' # Load the data Peak-Count matrix
out_name = substr(file_name, start = 13, stop = nchar(file_name) - 11)

countmatrix = read.table(file_name)

motif_file = readRDS('motif_files/hocomoco_pcm.rds') # hocomoco motif file

peaksdf = read.table(text=rownames(countmatrix), sep="_")
colnames(peaksdf) = c("chr", "start", "end")
peakgranges = makeGRangesFromDataFrame(peaksdf)
countmatrix <- as.matrix(countmatrix)
fragment_counts <- SummarizedExperiment(assays = 
                                          list(counts = countmatrix),
                                        rowRanges = peakgranges, colData = DataFrame(celltype = colnames(countmatrix)))
fragment_counts <- addGCBias(fragment_counts, genome = BSgenome.Hsapiens.UCSC.hg38)
fragment_counts_filtered = filterPeaks(fragment_counts, non_overlapping = TRUE)

motif_ix <- matchMotifs(motif_file, fragment_counts_filtered, genome = BSgenome.Hsapiens.UCSC.hg38)
dev <- computeDeviations(object = fragment_counts_filtered, annotations = motif_ix)
barcode_deviations = t(dev@assays@data$deviations)
z_score = t(dev@assays@data$z)

write.csv(barcode_deviations, paste(out_name, '_adjusted_v2.csv', sep=""), sep = '\t')
write.csv(z_score, paste('chromvar_result/', out_name, '_zscore_v2.csv', sep=""), sep = '\t')

