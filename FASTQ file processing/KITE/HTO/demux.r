library(Seurat)
library(stringr)
args = commandArgs(trailingOnly=TRUE)

hto_path = gsub("/$", "", args[1])
adt_path = gsub("/$", "", args[2])
run_id = gsub("/$", "", args[3])

paths.hto = paste0('KITE/HTO/', hto_path)
paths.adt = paste0('KITE/ADT/', adt_path)
paths.atac_gex = paste0('CellRanger/', run_id, '-out/outs/filtered_feature_bc_matrix.h5')
paths.output = paste0(run_id, '-demuxed_hashes.csv')

print(paths.hto)
print(paths.adt)
print(paths.atac_gex)

pbmc_htos = Read10X(paths.hto, gene.column=1)
colnames(pbmc_htos) = gsub('$', '-1', colnames(pbmc_htos))

pbmc_adts = Read10X(paths.adt, gene.column=1)
colnames(pbmc_adts) = gsub('$', '-1', colnames(pbmc_adts))

pbmc = Read10X_h5(paths.atac_gex)

# pbmc_atac = pbmc[["Peaks"]]
pbmc_umis = pbmc[["Gene Expression"]]

# stopifnot(unique(str_sub(colnames(pbmc_atac), start=-2))=="-1")
stopifnot(unique(str_sub(colnames(pbmc_umis), start=-2))=="-1")

joint_bcs = intersect(colnames(pbmc_umis), colnames(pbmc_htos))
joint_bcs = intersect(joint_bcs, colnames(pbmc_adts))

pbmc_htos = as.matrix(pbmc_htos[, joint_bcs])



pbmc_hashtag <- CreateSeuratObject(counts = pbmc_htos, assay = 'HTO')
# Normalize HTO data, here we use centered log-ratio (CLR) transformation
pbmc_hashtag = NormalizeData(pbmc_hashtag, assay = "HTO", normalization.method = "CLR")


# # If you have a very large dataset we suggest using k_function = 'clara'. This is a k-medoid
# # clustering function for large applications You can also play with additional parameters (see
# # documentation for HTODemux()) to adjust the threshold for classification Here we are using
# # the default settings
pbmc_hashtag = HTODemux(pbmc_hashtag, assay = "HTO", positive.quantile = 0.99)

Idents(pbmc_hashtag) <- "HTO_maxID"

pdf(paste0(run_id, '-hto_ridge_plot.pdf'), height = 10, width = 30)
ridge_plot = RidgePlot(pbmc_hashtag, assay = "HTO", features = rownames(pbmc_hashtag[["HTO"]]), ncol = 4)
print(ridge_plot)
dev.off()

print(table(pbmc_hashtag$hash.ID))

write.csv(pbmc_hashtag[['hash.ID']], paths.output)
