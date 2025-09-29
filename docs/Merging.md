---
title: "HumanHeartProject - Merging"
author: "I. Hanka"
date: "2025-09-24"
output: 
  html_document:
    keep_md: true
    toc: true
editor_options: 
  chunk_output_type: inline
---




``` r
#MaxVSize manuell festlegen
options(future.globals.maxVSize=45000*1024^2)
#options("maxVSize"=)
mem.maxVSize()
```

```
## [1] 153600
```

## load packages


#start pre-processing

## load files and merge

``` r
### load and merge all 
basedir <- "/Users/immbio/Desktop/HumanHeartCarTrans2/data/seurat files/"
fileNamList <- list.files(path = basedir)

for(i in 1:length(fileNamList)){
  seuratS <- readRDS(paste0(basedir, fileNamList[i]))
  seuratS@assays$RNA@layers$scale.data <- NULL ##remove scale data slot to reduce size of objects
  if(exists("seuratM")){
    seuratM <- merge(x = seuratM, y = seuratS)
    cat("Merged", i, "of", length(fileNamList), "files - Total cells:", ncol(seuratM), "\n")
  }else{
    seuratM <- seuratS
    cat("Initialized with first Seurat object:", fileNamList[i], "\n")
  }
}

remove(seuratS)
table(seuratM$dataset)
table(seuratM$orig.ident)

#join layers
seuratM <- JoinLayers(seuratM)

#rerun seurat
seuratM <- NormalizeData (object = seuratM)
seuratM <- FindVariableFeatures(object = seuratM)
seuratM <- ScaleData(object = seuratM, verbose = TRUE)
seuratM <- RunPCA(object=seuratM, npcs = 30, verbose = FALSE)
seuratM <- RunTSNE(object=seuratM, reduction="pca", dims = 1:20)
seuratM <- RunUMAP(object=seuratM, reduction="pca", dims = 1:20)
seuratM <- FindNeighbors(object = seuratM, reduction = "pca", dims= 1:20)

res <- c(0.25, 0.6, 0.8, 0.4)
for (i in 1:length(res)) {
  seuratM <- FindClusters(object = seuratM, resolution = res[i], random.seed = 1234)
}
```

## save object

``` r
### save seurat object
saveRDS(seuratM, file="/Users/immbio/Desktop/HumanHeartCarTrans2/data/Human_heart_allmerged_seurat.rds")
```

