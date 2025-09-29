---
title: "ExpDH_GSEA"
author: "I. Hanka"
date: "2025-09-29"
output: 
  html_document:
    keep_md: true
    toc: true
editor_options: 
  chunk_output_type: inline
---





## load packages



## load object

``` r
fileNam <- "/Users/immbio/Desktop/HumanHeartCarTrans2/data/Human_heart_ExpDH.rds"
seuratExpDH <- readRDS(fileNam)
```

# Run marker Genes

``` r
##cluster marker
Idents(seuratExpDH) <- seuratExpDH$RNA_snn_res.0.4
markerGenesExpDH <- FindAllMarkers(seuratExpDH, only.pos=T) %>% 
  dplyr::filter(p_val_adj < 0.01)

#save table
write.table(markerGenesExpDH, 
            file= "/Users/immbio/Desktop/HumanHeartCarTrans2/analysis/ExpDHmarkerGenesRNA_snn_res.0.4",
            sep="\t",
            quote=F,
            row.names=F,
            col.names=T)
```



``` r
#dittoSeq::dittoDotPlot(seuratExpDH,vars = unique(markerGenesExpDH), group.by = "clusterName")
```


#Differential analysis using muscat
https://github.com/HelenaLC/muscat

``` r
#Preparation: Set as ScExperiment
sceExpDH <- as.SingleCellExperiment(seuratExpDH)
#Compute sum factors and normalize
sceExpDH <- computeLibraryFactors(sceExpDH)
sceExpDH <- logNormCounts(sceExpDH)

# Step 1: Add metadata to colData
colData(sceExpDH)$sample_id <- sceExpDH$Sample
colData(sceExpDH)$cluster_id <- sceExpDH$clusterName
colData(sceExpDH)$group_id <- sceExpDH$diseaseCond

# Step 2: Standardize structure
sceExpDH <- prepSCE(sceExpDH,
    kid = "cluster_id",
    gid = "group_id",
    sid = "sample_id")

# Step 3: Pseudobulk aggregation
pb <- aggregateData(sceExpDH, 
    assay = "counts", fun = "sum",
    by = c("cluster_id", "sample_id"))


# run pseudobulk (aggregation-based) DS analysis (analyze differential gene expression in specific cell subpopulations or states, comparing conditions)
ds_pb <- pbDS(pb, method = "edgeR")
```

```
##   |                                                                                                  |                                                                                          |   0%  |                                                                                                  |=====                                                                                     |   5%  |                                                                                                  |=========                                                                                 |  11%  |                                                                                                  |==============                                                                            |  16%  |                                                                                                  |===================                                                                       |  21%  |                                                                                                  |========================                                                                  |  26%  |                                                                                                  |============================                                                              |  32%  |                                                                                                  |=================================                                                         |  37%  |                                                                                                  |======================================                                                    |  42%  |                                                                                                  |===========================================                                               |  47%  |                                                                                                  |===============================================                                           |  53%  |                                                                                                  |====================================================                                      |  58%  |                                                                                                  |=========================================================                                 |  63%  |                                                                                                  |==============================================================                            |  68%  |                                                                                                  |==================================================================                        |  74%  |                                                                                                  |=======================================================================                   |  79%  |                                                                                                  |============================================================================              |  84%  |                                                                                                  |=================================================================================         |  89%  |                                                                                                  |=====================================================================================     |  95%  |                                                                                                  |==========================================================================================| 100%
```

``` r
#print(ds_pb)

#run pseudobulk multidimensional scaling (explore and visualize global patterns of expression similarity/differences across samples or groups)
pb_mds <- pbMDS(pb)
# use very distinctive shaping of groups & change cluster colors
#pb_mds <- pb_mds + 
#  scale_shape_manual(values = c(17, 4)) +
#  scale_color_manual(values = RColorBrewer::brewer.pal(8, "Set2"))
# change point size & alpha
  #pb_mds$layers[[1]]$aes_params$size <- 5
  #pb_mds$layers[[1]]$aes_params$alpha <- 0.6
#print(pb_mds)
```


## session info

``` r
date()
```

```
## [1] "Mon Sep 29 21:42:22 2025"
```

``` r
sessionInfo()
```

```
## R version 4.5.1 (2025-06-13)
## Platform: aarch64-apple-darwin20
## Running under: macOS Sequoia 15.1
## 
## Matrix products: default
## BLAS:   /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/lib/libRblas.0.dylib 
## LAPACK: /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.1
## 
## locale:
## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
## 
## time zone: Europe/Zurich
## tzcode source: internal
## 
## attached base packages:
## [1] grid      stats4    stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
##  [1] NCmisc_1.2.0                VennDiagram_1.7.3           futile.logger_1.4.3        
##  [4] ggupset_0.4.1               gridExtra_2.3               DOSE_4.2.0                 
##  [7] enrichplot_1.28.4           msigdbr_25.1.1              org.Hs.eg.db_3.21.0        
## [10] AnnotationDbi_1.70.0        clusterProfiler_4.16.0      multtest_2.64.0            
## [13] metap_1.12                  scater_1.35.0               scuttle_1.18.0             
## [16] destiny_3.22.0              circlize_0.4.16             edgeR_4.6.3                
## [19] limma_3.64.3                muscat_1.22.0               viridis_0.6.5              
## [22] viridisLite_0.4.2           lubridate_1.9.4             forcats_1.0.1              
## [25] stringr_1.5.2               purrr_1.1.0                 readr_2.1.5                
## [28] tidyr_1.3.1                 tibble_3.3.0                tidyverse_2.0.0            
## [31] dplyr_1.1.4                 SingleCellExperiment_1.30.1 SummarizedExperiment_1.38.1
## [34] Biobase_2.68.0              GenomicRanges_1.60.0        GenomeInfoDb_1.44.3        
## [37] IRanges_2.42.0              S4Vectors_0.46.0            BiocGenerics_0.54.0        
## [40] generics_0.1.4              MatrixGenerics_1.20.0       matrixStats_1.5.0          
## [43] pheatmap_1.0.13             ggpubr_0.6.1                ggplot2_4.0.0              
## [46] Seurat_5.3.0                SeuratObject_5.2.0          sp_2.2-0                   
## 
## loaded via a namespace (and not attached):
##   [1] R.methodsS3_1.8.2        progress_1.2.3           nnet_7.3-20             
##   [4] goftest_1.2-3            Biostrings_2.76.0        TH.data_1.1-4           
##   [7] vctrs_0.6.5              ggtangle_0.0.7           spatstat.random_3.4-2   
##  [10] digest_0.6.37            png_0.1-8                corpcor_1.6.10          
##  [13] shape_1.4.6.1            proxy_0.4-27             ggrepel_0.9.6           
##  [16] deldir_2.0-4             parallelly_1.45.1        MASS_7.3-65             
##  [19] reshape2_1.4.4           qvalue_2.40.0            httpuv_1.6.16           
##  [22] foreach_1.5.2            withr_3.0.2              ggfun_0.2.0             
##  [25] xfun_0.53                survival_3.8-3           memoise_2.0.1           
##  [28] hexbin_1.28.5            ggbeeswarm_0.7.2         gson_0.1.0              
##  [31] tidytree_0.4.6           zoo_1.8-14               GlobalOptions_0.1.2     
##  [34] gtools_3.9.5             pbapply_1.7-4            R.oo_1.27.1             
##  [37] DEoptimR_1.1-4           Formula_1.2-5            prettyunits_1.2.0       
##  [40] KEGGREST_1.48.1          promises_1.3.3           scatterplot3d_0.3-44    
##  [43] httr_1.4.7               rstatix_0.7.2            globals_0.18.0          
##  [46] fitdistrplus_1.2-4       rstudioapi_0.17.1        UCSC.utils_1.4.0        
##  [49] miniUI_0.1.2             babelgene_22.9           curl_7.0.0              
##  [52] ScaledMatrix_1.16.0      polyclip_1.10-7          TFisher_0.2.0           
##  [55] GenomeInfoDbData_1.2.14  SparseArray_1.8.1        RcppEigen_0.3.4.0.2     
##  [58] xtable_1.8-4             doParallel_1.0.17        evaluate_1.0.5          
##  [61] S4Arrays_1.8.1           hms_1.1.3                irlba_2.3.5.1           
##  [64] colorspace_2.1-2         ROCR_1.0-11              reticulate_1.43.0       
##  [67] spatstat.data_3.1-8      magrittr_2.0.4           lmtest_0.9-40           
##  [70] ggtree_3.16.3            later_1.4.4              lattice_0.22-7          
##  [73] spatstat.geom_3.6-0      future.apply_1.20.0      robustbase_0.99-6       
##  [76] scattermore_1.2          cowplot_1.2.0            RcppAnnoy_0.0.22        
##  [79] xts_0.14.1               class_7.3-23             pillar_1.11.1           
##  [82] nlme_3.1-168             iterators_1.0.14         caTools_1.18.3          
##  [85] compiler_4.5.1           beachmat_2.24.0          RSpectra_0.16-2         
##  [88] stringi_1.8.7            tensor_1.5.1             minqa_1.2.8             
##  [91] plyr_1.8.9               crayon_1.5.3             abind_1.4-8             
##  [94] blme_1.0-6               gridGraphics_0.5-1       sn_2.1.1                
##  [97] locfit_1.5-9.12          bit_4.6.0                mathjaxr_1.8-0          
## [100] sandwich_3.1-1           pcaMethods_2.0.0         fastmatch_1.1-6         
## [103] multcomp_1.4-28          codetools_0.2-20         BiocSingular_1.24.0     
## [106] TTR_0.24.4               bslib_0.9.0              e1071_1.7-16            
## [109] GetoptLong_1.0.5         ggplot.multistats_1.0.1  plotly_4.11.0           
## [112] remaCor_0.0.20           mime_0.13                splines_4.5.1           
## [115] Rcpp_1.1.0               fastDummies_1.7.5        blob_1.2.4              
## [118] knitr_1.50               clue_0.3-66              lme4_1.1-37             
## [121] fs_1.6.6                 listenv_0.9.1            Rdpack_2.6.4            
## [124] ggplotify_0.1.3          ggsignif_0.6.4           Matrix_1.7-4            
## [127] statmod_1.5.0            tzdb_0.5.0               fANCOVA_0.6-1           
## [130] pkgconfig_2.0.3          tools_4.5.1              cachem_1.1.0            
## [133] RSQLite_2.4.3            RhpcBLASctl_0.23-42      rbibutils_2.3           
## [136] DBI_1.2.3                smoother_1.3             numDeriv_2016.8-1.1     
## [139] fastmap_1.2.0            rmarkdown_2.29           scales_1.4.0            
## [142] ica_1.0-3                broom_1.0.10             sass_0.4.10             
## [145] patchwork_1.3.2          dotCall64_1.2            carData_3.0-5           
## [148] RANN_2.6.2               farver_2.1.2             reformulas_0.4.1        
## [151] aod_1.3.3                mgcv_1.9-3               yaml_2.3.10             
## [154] ggthemes_5.1.0           cli_3.6.5                lifecycle_1.0.4         
## [157] uwot_0.2.3               glmmTMB_1.1.12           mvtnorm_1.3-3           
## [160] lambda.r_1.2.4           backports_1.5.0          BiocParallel_1.42.2     
## [163] timechange_0.3.0         gtable_0.3.6             rjson_0.2.23            
## [166] ggridges_0.5.7           progressr_0.16.0         ape_5.8-1               
## [169] parallel_4.5.1           jsonlite_2.0.0           RcppHNSW_0.6.0          
## [172] bitops_1.0-9             assertthat_0.2.1         bit64_4.6.0-1           
## [175] qqconf_1.3.2             Rtsne_0.17               yulab.utils_0.2.1       
## [178] spatstat.utils_3.2-0     BiocNeighbors_2.2.0      ranger_0.17.0           
## [181] mutoss_0.1-13            futile.options_1.0.1     jquerylib_0.1.4         
## [184] GOSemSim_2.34.0          R.utils_2.13.0           spatstat.univar_3.1-4   
## [187] pbkrtest_0.5.5           lazyeval_0.2.2           shiny_1.11.1            
## [190] htmltools_0.5.8.1        GO.db_3.21.0             sctransform_0.4.2       
## [193] formatR_1.14             rappdirs_0.3.3           glue_1.8.0              
## [196] spam_2.11-1              XVector_0.48.0           VIM_6.2.6               
## [199] treeio_1.32.0            mnormt_2.1.1             EnvStats_3.1.0          
## [202] boot_1.3-32              igraph_2.1.4             variancePartition_1.38.1
## [205] TMB_1.9.17               R6_2.6.1                 DESeq2_1.48.2           
## [208] gplots_3.2.0             vcd_1.4-13               cluster_2.1.8.1         
## [211] aplot_0.2.9              nloptr_2.2.1             plotrix_3.8-4           
## [214] DelayedArray_0.34.1      tidyselect_1.2.1         vipor_0.4.7             
## [217] car_3.1-3                future_1.67.0            rsvd_1.0.5              
## [220] KernSmooth_2.23-26       S7_0.2.0                 laeken_0.5.3            
## [223] data.table_1.17.8        fgsea_1.34.2             htmlwidgets_1.6.4       
## [226] ComplexHeatmap_2.24.1    RColorBrewer_1.1-3       rlang_1.1.6             
## [229] spatstat.sparse_3.1-0    spatstat.explore_3.5-3   lmerTest_3.1-3          
## [232] beeswarm_0.4.0
```

