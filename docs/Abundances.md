---
title: "Abundances1_3"
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
fileNam <- "/Users/immbio/Desktop/HumanHeartCarTrans2/data/Human_heart_allmerged_seurat.rds"
seuratM <- readRDS(fileNam)
```

##set color vectors 

``` r
seuratM$clusterName <- "clusterName"
seuratM$clusterName[which(seuratM$RNA_snn_res.0.4 %in% "0" )] <- "Fb1"
seuratM$clusterName[which(seuratM$RNA_snn_res.0.4 %in% "1" )] <- "PerivFb1"
seuratM$clusterName[which(seuratM$RNA_snn_res.0.4 %in% "2" )] <- "Mph2"
seuratM$clusterName[which(seuratM$RNA_snn_res.0.4 %in% "3" )] <- "BEC1"
seuratM$clusterName[which(seuratM$RNA_snn_res.0.4 %in% "4" )] <- "Fb2"
seuratM$clusterName[which(seuratM$RNA_snn_res.0.4 %in% "5" )] <- "CM"
seuratM$clusterName[which(seuratM$RNA_snn_res.0.4 %in% "6" )] <- "Tcell1"
seuratM$clusterName[which(seuratM$RNA_snn_res.0.4 %in% "7" )] <- "BEC2"
seuratM$clusterName[which(seuratM$RNA_snn_res.0.4 %in% "8" )] <- "VSMC"
seuratM$clusterName[which(seuratM$RNA_snn_res.0.4 %in% "9" )] <- "Mph1"
seuratM$clusterName[which(seuratM$RNA_snn_res.0.4 %in% "10" )] <- "BEC3"
seuratM$clusterName[which(seuratM$RNA_snn_res.0.4 %in% "11" )] <- "NC"
seuratM$clusterName[which(seuratM$RNA_snn_res.0.4 %in% "12" )] <- "BaroRec"
seuratM$clusterName[which(seuratM$RNA_snn_res.0.4 %in% "13" )] <- "Bcell"
seuratM$clusterName[which(seuratM$RNA_snn_res.0.4 %in% "14" )] <- "Fb3"
seuratM$clusterName[which(seuratM$RNA_snn_res.0.4 %in% "15" )] <- "Tcell2"
seuratM$clusterName[which(seuratM$RNA_snn_res.0.4 %in% "16" )] <- "LEC"
seuratM$clusterName[which(seuratM$RNA_snn_res.0.4 %in% "17" )] <- "PerivFb2"
seuratM$clusterName[which(seuratM$RNA_snn_res.0.4 %in% "18" )] <- "Adipoc"

colclusterName <- c("#67001f", "#f4a582","#D53E4F", "#B45B5C","#003c30","#01665e","#66C2A5", "#BEAEF8","#BEAED4", "#c7eae5", "#B09C85", "#4e5a4c","#393A3F","pink","#4588CA","#3299CA","#FCC80B","#FEE60B","#628395")
names(colclusterName) <- c("CM","Fb1","Fb2","Fb3","PerivFb1","PerivFb2","VSMC","BEC1","BEC2","BEC3","LEC","NC","BaroRec","Adipoc","Mph1","Mph2","Tcell1","Tcell2","Bcell")

coldiseaseCond <- c("#dfc27d","#BE3144","#f4a582","#B45B5C","#8c510a","#202547","#355C7D","#779d8d", "#01665e", "#3288BD", "#BEAED4") 
names(coldiseaseCond) <- c("donorheart", "explant", "visit1", "visit2", "visit3", "visit4", "visit5", "visitX1", "visitX2", "visitX3", "visitX4")

#Patient
seuratM$patient <- "pat_nr"
seuratM$patient[grepl("HTx001|EMB001", seuratM$dataset)] <- "CarTransPat01"
seuratM$patient[grepl("HTx002|EMB002", seuratM$dataset)] <- "CarTransPat02"
seuratM$patient[grepl("HTx003|EMB003", seuratM$dataset)] <- "CarTransPat03"
seuratM$patient[grepl("HTx004|EMB004", seuratM$dataset)] <- "CarTransPat04"
seuratM$patient[grepl("HTx005|EMB005", seuratM$dataset)] <- "CarTransPat05"
seuratM$patient[grepl("HTx006|EMB006", seuratM$dataset)] <- "CarTransPat06"
seuratM$patient[grepl("HTx007|EMB007", seuratM$dataset)] <- "CarTransPat07"
seuratM$patient[grepl("HTx008|EMB008", seuratM$dataset)] <- "CarTransPat08"
seuratM$patient[grepl("HTx009|EMB009", seuratM$dataset)] <- "CarTransPat09"
seuratM$patient[grepl("HTx010|EMB010", seuratM$dataset)] <- "CarTransPat10"
seuratM$patient[grepl("HTx011|EMB011", seuratM$dataset)] <- "CarTransPat11"
seuratM$patient[grepl("HTx012|EMB012", seuratM$dataset)] <- "CarTransPat12"
seuratM$patient[grepl("HTx013|EMB013", seuratM$dataset)] <- "CarTransPat13"
seuratM$patient[grepl("HTx014|EMB014", seuratM$dataset)] <- "CarTransPat14"
seuratM$patient[grepl("HTx015|EMB015", seuratM$dataset)] <- "CarTransPat15"
seuratM$patient[grepl("HTx016|EMB016", seuratM$dataset)] <- "CarTransPat16"
seuratM$patient[grepl("HTx018|EMB018", seuratM$dataset)] <- "CarTransPat18"
seuratM$patient[grepl("HTx019|EMB019", seuratM$dataset)] <- "CarTransPat19"
seuratM$patient[grepl("HTx024|EMB024", seuratM$dataset)] <- "CarTransPat24"

seuratM$patient[which(seuratM$dataset == "o28576_1_08-8_20220525_Hu_nucseq_Graz_8_HH_GEM")] <- "DH01"
seuratM$patient[which(seuratM$dataset == "o28576_1_10-10_20220525_Hu_nucseq_Graz_10_HH_GEM")] <- "DH02"
seuratM$patient[which(seuratM$dataset == "o28576_1_11-11_20220525_Hu_nucseq_Graz_11_HH_GEM")] <- "DH03"
seuratM$patient[which(seuratM$dataset == "o28576_1_12-12_20220525_Hu_nucseq_Graz_12_HH_GEM")] <- "DH04"
seuratM$patient[which(seuratM$dataset =="o292731_1-1_20220818_Hu_nucseq_Graz_9_HH_GEM")] <- "DH05"
seuratM$patient[which(seuratM$dataset =="o292731_2-2_20220818_Hu_nucseq_Graz_13_HH_GEM")] <- "DH06"
seuratM$patient[which(seuratM$dataset == "o294781_01-1_20220912_Hu_nucseq_Graz_21_HH_GEM")] <- "DH07"
seuratM$patient[which(seuratM$dataset == "o294781_02-2_20220912_Hu_nucseq_Graz_22_HH_GEM")] <- "DH08"
seuratM$patient[which(seuratM$dataset == "o294781_03-3_20220912_Hu_nucseq_Graz_23_HH_GEM")] <- "DH09"
seuratM$patient[which(seuratM$dataset == "o294781_04-4_20220912_Hu_nucseq_Graz_24_HH_GEM")] <- "DH10"

ordpatients <- c("DH01", "DH02", "DH03", "DH04", "DH05", "DH06", "DH07", "DH08", "DH09", "DH10", "CarTransPat01", "CarTransPat02", "CarTransPat03", "CarTransPat04", "CarTransPat05", "CarTransPat06", "CarTransPat07", "CarTransPat08","CarTransPat10", "CarTransPat11", "CarTransPat12", "CarTransPat13", "CarTransPat14", "CarTransPat15", "CarTransPat16", "CarTransPat18", "CarTransPat19", "CarTransPat24")

#Disease Condition
seuratM$diseaseCond <- "diseaseCond"
seuratM$diseaseCond[grepl("V1", seuratM$dataset)] <- "visit1"
seuratM$diseaseCond[grepl("V2|353921_12-12_20240515_Hu_nucseq_USZ_EMB010_V1_2", seuratM$dataset)] <- "visit2"
seuratM$diseaseCond[grepl("V3", seuratM$dataset)] <- "visit3"
seuratM$diseaseCond[grepl("V4", seuratM$dataset)] <- "visit4"
seuratM$diseaseCond[grepl("V5", seuratM$dataset)] <- "visit5"
seuratM$diseaseCond[grepl("VX1", seuratM$dataset)] <- "visitX1"
seuratM$diseaseCond[grepl("VX2", seuratM$dataset)] <- "visitX2"
seuratM$diseaseCond[grepl("VX3", seuratM$dataset)] <- "visitX3"
seuratM$diseaseCond[grepl("VX4", seuratM$dataset)] <- "visitX4"
seuratM$diseaseCond[grepl("HH", seuratM$dataset)] <- "donorheart"
seuratM$diseaseCond[grepl("RV|LV|expLV|expRV|331571_3-5_20231012_Hu_nucseq_USZ_HTx001|331571_4-6_20231012_Hu_nucseq_USZ_HTx002", seuratM$dataset)] <- "explant"
orddiseaseCond <- c("donorheart","visit1", "visit2" ,"visit3", "visit4", "visit5", "visitX1", "visitX2", "visitX3", "visitX4", "explant")

seuratM$patient_diseaseCond <- paste0(seuratM$patient, '_', seuratM$diseaseCond)
seuratM$patient_clusterName <- paste0(seuratM$patient, '_', seuratM$clusterName)
```

## abundance plots

``` r
order_keywords <- c("donorheart", 
                    "V1", "V2","V3",
                    "RV|LV|expLV|expRV|331571_3-5_20231012_Hu_nucseq_USZ_HTx001|331571_4-6_20231012_Hu_nucseq_USZ_HTx002")
files <- unique(seuratM$dataset)
ordered_files <- c()
for (key in order_keywords) {
  ordered_files <- c(ordered_files, files[grepl(key, files)])
}

###dataset
datList <- NULL
for(con in unique(seuratM$dataset)){
  seuratSub <- subset(seuratM, dataset==con)
  print(dim(seuratSub))
  dat_con <- as.data.frame(table(seuratSub$clusterName)) %>%
  mutate(percent=Freq/ncol(seuratSub)) %>% mutate(dataset=con)
  datList[[con]] <- dat_con
}
```

```
## [1] 45702  3400
## [1] 45702  2539
## [1] 45702  4528
## [1] 45702  4179
## [1] 45702  4217
## [1] 45702  3958
## [1] 45702  5294
## [1] 45702  2917
## [1] 45702  1894
## [1] 45702  2566
## [1] 45702  3532
## [1] 45702  2129
## [1] 45702  5295
## [1] 45702   509
## [1] 45702  3936
## [1] 45702  1132
## [1] 45702   493
## [1] 45702  1524
## [1] 45702  1941
## [1] 45702  1987
## [1] 45702  1688
## [1] 45702  1838
## [1] 45702   589
## [1] 45702    71
## [1] 45702   368
## [1] 45702   175
## [1] 45702   857
## [1] 45702  2057
## [1] 45702   258
## [1] 45702   403
## [1] 45702   193
## [1] 45702    55
## [1] 45702  2800
## [1] 45702  2581
## [1] 45702  1184
## [1] 45702  2976
## [1] 45702   418
## [1] 45702   638
## [1] 45702   719
## [1] 45702   863
## [1] 45702  3880
## [1] 45702  6662
## [1] 45702  2231
## [1] 45702  3143
## [1] 45702  1954
## [1] 45702  1433
## [1] 45702   857
## [1] 45702   245
## [1] 45702  1393
## [1] 45702  3138
## [1] 45702   826
## [1] 45702   743
## [1] 45702  1909
## [1] 45702  2673
## [1] 45702  1118
## [1] 45702   292
## [1] 45702  2227
## [1] 45702   672
## [1] 45702  2740
## [1] 45702  1491
## [1] 45702   541
## [1] 45702   139
## [1] 45702   401
## [1] 45702  2538
## [1] 45702  2606
## [1] 45702  1106
## [1] 45702   360
## [1] 45702  6256
## [1] 45702  7622
## [1] 45702   437
## [1] 45702   178
## [1] 45702   216
## [1] 45702   477
## [1] 45702   607
## [1] 45702   107
## [1] 45702   295
## [1] 45702  1044
## [1] 45702   453
## [1] 45702   431
## [1] 45702   432
## [1] 45702  1166
## [1] 45702   104
## [1] 45702   321
## [1] 45702   114
## [1] 45702   340
## [1] 45702   366
## [1] 45702    62
## [1] 45702   570
## [1] 45702   445
## [1] 45702    60
## [1] 45702   382
## [1] 45702    42
## [1] 45702   430
## [1] 45702   199
## [1] 45702  4026
## [1] 45702   294
## [1] 45702   480
## [1] 45702   706
## [1] 45702   907
## [1] 45702   472
## [1] 45702   408
## [1] 45702   642
## [1] 45702   737
## [1] 45702   186
## [1] 45702   402
## [1] 45702    64
## [1] 45702   542
## [1] 45702    93
## [1] 45702  4005
## [1] 45702  3922
## [1] 45702  4265
## [1] 45702  3853
## [1] 45702  6434
## [1] 45702 11568
## [1] 45702  1465
## [1] 45702  2064
## [1] 45702   866
## [1] 45702  2181
```

``` r
dat_all <- do.call("rbind", datList)

## plot abundance
ggbarplot(dat_all, x= "dataset", y= "percent", fill = "Var1", palette = colclusterName, legend = "right", legend.titel = "cluster", ylab = "frequency")  + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + scale_x_discrete(limits=ordered_files)
```

![](Abundances_files/figure-html/abundance dataset-1.png)<!-- -->


``` r
order_keywords <- c("donorheart", "visit1", "visit2", "visit3", "explant")
files <- unique(seuratM$patient_diseaseCond)
ordered_files <- c()
for (key in order_keywords) {
  ordered_files <- c(ordered_files, files[grepl(key, files)])
}

###patient_diseaseCond
datList <- NULL
for(con in unique(seuratM$patient_diseaseCond)){
  seuratSub <- subset(seuratM, patient_diseaseCond==con)
  print(dim(seuratSub))
  dat_con <- as.data.frame(table(seuratSub$clusterName)) %>%
  mutate(percent=Freq/ncol(seuratSub)) %>% mutate(patient_diseaseCond=con)
  datList[[con]] <- dat_con
}
```

```
## [1] 45702  5529
## [1] 45702  7067
## [1] 45702  8137
## [1] 45702  9511
## [1] 45702  4811
## [1] 45702  6098
## [1] 45702  9231
## [1] 45702   509
## [1] 45702  1132
## [1] 45702   493
## [1] 45702  1524
## [1] 45702  1941
## [1] 45702  1987
## [1] 45702  1688
## [1] 45702  1838
## [1] 45702   589
## [1] 45702    71
## [1] 45702   368
## [1] 45702   175
## [1] 45702   857
## [1] 45702  2057
## [1] 45702   258
## [1] 45702   403
## [1] 45702   193
## [1] 45702    55
## [1] 45702  5381
## [1] 45702  1184
## [1] 45702  2976
## [1] 45702   418
## [1] 45702   638
## [1] 45702   719
## [1] 45702   863
## [1] 45702 10542
## [1] 45702  5374
## [1] 45702  3387
## [1] 45702   857
## [1] 45702   245
## [1] 45702  4531
## [1] 45702   826
## [1] 45702   743
## [1] 45702  4582
## [1] 45702  1118
## [1] 45702   292
## [1] 45702  3718
## [1] 45702   672
## [1] 45702  2740
## [1] 45702   541
## [1] 45702   139
## [1] 45702   401
## [1] 45702  5144
## [1] 45702  1106
## [1] 45702   360
## [1] 45702 13878
## [1] 45702   437
## [1] 45702   178
## [1] 45702   216
## [1] 45702   477
## [1] 45702   607
## [1] 45702   107
## [1] 45702   295
## [1] 45702  1044
## [1] 45702   453
## [1] 45702   431
## [1] 45702   432
## [1] 45702  1166
## [1] 45702   104
## [1] 45702   321
## [1] 45702   114
## [1] 45702   340
## [1] 45702   366
## [1] 45702    62
## [1] 45702  1015
## [1] 45702    60
## [1] 45702   382
## [1] 45702    42
## [1] 45702   430
## [1] 45702   199
## [1] 45702  4026
## [1] 45702   294
## [1] 45702   480
## [1] 45702   706
## [1] 45702   907
## [1] 45702   472
## [1] 45702   408
## [1] 45702   642
## [1] 45702   737
## [1] 45702   186
## [1] 45702   402
## [1] 45702    64
## [1] 45702   542
## [1] 45702    93
## [1] 45702  4005
## [1] 45702  3922
## [1] 45702  4265
## [1] 45702  3853
## [1] 45702  6434
## [1] 45702 11568
## [1] 45702  1465
## [1] 45702  2064
## [1] 45702   866
## [1] 45702  2181
```

``` r
dat_all <- do.call("rbind", datList)

## plot abundance
ggbarplot(dat_all, x= "patient_diseaseCond", 
          y= "percent", fill = "Var1", 
          palette = colclusterName, legend = "right", 
          legend.titel = "cluster", 
          ylab = "frequency",
          )  + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + scale_x_discrete(limits=ordered_files)
```

![](Abundances_files/figure-html/abundance patient_diseaseCond-1.png)<!-- -->


``` r
orddiseaseCond <- c("donorheart","visit1", "visit2" ,"visit3","explant")

###diseaseCond
datList <- NULL
for(con in unique(seuratM$diseaseCond)){
  seuratSub <- subset(seuratM, diseaseCond==con)
  print(dim(seuratSub))
  dat_con <- as.data.frame(table(seuratSub$clusterName)) %>%
  mutate(percent=Freq/ncol(seuratSub)) %>% mutate(diseaseCond=con)
  datList[[con]] <- dat_con
}
```

```
## [1]  45702 107936
## [1] 45702 15328
## [1] 45702 14613
## [1] 45702  7921
## [1] 45702  2423
## [1] 45702  2784
## [1] 45702  5208
## [1] 45702   512
## [1] 45702  2007
## [1] 45702   402
## [1] 45702 40623
```

``` r
dat_all <- do.call("rbind", datList)

## plot abundance
ggbarplot(dat_all, x= "diseaseCond", y= "percent", fill = "Var1", palette = colclusterName, legend = "right", legend.titel = "cluster", ylab = "frequency")  + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + scale_x_discrete(limits=orddiseaseCond)
```

![](Abundances_files/figure-html/abundance diseaseCond-1.png)<!-- -->

## fractions according to patients&disease cond

``` r
##set order
seuratM$diseaseCond <- factor(seuratM$diseaseCond, levels=c("donorheart", "visit1", "visit2", "visit3", "explant"))

## 1. create data.frame with cluster counts per patient
## change "RNA_snn_res.0.25" to subset/cluster you're interested in ...
datFrac <- data.frame(table(seuratM$patient_diseaseCond, seuratM$clusterName))
colnames(datFrac) <- c("patient_diseaseCond", "subset", "cnt")

## 2. get total counts per patient to compute relative abundances from
## I added cond here as grouping variable for the plotting later ...
datSumPat <- data.frame(table(seuratM$patient_diseaseCond, seuratM$diseaseCond)) %>% 
  filter(Freq >0)
colnames(datSumPat) <- c("patient_diseaseCond", "diseaseCond", "cntPatTot")

## 3. join data.frames to compute rel abundances per patient
datFracSum <- datFrac %>% left_join(., datSumPat, by = "patient_diseaseCond") %>% 
  mutate(relCnt = cnt/cntPatTot)

## plot barplot with abundances for each subset grouped by cond
ggbarplot(datFracSum, x = "subset", y = "relCnt",
          fill = "diseaseCond", color = "diseaseCond",
          palette = coldiseaseCond,
          add = c("mean_se", "dotplot"),
          add.params = list(color="black", fill="diseaseCond", size=0.2),
          position = position_dodge(0.9),
          xlab = "subset",
          ylab = "relative abundance",
          legend = "right",
          legend.title = "") +
  rotate_x_text(angle = 90) 
```

![](Abundances_files/figure-html/fractions-1.png)<!-- -->

``` r
## plot barplot with abundances for individual subsets
clusterVec <- levels(seuratM)
createClusterPlot <- function(cluster) {
  datFracSumC <- datFracSum %>% filter(subset == cluster)

  ggbarplot(datFracSumC, x = "diseaseCond", y = "relCnt",
            fill = "diseaseCond", color = "diseaseCond",
            palette = coldiseaseCond,
            add = c("mean_se", "dotplot"),
            size = 5,
            add.params = list(color = "black", fill = "diseaseCond"),
            position = position_dodge(0.9),
            xlab = cluster,
            ylab = "relative abundance",
            legend = "right",
            legend.title = "") +
    stat_compare_means(method = "kruskal.test", label.y = 0.0)
}
lapply(clusterVec, createClusterPlot)
```

```
## [[1]]
```

![](Abundances_files/figure-html/fractions-2.png)<!-- -->

```
## 
## [[2]]
```

![](Abundances_files/figure-html/fractions-3.png)<!-- -->

```
## 
## [[3]]
```

![](Abundances_files/figure-html/fractions-4.png)<!-- -->

```
## 
## [[4]]
```

![](Abundances_files/figure-html/fractions-5.png)<!-- -->

```
## 
## [[5]]
```

![](Abundances_files/figure-html/fractions-6.png)<!-- -->

```
## 
## [[6]]
```

![](Abundances_files/figure-html/fractions-7.png)<!-- -->

```
## 
## [[7]]
```

![](Abundances_files/figure-html/fractions-8.png)<!-- -->

```
## 
## [[8]]
```

![](Abundances_files/figure-html/fractions-9.png)<!-- -->

```
## 
## [[9]]
```

![](Abundances_files/figure-html/fractions-10.png)<!-- -->

```
## 
## [[10]]
```

![](Abundances_files/figure-html/fractions-11.png)<!-- -->

```
## 
## [[11]]
```

![](Abundances_files/figure-html/fractions-12.png)<!-- -->

```
## 
## [[12]]
```

![](Abundances_files/figure-html/fractions-13.png)<!-- -->

```
## 
## [[13]]
```

![](Abundances_files/figure-html/fractions-14.png)<!-- -->

```
## 
## [[14]]
```

![](Abundances_files/figure-html/fractions-15.png)<!-- -->

```
## 
## [[15]]
```

![](Abundances_files/figure-html/fractions-16.png)<!-- -->

```
## 
## [[16]]
```

![](Abundances_files/figure-html/fractions-17.png)<!-- -->

```
## 
## [[17]]
```

![](Abundances_files/figure-html/fractions-18.png)<!-- -->

```
## 
## [[18]]
```

![](Abundances_files/figure-html/fractions-19.png)<!-- -->

```
## 
## [[19]]
```

![](Abundances_files/figure-html/fractions-20.png)<!-- -->

## session info

``` r
date()
```

```
## [1] "Mon Sep 29 19:53:29 2025"
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
## [16] destiny_3.22.0              circlize_0.4.16             muscat_1.22.0              
## [19] viridis_0.6.5               viridisLite_0.4.2           lubridate_1.9.4            
## [22] forcats_1.0.1               stringr_1.5.2               purrr_1.1.0                
## [25] readr_2.1.5                 tidyr_1.3.1                 tibble_3.3.0               
## [28] tidyverse_2.0.0             dplyr_1.1.4                 SingleCellExperiment_1.30.1
## [31] SummarizedExperiment_1.38.1 Biobase_2.68.0              GenomicRanges_1.60.0       
## [34] GenomeInfoDb_1.44.3         IRanges_2.42.0              S4Vectors_0.46.0           
## [37] BiocGenerics_0.54.0         generics_0.1.4              MatrixGenerics_1.20.0      
## [40] matrixStats_1.5.0           pheatmap_1.0.13             ggpubr_0.6.1               
## [43] ggplot2_4.0.0               Seurat_5.3.0                SeuratObject_5.2.0         
## [46] sp_2.2-0                   
## 
## loaded via a namespace (and not attached):
##   [1] igraph_2.1.4             ica_1.0-3                plotly_4.11.0           
##   [4] Formula_1.2-5            tidyselect_1.2.1         bit_4.6.0               
##   [7] doParallel_1.0.17        clue_0.3-66              lattice_0.22-7          
##  [10] rjson_0.2.23             blob_1.2.4               S4Arrays_1.8.1          
##  [13] pbkrtest_0.5.5           parallel_4.5.1           png_0.1-8               
##  [16] plotrix_3.8-4            cli_3.6.5                ggplotify_0.1.3         
##  [19] goftest_1.2-3            VIM_6.2.6                variancePartition_1.38.1
##  [22] BiocNeighbors_2.2.0      uwot_0.2.3               curl_7.0.0              
##  [25] mime_0.13                evaluate_1.0.5           tidytree_0.4.6          
##  [28] ComplexHeatmap_2.24.1    stringi_1.8.7            backports_1.5.0         
##  [31] lmerTest_3.1-3           qqconf_1.3.2             httpuv_1.6.16           
##  [34] magrittr_2.0.4           rappdirs_0.3.3           splines_4.5.1           
##  [37] sctransform_0.4.2        ggbeeswarm_0.7.2         DBI_1.2.3               
##  [40] jquerylib_0.1.4          smoother_1.3             withr_3.0.2             
##  [43] corpcor_1.6.10           reformulas_0.4.1         class_7.3-23            
##  [46] lmtest_0.9-40            formatR_1.14             htmlwidgets_1.6.4       
##  [49] fs_1.6.6                 ggrepel_0.9.6            labeling_0.4.3          
##  [52] fANCOVA_0.6-1            SparseArray_1.8.1        DESeq2_1.48.2           
##  [55] ranger_0.17.0            DEoptimR_1.1-4           reticulate_1.43.0       
##  [58] hexbin_1.28.5            zoo_1.8-14               XVector_0.48.0          
##  [61] knitr_1.50               ggplot.multistats_1.0.1  UCSC.utils_1.4.0        
##  [64] RhpcBLASctl_0.23-42      timechange_0.3.0         foreach_1.5.2           
##  [67] patchwork_1.3.2          caTools_1.18.3           data.table_1.17.8       
##  [70] ggtree_3.16.3            R.oo_1.27.1              RSpectra_0.16-2         
##  [73] irlba_2.3.5.1            fastDummies_1.7.5        gridGraphics_0.5-1      
##  [76] lazyeval_0.2.2           yaml_2.3.10              survival_3.8-3          
##  [79] scattermore_1.2          crayon_1.5.3             RcppAnnoy_0.0.22        
##  [82] RColorBrewer_1.1-3       progressr_0.16.0         later_1.4.4             
##  [85] ggridges_0.5.7           codetools_0.2-20         GlobalOptions_0.1.2     
##  [88] aod_1.3.3                KEGGREST_1.48.1          Rtsne_0.17              
##  [91] shape_1.4.6.1            limma_3.64.3             pkgconfig_2.0.3         
##  [94] TMB_1.9.17               spatstat.univar_3.1-4    mathjaxr_1.8-0          
##  [97] EnvStats_3.1.0           aplot_0.2.9              scatterplot3d_0.3-44    
## [100] spatstat.sparse_3.1-0    ape_5.8-1                xtable_1.8-4            
## [103] car_3.1-3                plyr_1.8.9               httr_1.4.7              
## [106] rbibutils_2.3            tools_4.5.1              globals_0.18.0          
## [109] beeswarm_0.4.0           broom_1.0.10             nlme_3.1-168            
## [112] lambda.r_1.2.4           assertthat_0.2.1         lme4_1.1-37             
## [115] digest_0.6.37            numDeriv_2016.8-1.1      Matrix_1.7-4            
## [118] farver_2.1.2             tzdb_0.5.0               remaCor_0.0.20          
## [121] reshape2_1.4.4           yulab.utils_0.2.1        glue_1.8.0              
## [124] cachem_1.1.0             polyclip_1.10-7          Biostrings_2.76.0       
## [127] mvtnorm_1.3-3            parallelly_1.45.1        mnormt_2.1.1            
## [130] statmod_1.5.0            RcppHNSW_0.6.0           ScaledMatrix_1.16.0     
## [133] carData_3.0-5            minqa_1.2.8              pbapply_1.7-4           
## [136] spam_2.11-1              gson_0.1.0               gtools_3.9.5            
## [139] ggsignif_0.6.4           RcppEigen_0.3.4.0.2      shiny_1.11.1            
## [142] GenomeInfoDbData_1.2.14  glmmTMB_1.1.12           R.utils_2.13.0          
## [145] memoise_2.0.1            rmarkdown_2.29           scales_1.4.0            
## [148] R.methodsS3_1.8.2        future_1.67.0            RANN_2.6.2              
## [151] spatstat.data_3.1-8      rstudioapi_0.17.1        cluster_2.1.8.1         
## [154] mutoss_0.1-13            spatstat.utils_3.2-0     hms_1.1.3               
## [157] fitdistrplus_1.2-4       cowplot_1.2.0            colorspace_2.1-2        
## [160] rlang_1.1.6              xts_0.14.1               dotCall64_1.2           
## [163] ggtangle_0.0.7           laeken_0.5.3             mgcv_1.9-3              
## [166] xfun_0.53                e1071_1.7-16             TH.data_1.1-4           
## [169] iterators_1.0.14         abind_1.4-8              GOSemSim_2.34.0         
## [172] treeio_1.32.0            futile.options_1.0.1     bitops_1.0-9            
## [175] Rdpack_2.6.4             promises_1.3.3           RSQLite_2.4.3           
## [178] qvalue_2.40.0            sandwich_3.1-1           fgsea_1.34.2            
## [181] DelayedArray_0.34.1      proxy_0.4-27             GO.db_3.21.0            
## [184] compiler_4.5.1           prettyunits_1.2.0        boot_1.3-32             
## [187] beachmat_2.24.0          listenv_0.9.1            Rcpp_1.1.0              
## [190] edgeR_4.6.3              BiocSingular_1.24.0      tensor_1.5.1            
## [193] MASS_7.3-65              progress_1.2.3           BiocParallel_1.42.2     
## [196] babelgene_22.9           spatstat.random_3.4-2    R6_2.6.1                
## [199] fastmap_1.2.0            multcomp_1.4-28          fastmatch_1.1-6         
## [202] rstatix_0.7.2            vipor_0.4.7              TTR_0.24.4              
## [205] ROCR_1.0-11              TFisher_0.2.0            rsvd_1.0.5              
## [208] vcd_1.4-13               nnet_7.3-20              gtable_0.3.6            
## [211] KernSmooth_2.23-26       miniUI_0.1.2             deldir_2.0-4            
## [214] htmltools_0.5.8.1        ggthemes_5.1.0           bit64_4.6.0-1           
## [217] spatstat.explore_3.5-3   lifecycle_1.0.4          blme_1.0-6              
## [220] S7_0.2.0                 nloptr_2.2.1             sass_0.4.10             
## [223] vctrs_0.6.5              robustbase_0.99-6        spatstat.geom_3.6-0     
## [226] sn_2.1.1                 ggfun_0.2.0              future.apply_1.20.0     
## [229] bslib_0.9.0              pillar_1.11.1            gplots_3.2.0            
## [232] pcaMethods_2.0.0         locfit_1.5-9.12          jsonlite_2.0.0          
## [235] GetoptLong_1.0.5
```
