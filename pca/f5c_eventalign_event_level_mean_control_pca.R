
assign(".lib.loc", "/projects/yekim_prj/scratch/software/R/x86_64-pc-linux-gnu-library/4.1", envir = environment(.libPaths))

library("tidyverse")

# read f5c eventalign GGACT signal ####################################################################

control <- data.table::fread("/projects/ly_vu_direct_rna/MetaCompore/results/xpore/xpore_eventalign/control_0_data_GGACT.tsv")
control2 <- data.table::fread("/projects/ly_vu_direct_rna/MetaCompore/results/xpore/xpore_eventalign/control_1_data_GGACT.tsv")

colnames(control) <-    c("contig","position","reference_kmer","event_level_mean")
colnames(control2) <-   c("contig","position","reference_kmer","event_level_mean")

control$Sample <- "WT:1"
control2$Sample <- "WT:2"

# PCA sample ####################################################################

pca <- rbind(control, control2)
pca$index <- paste0(pca$contig,";",pca$position,";",pca$reference_kmer)
pca <- pca %>% select(c("index", "Sample", "event_level_mean"))
pca <- pca %>% reshape2::dcast(index ~ Sample, value.var="event_level_mean", fun.aggregate = median)

rownames(pca) <- pca$index
pca <- pca[,-1]
pca <- pca %>% drop_na() # 10379 rows > 9768 rows 

pca2 <- prcomp((as.matrix(pca)), center = T, scale=F)

## make a scree plot
pca2.var <- pca2$sdev^2
pca2.var.per <- round(pca2.var/sum(pca2.var)*100, 1)


## now make a fancy looking plot that shows the PCs and the variation:
pca2.data <- data.frame(ENST=substr(rownames(pca2$x),0,15),
                        X=pca2$x[,1],
                        Y=pca2$x[,2])
pca2.data
str(pca2.data)

pca3 <- 
  pca2.data %>% 
  ggplot(aes(x=X, y=Y)) + 
  geom_point() +
  xlab(paste("PC1 - ", pca2.var.per[1], "%", sep="")) +
  ylab(paste("PC2 - ", pca2.var.per[2], "%", sep="")) +
  theme_bw() +
  ggtitle("f5c eventalign event_level_mean PCA")

pdf(file = "f5c_eventalign_event_level_mean_control_pca.pdf",
    width = 5,
    height = 4)
pca3
dev.off()
