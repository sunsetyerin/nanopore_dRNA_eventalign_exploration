
assign(".lib.loc", "/projects/yekim_prj/scratch/software/R/x86_64-pc-linux-gnu-library/4.1", envir = environment(.libPaths))

library("tidyverse")

# read f5c eventalign GGACT signal ####################################################################

cnot3_kd <- data.table::fread("/projects/ly_vu_direct_rna/MetaCompore/results/xpore/xpore_eventalign/cnot3_test_0_data_GGACT.tsv")
cnot3_kd2 <- data.table::fread("/projects/ly_vu_direct_rna/MetaCompore/results/xpore/xpore_eventalign/cnot3_test_1_data_GGACT.tsv")
control <- data.table::fread("/projects/ly_vu_direct_rna/MetaCompore/results/xpore/xpore_eventalign/control_0_data_GGACT.tsv")
control2 <- data.table::fread("/projects/ly_vu_direct_rna/MetaCompore/results/xpore/xpore_eventalign/control_1_data_GGACT.tsv")
mettl3_kd <- data.table::fread("/projects/ly_vu_direct_rna/MetaCompore/results/resquiggling/f5c_eventalign/test_0_data_GGACT.tsv")
mettl3_kd2 <- data.table::fread("/projects/ly_vu_direct_rna/MetaCompore/results/resquiggling/f5c_eventalign/test_1_data_GGACT.tsv")

colnames(cnot3_kd) <-   c("contig","position","reference_kmer","event_level_mean")
colnames(cnot3_kd2) <-  c("contig","position","reference_kmer","event_level_mean")
colnames(control) <-    c("contig","position","reference_kmer","event_level_mean")
colnames(control2) <-   c("contig","position","reference_kmer","event_level_mean")
colnames(mettl3_kd) <-  c("contig","position","reference_kmer","event_level_mean")
colnames(mettl3_kd2) <- c("contig","position","reference_kmer","event_level_mean")

# calculate density ####################################################################

D  <- density(cnot3_kd$event_level_mean)
D2 <- density(cnot3_kd2$event_level_mean)
D3 <- density(control$event_level_mean)
D4 <- density(control2$event_level_mean)
D5 <- density(mettl3_kd$event_level_mean)
D6 <- density(mettl3_kd2$event_level_mean)


pdf("f5c_eventalign_mean_level_density.pdf",8,8)
plot(D$x, D$y, 
     main="CNOT3 KDs vs controls vs METTL3 KDs",
     type='l', col="blue",ylab="density", xlab="event_level_mean")
lines(D2$x, D2$y, lty=1, col="blue")
lines(D3$x, D3$y, lty=1, col="red")
lines(D4$x, D4$y, lty=1, col="red")
lines(D5$x, D5$y, lty=1, col="green")
lines(D6$x, D6$y, lty=1, col="green")
legend("topright",
       c("CNOT3_KD","Control","METTL3_KD"),
       fill=c("blue","red","green"))
dev.off()

# scatter ####################################################################

cnot3_kd$Sample <- "CNOT3_KD:1"
cnot3_kd2$Sample <- "CNOT3_KD:2"
control$Sample <- "WT:1"
control2$Sample <- "WT:2"
mettl3_kd$Sample <- "METTL3_KD:1"
mettl3_kd2$Sample <- "METTL3_KD:2"

event <- rbind(cnot3_kd, cnot3_kd2, control, control2, mettl3_kd, mettl3_kd2)
event$index <- paste0(event$contig,";",event$position,";",event$reference_kmer)

event <- event %>% 
  group_by(Sample) %>% 
  separate(Sample, c("Condition", "Rep"),":") 

event2 <- event %>%
  reshape2::dcast(index+Condition~Rep, value.var="event_level_mean", fun.aggregate = median)

png("f5c_eventalign_scatter.png",type="cairo")
event2 %>% 
  ggplot(aes(x=`1`, y=`2`, colour=Condition)) + geom_point(size=0.5, alpha=0.3) + 
  geom_abline(slope=1, intercept=0) + 
  theme_bw(20) + 
  # scale_x_log10() +
  # scale_y_log10() +
  xlab("Replicate 1") + ylab("Replicate 2") + 
  geom_smooth(method=lm) +
  theme(legend.position="bottom")
dev.off()

# PCA sample ####################################################################

pca <- rbind(cnot3_kd, cnot3_kd2, control, control2, mettl3_kd, mettl3_kd2)
pca$index <- paste0(pca$contig,";",pca$position,";",pca$reference_kmer)
pca <- pca %>% select(c("index", "Sample", "event_level_mean"))
pca <- pca %>% reshape2::dcast(index ~ Sample, value.var="event_level_mean", fun.aggregate = median)

rownames(pca) <- pca$index
pca <- pca[,-1]
pca <- pca %>% drop_na() # 11150 rows > 9356 rows 

pca2 <- prcomp(t(as.matrix(pca)), center = T, scale=F)

## make a scree plot
pca2.var <- pca2$sdev^2
pca2.var.per <- round(pca2.var/sum(pca2.var)*100, 1)


## now make a fancy looking plot that shows the PCs and the variation:
pca2.data <- data.frame(SAMPLE_ID=rownames(pca2$x),
                           X=pca2$x[,1],
                           Y=pca2$x[,2])
pca2.data
str(pca2.data)
pca2.data$condition <- c("CNOT3_KD","CNOT3_KD","METTL3_KD","METTL3_KD","WT","WT")

pca3 <- 
  pca2.data %>% 
  ggplot(aes(x=X, y=Y, color = condition)) + #"f5c_eventalign_event_level_mean_pca.pdf"
  # ggplot(aes(x=X, y=Y, color = SAMPLE_ID)) + #"f5c_eventalign_event_level_mean_pca_replicate.pdf"
  geom_point() +
  xlab(paste("PC1 - ", pca2.var.per[1], "%", sep="")) +
  ylab(paste("PC2 - ", pca2.var.per[2], "%", sep="")) +
  theme_bw() +
  ggtitle("f5c eventalign event_level_mean PCA")

pdf(file = "f5c_eventalign_event_level_mean_pca.pdf",
    width = 5,
    height = 4)
pca3
dev.off()
