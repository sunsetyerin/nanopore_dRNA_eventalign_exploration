
assign(".lib.loc", "/projects/yekim_prj/scratch/software/R/x86_64-pc-linux-gnu-library/4.1", envir = environment(.libPaths))

library("tidyverse")
library("cluster") # for pam
library("Rtsne")
library("biomaRt")
library("dplyr")

####################################################################
# read f5c eventalign GGACT signal (control n=4) 
####################################################################

control <-  data.table::fread("/projects/ly_vu_direct_rna/MetaCompore/results/xpore/xpore_eventalign/GGACT/control_0_data_GGACT.tsv")
control2 <- data.table::fread("/projects/ly_vu_direct_rna/MetaCompore/results/xpore/xpore_eventalign/GGACT/control_1_data_GGACT.tsv")
control3 <- data.table::fread("/projects/ly_vu_direct_rna/MetaCompore/results/xpore/xpore_eventalign/GGACT/control_2_data_GGACT.tsv")
control4 <- data.table::fread("/projects/ly_vu_direct_rna/MetaCompore/results/xpore/xpore_eventalign/GGACT/control_3_data_GGACT.tsv")

colnames(control) <-   c("contig","position","reference_kmer","event_level_mean")
colnames(control2) <-  c("contig","position","reference_kmer","event_level_mean")
colnames(control3) <-  c("contig","position","reference_kmer","event_level_mean")
colnames(control4) <-  c("contig","position","reference_kmer","event_level_mean")

control$Sample <- "WT:1"
control2$Sample <- "WT:2"
control3$Sample <- "WT:3"
control4$Sample <- "WT:4"

# t-SNE ####################################################################

event <- rbind(control, control2, control3, control4)
event$index <- paste0(event$contig,";",event$position,";",event$reference_kmer)
event <- event %>% select(c("index", "Sample", "event_level_mean"))
event <- event %>% reshape2::dcast(index ~ Sample, value.var="event_level_mean", fun.aggregate = median)

rownames(event) <- event$index
event <- event[,-1]
event <- event %>% drop_na() # 10535 > 9640

# scale the data first prior to running t-SNE
event2 <- scale(event)

library("cluster") # for pam
cl <- pam(event2, 5) # PAM algorithm (Partitioning Around Medoids)

####################################################################

require(Rtsne)
set.seed(1234)

tsne <- Rtsne(event2, check_duplicates = F, pca = T, perplexity=30, theta=0, dims=2)

t.df <- as.data.frame(tsne$Y)
colnames(t.df) <- c("V1", "V2")
t.df <- cbind(t.df, Cluster = factor(cl$clustering))

png(file = "f5c_eventalign_event_level_mean_control_t-sne3.png", width = 1000, height = 1000)
ggplot(t.df, aes(x = V1, y = V2, color=Cluster)) + geom_point(size=1, alpha=0.5)
dev.off()

####################################################################

t.df$index <- rownames(t.df)
t.df <- t.df %>% separate(index, c("transcript", "pos","kmer"),";")

# Annotate gene names
ensembl = useMart("ensembl", dataset="hsapiens_gene_ensembl")
id2name <- getBM(attributes=c('ensembl_transcript_id_version', 'hgnc_symbol'), 
                 filters = 'ensembl_transcript_id_version', 
                 values = unique(t.df$transcript), 
                 mart = ensembl)

t.df <- left_join(t.df, id2name, by=c("transcript"="ensembl_transcript_id_version"))
# write.csv(t.df, "f5c_eventalign_event_level_mean_control_GGACT_PAM_cluster.csv",row.names = F)

####################################################################
# read f5c eventalign GGACT signal (METTL3 KD n=2)
####################################################################

mettl3_kd <-  data.table::fread("/projects/ly_vu_direct_rna/MetaCompore/results/resquiggling/f5c_eventalign/GGACT/test_0_data_GGACT.tsv")
mettl3_kd2 <- data.table::fread("/projects/ly_vu_direct_rna/MetaCompore/results/resquiggling/f5c_eventalign/GGACT/test_1_data_GGACT.tsv")

colnames(mettl3_kd) <-   c("contig","position","reference_kmer","event_level_mean")
colnames(mettl3_kd2) <-  c("contig","position","reference_kmer","event_level_mean")

mettl3_kd$Sample <- "METTL3_KD:1"
mettl3_kd2$Sample <- "METTL3_KD:2"

# t-SNE ####################################################################

event <- rbind(mettl3_kd, mettl3_kd2)
event$index <- paste0(event$contig,";",event$position,";",event$reference_kmer)
event <- event %>% dplyr::select(c("index", "Sample", "event_level_mean"))
event <- event %>% reshape2::dcast(index ~ Sample, value.var="event_level_mean", fun.aggregate = median)

rownames(event) <- event$index
event <- event[,-1]
event <- event %>% drop_na() # 10738 > 10208

# scale the data first prior to running t-SNE
event2 <- scale(event)

library("cluster") # for pam
cl2 <- pam(event2, 5) # PAM algorithm (Partitioning Around Medoids)

####################################################################

require(Rtsne)
set.seed(1234)

tsne2 <- Rtsne(event2, check_duplicates = F, pca = T, perplexity=30, theta=0, dims=2)

t.df <- as.data.frame(tsne2$Y)
colnames(t.df) <- c("V1", "V2")
t.df <- cbind(t.df, Cluster = factor(cl2$clustering))

png(file = "f5c_eventalign_event_level_mean_mettl3KD_t-sne2.png", width = 800, height = 800)
ggplot(t.df, aes(x = V1, y = V2, color=Cluster)) + geom_point(size=0.5, alpha=0.5)
dev.off()

####################################################################

t.df$index <- rownames(t.df)
t.df <- t.df %>% separate(index, c("transcript", "pos","kmer"),";")

# Annotate gene names
ensembl = useMart("ensembl", dataset="hsapiens_gene_ensembl")
id2name <- getBM(attributes=c('ensembl_transcript_id_version', 'hgnc_symbol'), 
                 filters = 'ensembl_transcript_id_version', 
                 values = unique(t.df$transcript), 
                 mart = ensembl)

t.df <- left_join(t.df, id2name, by=c("transcript"="ensembl_transcript_id_version"))
# write.csv(t.df, "f5c_eventalign_event_level_mean_mettl3KD_GGACT_PAM_cluster.csv",row.names = F)

####################################################################
# read f5c eventalign GGACT signal (control n=4) 
####################################################################

cnot3_kd <-  data.table::fread("/projects/ly_vu_direct_rna/MetaCompore/results/xpore/xpore_eventalign/GGACT/cnot3_test_0_data_GGACT.tsv")
cnot3_kd2 <-  data.table::fread("/projects/ly_vu_direct_rna/MetaCompore/results/xpore/xpore_eventalign/GGACT/cnot3_test_1_data_GGACT.tsv")
cnot3_kd3 <-  data.table::fread("/projects/ly_vu_direct_rna/MetaCompore/results/xpore/xpore_eventalign/GGACT/cnot3_test_2_data_GGACT.tsv")
cnot3_kd4 <-  data.table::fread("/projects/ly_vu_direct_rna/MetaCompore/results/xpore/xpore_eventalign/GGACT/cnot3_test_3_data_GGACT.tsv")

colnames(cnot3_kd) <-   c("contig","position","reference_kmer","event_level_mean")
colnames(cnot3_kd2) <-  c("contig","position","reference_kmer","event_level_mean")
colnames(cnot3_kd3) <-  c("contig","position","reference_kmer","event_level_mean")
colnames(cnot3_kd4) <-  c("contig","position","reference_kmer","event_level_mean")

cnot3_kd$Sample <- "CNOT3_KD:1"
cnot3_kd2$Sample <- "CNOT3_KD:2"
cnot3_kd3$Sample <- "CNOT3_KD:3"
cnot3_kd4$Sample <- "CNOT3_KD:4"

# t-SNE ####################################################################

event <- rbind(cnot3_kd, cnot3_kd2, cnot3_kd3, cnot3_kd4)
event$index <- paste0(event$contig,";",event$position,";",event$reference_kmer)
event <- event %>% dplyr::select(c("index", "Sample", "event_level_mean"))
event <- event %>% reshape2::dcast(index ~ Sample, value.var="event_level_mean", fun.aggregate = median)

rownames(event) <- event$index
event <- event[,-1]
event <- event %>% drop_na() # 10577 > 9783

# scale the data first prior to running t-SNE
event2 <- scale(event)

library("cluster") # for pam
cl3 <- pam(event2, 5) # PAM algorithm (Partitioning Around Medoids)

####################################################################

require(Rtsne)
set.seed(1234)

tsne3 <- Rtsne(event2, check_duplicates = F, pca = T, perplexity=30, theta=0.5, dims=2)

t.df <- as.data.frame(tsne3$Y)
colnames(t.df) <- c("V1", "V2")
t.df <- cbind(t.df, Cluster = factor(cl3$clustering))

png(file = "f5c_eventalign_event_level_mean_cnot3kd_t-sne2.png", width = 800, height = 800)
ggplot(t.df, aes(x = V1, y = V2, color=Cluster)) + geom_point(size=0.5, alpha=0.5)
dev.off()

####################################################################

t.df$index <- rownames(t.df)
t.df <- t.df %>% separate(index, c("transcript", "pos","kmer"),";")

# Annotate gene names
ensembl = useMart("ensembl", dataset="hsapiens_gene_ensembl")
id2name <- getBM(attributes=c('ensembl_transcript_id_version', 'hgnc_symbol'), 
                 filters = 'ensembl_transcript_id_version', 
                 values = unique(t.df$transcript), 
                 mart = ensembl)

t.df <- left_join(t.df, id2name, by=c("transcript"="ensembl_transcript_id_version"))
# write.csv(t.df, "f5c_eventalign_event_level_mean_cnot3kd_GGACT_PAM_cluster.csv",row.names = F)
