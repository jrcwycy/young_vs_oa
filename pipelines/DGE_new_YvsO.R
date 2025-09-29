library(BiocManager)
library("DESeq2")
library(dplyr)

# directory to HTSeq-count output files 
datadir <- "/scratch/indikar_root/indikar1/shared_data/chondroform/deseq_input/"

# create dataframe with file names
sampleFiles <- data.frame(fname=list.files(path=datadir, pattern="*.tsv"), stringsAsFactors = FALSE)

# add column for file basename (e.g. bc1 - bc15)
sampleFiles <- sampleFiles %>% transmute(sample=gsub("\\.counts\\.tsv", "", fname),fname)
head(sampleFiles)

# add sample info 
sampleFiles <- sampleFiles %>%
  mutate(Condition = c("Young", rep("OA", 2), rep("Young", 2), rep("OA", 2),
                       rep("Young", 2), rep("OA", 2), rep("Young", 2)),
         Time = c(rep("D28", 3), "D14", "D28", "D14", "D28", rep("D14", 4), rep("D28", 2)),
         Batch = c("batch1", rep("batch3", 2), rep("batch4", 4),
                   rep("batch2", 4), rep("batch3", 2)))

print(sampleFiles)

# with interaction term
dds <- DESeqDataSetFromHTSeqCount(sampleTable = sampleFiles, directory = datadir,
                                  design = ~ Batch + Time + Condition + Time:Condition)

dds

# without interaction term
dds <- DESeqDataSetFromHTSeqCount(sampleTable = sampleFiles, directory = datadir,
                                  design = ~ Batch + Time + Condition)

dds

# filter for genes with at least 10 counts in at least 3 samples
smallestGroupSize <- 3
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]

# set reference level as Young D14
dds$Condition <- relevel(dds$Condition, ref = "Young")
dds$Time <- relevel(dds$Time, ref = "D14")


dds <- DESeq(dds)
resultsNames(dds)

# Main effect of Condition
res <- results(dds, name = "Condition_OA_vs_Young", alpha = 0.05)

#res <- results(dds,
#               contrast=c("Condition", "OA", "Young"),
#               alpha=0.05)
res
summary(res)

# Main effect of Time
res_time <- results(dds, name = "Time_D28_vs_D14", alpha = 0.05)
summary(res_time)

# Interaction effect
res_int <- results(dds, name = "TimeD28.ConditionOA", alpha = 0.05)
summary(res_int)


# # Time effect for Young
# res_y_time <- results(dds,
#                       #contrast = list(c("Time_D28_vs_D14", "Time.Day28.ConditionYoung")),
#                       contrast = list(c("Time_D28_vs_D14")),
#                       alpha=0.05)
# res_y_time
# summary(res_y_time)
# 
# # Time effect for OA
# res_oa_time <- results(dds,
#                        contrast = list(c("Time_D28_vs_D14", "TimeD28.ConditionOA")),
#                        alpha=0.05)
# res_oa_time
# summary(res_oa_time)
# 
# # Interaction between condition and time
# res_int <- results(dds, name = "TimeD28.ConditionOA", alpha=0.05)
# res_int
# summary(res_int)


# Write results to dataframe for Python
res_df <- as.data.frame(res)
write.csv(as.data.frame(res_df), file="Results_Condition_OvsY.csv")

# for no interaction term data
# write.csv(as.data.frame(res_df), file="NewResults_OvsY_NoInt.csv")

res_df1 <- as.data.frame(res_time)
write.csv(as.data.frame(res_df1), file="Results_Time_D28vsD14.csv")

# res_df2 <- as.data.frame(res_oa_time)
# write.csv(as.data.frame(res_df2), file="NewResults_OvsO.csv")

res_df3 <- as.data.frame(res_int)
write.csv(as.data.frame(res_df3), file="Results_interaction.csv")


### Get normalized expression table ###
normalized_counts <- counts(dds, normalized = TRUE)
write.csv(as.data.frame(normalized_counts), file="NormCounts.csv")



### Check SD ###
rld <- rlog(dds, blind = FALSE) # regularized-log transform
library("vsn")
meanSdPlot(assay(rld))


### Plot heatmap of sample to sample distances ###
sampleDists <- dist(t(assay(rld)))

library("pheatmap")
library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(rld$Condition, rld$Time, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

### Plot PCA ###
library(ggplot2)

pcaData <- plotPCA(rld, intgroup=c("Condition", "Time"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
#write.csv(as.data.frame(pcaData), file = "pcaData_OvsY.csv")

ggplot(pcaData, aes(PC1, PC2, color=Condition, shape=Time)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() +
  scale_shape_manual(values = c(17, 19)) +
  scale_color_manual(values = c("#55a0fbff", "#ff8080ff")) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.line.x.top = element_line(colour = "black"),
    axis.line.y.right = element_line(colour = "black"),
    axis.ticks.length = unit(0.15, "cm"),
    axis.ticks.x.top = element_blank(),  # Remove ticks from the top axis
    axis.ticks.y.right = element_blank(),  # Remove ticks from the right axis
    axis.text.x.top = element_blank(),  # Remove text from the top axis
    axis.text.y.right = element_blank(),
    axis.ticks = element_line(color = "black"),
    axis.text.x = element_text(margin = margin(t = 6)),  # Increase space for x-axis tick labels
    axis.text.y = element_text(margin = margin(r = 6)),
    panel.background = element_rect(fill = "white", color = "white")
  ) +
  scale_x_continuous(sec.axis = dup_axis(name = NULL)) +  
  scale_y_continuous(sec.axis = dup_axis(name = NULL))

# Check batch effect
plotPCA(rld, intgroup=c("Batch"))

pcaData <- plotPCA(rld, intgroup=c("Batch"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
#write.csv(as.data.frame(pcaData), file = "pcaData_Batch.csv")

ggplot(pcaData, aes(PC1, PC2, shape=Batch)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() +
  scale_shape_manual(values = c(3, 7, 8, 6)) +
  #scale_color_manual(values = c("#55a0fbff", "#ff8080ff")) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.line.x.top = element_line(colour = "black"),
    axis.line.y.right = element_line(colour = "black"),
    axis.ticks.length = unit(0.15, "cm"),
    axis.ticks.x.top = element_blank(),  # Remove ticks from the top axis
    axis.ticks.y.right = element_blank(),  # Remove ticks from the right axis
    axis.text.x.top = element_blank(),  # Remove text from the top axis
    axis.text.y.right = element_blank(),
    axis.ticks = element_line(color = "black"),
    axis.text.x = element_text(margin = margin(t = 6)),  # Increase space for x-axis tick labels
    axis.text.y = element_text(margin = margin(r = 6)),
    panel.background = element_rect(fill = "white", color = "white")
  ) +
  scale_x_continuous(sec.axis = dup_axis(name = NULL)) +  
  scale_y_continuous(sec.axis = dup_axis(name = NULL))

# Get PCA data for Python plotting
allPcaData <- plotPCA(rld, intgroup=c("Condition", "Time", "Batch"), returnData=TRUE)
write.csv(as.data.frame(allPcaData), file = "pcaData_All.csv")

### Merge with gene_map to get gene_names ###
map_path <- "/home/jrcwycy/iChon2024/gene_map.csv"

gene_map <- read.csv(map_path, sep="\t", header = FALSE)
colnames(gene_map) <- c('gene_id', 'gene_name')
rownames(gene_map) <- gene_map$gene_id
gene_map <- gene_map %>% select(-gene_id)
head(gene_map)

res$gene_name <- gene_map[rownames(res), 'gene_name']
head(res)

res_time$gene_name <- gene_map[rownames(res_time), 'gene_name']
head(res_time)

res_int$gene_name <- gene_map[rownames(res_int), 'gene_name']
head(res_int)

### Write results to dataframes and save for Python analysis ###
# res_df <- as.data.frame(res)
# 
# write.csv(as.data.frame(res_df), file="DEGs_OvsY_results.csv")



### Volcano plot with top expressed DEGs
library("ggrepel")
res_df <- as.data.frame(res)

# Remove unnamed genes
res_df <- res_df %>% filter(!is.na(gene_name) & gene_name != "")

# Create custom color and point size for significant genes
res_df$color <- ifelse(res_df$padj < 0.01 & abs(res_df$log2FoldChange) > 1, "significant", "not_significant")
res_df$point_size <- ifelse(res_df$padj < 0.01 & abs(res_df$log2FoldChange) > 1, 2, 1.5)  # Larger size for significant points
head(res_df)

# Filter for significant genes
significant_genes <- res_df[res_df$color == "significant", ]

# Get top sig DE genes (both up and down regulated)
top_pos <- significant_genes %>% arrange(desc(log2FoldChange)) %>% head(10)
top_neg <- significant_genes %>% arrange(log2FoldChange) %>% head(10)

# Create plot
ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) +
  
  # Color points based on significance
  geom_point(aes(color = color, size = point_size)) +
  
  # Highlight the top positive and negative DE genes
  geom_point(data = top_pos, aes(x = log2FoldChange, y = -log10(padj)), color = "black", fill = "#a33a3aff", size = 2, shape = 21) +
  geom_point(data = top_neg, aes(x = log2FoldChange, y = -log10(padj)), color = "black", fill = "#a33a3aff", size = 2, shape = 21) +
  
  # Add text labels for the top DE genes
  geom_text_repel(data = top_pos, aes(x = log2FoldChange, y = -log10(padj), label = gene_name), color = "black", max.overlaps = 15) +
  geom_text_repel(data = top_neg, aes(x = log2FoldChange, y = -log10(padj), label = gene_name), color = "black", max.overlaps = 15) +
  
  xlim(c(-10, 10)) +
  #xlim(c(-6, 6)) +
  #xlim(c(-7, 7)) +
  
  scale_color_manual(values = c("significant" = "#e0a14eff", "not_significant" = "#5281b9ff")) +

  theme_bw() +
  theme(
    panel.border = element_rect(color="black", fill=NA, linewidth = 1),
    panel.grid = element_blank(),         
    axis.ticks = element_line(color="black"),
    axis.text = element_text(color="black", size = 10),
  ) +
  scale_size(range = c(1.5, 2))







