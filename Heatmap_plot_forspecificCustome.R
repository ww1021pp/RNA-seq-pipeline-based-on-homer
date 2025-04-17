rm(list=ls())
setwd("~/Documents/DY_RNA_Seq_4_8_19/HFD_B6S129_metabolism/1.DE_res04162025/")

library(openxlsx)
library(dplyr)
library(limma)

data=read.xlsx("~/Documents/DY_RNA_Seq_4_8_19/HFD_B6S129_metabolism/UTF-8B6_129_NC_HFD_ZT10_Metabolites.xlsx")
sample.meta =  read.delim("~/Documents/DY_RNA_Seq_4_8_19/HFD_B6S129_metabolism/sampleInfo.txt") 

colnames(data)[1]="Metabolism_Name"


metabolite_count = data
metabolite_count <- tibble::column_to_rownames(metabolite_count,var="Metabolism_Name")
metabolite_count[1:ncol(metabolite_count)] <- sapply(metabolite_count[1:ncol(metabolite_count)],as.numeric)
col_name=gsub("129","Sv",gsub("-",".",names(metabolite_count)))
names(metabolite_count)=col_name
sample.meta = sample.meta 

##range(na.omit(metabolite_count))
##[1] 4.256370e-03 5.312282e+04
log_matrix <- log2(metabolite_count+0.01)
norm_matrix=log_matrix
raw_df=cbind(data$Metabolism_Name,metabolite_count)
colnames(raw_df)[1]="Metabolism"

##Missing values: 1.55 %
#####PCA analysis of all samples###
library(ggplot2)
# Transpose so samples are rows
pca <- prcomp(t(na.omit(norm_matrix)), scale. = TRUE)
# Create a data frame with PCA results and metadata
pca_df <- as.data.frame(pca$x)
pca_df$group <- sample.meta$Group
pca_df$sample <- sample.meta$Sample_ID
pca_df$Strain = gsub("Esrr.*","",sample.meta$Group)
# Plot PCA
ggplot(pca_df, aes(PC1, PC2, color = group, label = sample,shape = Strain)) +
  geom_point(size = 3) +
  geom_text(vjust = 1.5, size = 3, check_overlap = TRUE) +
  labs(title = "PCA of Metabolomics Data") +
  theme_bw() +
  theme(axis.line = element_blank(),  
        axis.text = element_text(color = "black", size = 8),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5), 
        axis.title = element_text(size = 8,colour = "black"),
        plot.title = element_text(size = 8,colour = "black"),
  )
ggsave("allSamples_PCA.pdf",width = 8,height = 8)

####all samples fo DE ZT10 vs ZT22#########
cons <- c("B6_NC","B6_HF","Sv_NC","Sv_HF")

sample.meta$group <- factor(gsub("129","Sv",gsub("-",".",sample.meta$Group)))
sample.meta$group2 <- factor(c(rep("B6",6),rep("Sv",5)))
design <- model.matrix(~0 + group, data = sample.meta)
colnames(design) <- levels(sample.meta$group)
design2 <-  model.matrix(~0 + group2, data = sample.meta)
colnames(design2) <- levels(sample.meta$group2)

# Fit the model
fit <- lmFit(norm_matrix, design)

# All pairwise comparisons
#combos <- combn(levels(sample.meta$group), 2, simplify = FALSE)

combos <- list(
  c("B6_HF","B6_NC"),
  c("Sv_NC","B6_NC"),
  c("Sv_HF","Sv_NC"),#### g2 is control group, while g1 is Case group
  c("Sv_HF","B6_HF"),
  c("Sv","B6")
  
  )

####rowmean of each condition###
cond_rawRowMean <- lapply(cons, function(x) {
  cols <- grep(x, colnames(raw_df), value = TRUE)
  rowMeans(raw_df[, cols], na.rm = TRUE)
})
cond_rawRowMean_Df <- as.data.frame(cond_rawRowMean)
colnames(cond_rawRowMean_Df) <- paste0(cons,".Mean")



results_allSam=list()
for (i in 1:4) {
  g1 <- combos[[i]][1]
  g2 <- combos[[i]][2]
  
  contrast <- makeContrasts(contrasts = paste0(g1, "-", g2), levels = design)
  fit2 <- contrasts.fit(fit, contrast)
  fit2 <- eBayes(fit2)
  result <- topTable(fit2, number = Inf, adjust = "fdr") %>% select(logFC,P.Value,adj.P.Val)
  result = result[rownames(raw_df),]
  # tmp=gsub("Esrrg","",paste0(g1,"-vs-",g2,"_",colnames(result)))
  #  names(result)=tmp
  results_allSam[[paste0(g1, "-", g2)]] = result
}
saveRDS(results_allSam,"DEanalysis.rds")

results= do.call(cbind,results_allSam)

res_df <- cbind(metabolite_meta[rownames(results),c(1:4,15:22)],results,cond_normRowMean_Df[rownames(results),],metabolite_count[rownames(results),])
res_df[is.na(res_df)] = NA
write.xlsx(res_df,"DE_meatabolism.remove1Sample_andcountable0415.xlsx",keepNA = TRUE)


B6_DEM <- results_allSam$`B6_HF-B6_NC` %>% filter(P.Value <0.05 & abs(logFC) >=log2(1.5))
Sv_DEM <- results_allSam$`Sv_HF-Sv_NC` %>% filter(P.Value <0.05 & abs(logFC) >=log2(1.5))


common <- intersect(rownames(B6_DEM),rownames(Sv_DEM))
B6_spe <- setdiff(rownames(B6_DEM),rownames(Sv_DEM))
S129_spe <- setdiff(rownames(Sv_DEM),rownames(B6_DEM))

res_list=list(
  B6_spe=B6_spe,common=common,
  Sv129_spe=S129_spe
)
openxlsx::write.xlsx(res_list,"ZT10_HFDvsNC_DEanalysis.xlsx")




########Strain###

fit11 <- lmFit(norm_matrix,design = design2)
i=5
g1 <- combos[[i]][1]
g2 <- combos[[i]][2]

contrast <- makeContrasts(contrasts = paste0(g1, "-", g2), levels = design2)
fit2 <- contrasts.fit(fit11, contrast)
fit2 <- eBayes(fit2)
result <- topTable(fit2, number = Inf, adjust = "fdr") %>% select(logFC,P.Value,adj.P.Val)
DE_result <- result %>% filter(P.Value <0.05 & abs(logFC) >=log2(1.5))
tmp=setdiff(rownames(DE_result),common)
results_allSam[[paste0(g1, "-", g2)]]=result[rownames(raw_df),]

write.table(tmp,"Strain_specific_DE.txt",row.names = F,col.names = F,quote = F)

DEList <- list(
  B6_DEM = results_allSam$`B6_HF-B6_NC` %>% filter(P.Value < 0.05 & abs(logFC)> log2(1.5)), 
  Sv_DEM = results_allSam$`Sv_HF-Sv_NC` %>% filter(P.Value < 0.05 & abs(logFC)> log2(1.5)),
  Strain_NC = results_allSam$`Sv_NC-B6_NC` %>% filter(P.Value < 0.05 & abs(logFC)> log2(1.5)),
  Strain_HF = results_allSam$`Sv_HF-B6_HF` %>% filter(P.Value < 0.05 & abs(logFC)> log2(1.5)),
  Sv_B6 = results_allSam$`Sv-B6` %>% filter(P.Value < 0.05 & abs(logFC)> log2(1.5))
)

tmp = list(
  B6_DEM=rownames(DEList$B6_DEM),
  Sv_DEM= rownames(DEList$Sv_DEM),
  Strain_NC = rownames(DEList$Strain_NC),
  Strain_HF = rownames(DEList$Strain_HF)
)

library(VennDiagram);
venn.diagram(
  x = list(
    B6_DEM=rownames(DEList$B6_DEM),
    Sv_DEM= rownames(DEList$Sv_DEM),
    Strain_NC = rownames(DEList$Strain_NC),
    Strain_HF = rownames(DEList$Strain_HF)
  ),
  filename = "ZT10.HFDvsNC.Strain.pdf",
  col = "black",
  lty = "solid",
  lwd = 4,
  fill = c("cornflowerblue", "green", "yellow", "darkorchid1"),
  alpha = 0.50,
  label.col = c("orange", "white", "darkorchid4", "white", "white", "white", "white", "white", "darkblue", "white", "white", "white", "white", "darkgreen", "white"),
  cex = 2.5,
  cat.pos = c(-20, 20, 20, -20),
  #cex=0, ##for num label
  fontfamily = "serif",
  fontface = "bold",
  cat.col = c("darkblue", "darkgreen", "orange", "darkorchid4"),
  cat.cex = 1.1,
  ##cat.cex=0, ###for catogory
  cat.fontfamily = "serif"
)


tmp = list(
  B6_DEM=rownames(DEList$B6_DEM),
  Sv_DEM= rownames(DEList$Sv_DEM),
  Sv_B6 = rownames(DEList$Sv_B6)
)

diet_common <- setdiff(intersect(tmp$B6_DEM,tmp$Sv_DEM),tmp$Sv_B6) 
B6_spe <- setdiff(setdiff(tmp$B6_DEM,tmp$Sv_DEM),tmp$Sv_B6)
Sv_spe <- setdiff(setdiff(tmp$Sv_DEM,tmp$B6_DEM),tmp$Sv_B6)
Strain_induced <- setdiff(setdiff(tmp$Sv_B6,tmp$B6_DEM),tmp$Sv_DEM)



order_B6spe <- raw_df[B6_spe,2:12]
order_common <- raw_df[diet_common,2:12]
order_Svspe <- raw_df[Sv_spe,2:12]
order_strain <- raw_df[Strain_induced,2:12]




library(pheatmap)
bk_FC = unique(c(seq(-1,1,length=20)))
hmcols_FC<- colorRampPalette(c("blue","black","red"))(length(bk_FC)-1)
#hmcols_FC<-  c(colorRampPalette(colors = c("blue", "white","red"))(length(bk_FC)/2), colorRampPalette(colors = c("white", "white"))(length(bk_FC)/2))

pdf("rawdata_heatmapHFDvsNC_B6Spe.pdf", width=4, height = 2) 
p1_B6spe <- pheatmap(as.matrix(order_B6spe), cluster_rows = T, clustering_distance_rows = "correlation", 
                     col= hmcols_FC,
                     breaks = bk_FC,cluster_cols = F, legend=TRUE, show_rownames=FALSE,scale ="row")
dev.off()

pdf("rawdata_heatmapHFDvsNC_dietCommon.pdf", width=4, height = 6) 
p1_common<-pheatmap(as.matrix(order_common), cluster_rows = T, clustering_distance_rows = "correlation", 
         col= hmcols_FC,
         breaks = bk_FC,cluster_cols = F, legend=TRUE, show_rownames=FALSE,scale ="row")
dev.off()

pdf("rawdata_heatmapHFDvsNC_SvSpe.pdf", width=4, height = 6) 
p1_Svspe <- pheatmap(as.matrix(order_Svspe), cluster_rows = T, clustering_distance_rows = "correlation", 
                     col= hmcols_FC,
                     breaks = bk_FC,cluster_cols = F, legend=TRUE, show_rownames=FALSE,scale ="row")
dev.off()


pdf("rawdata_heatmapHFDvsNC_Strain.pdf", width=4, height = 6) 
p1_strain <- pheatmap(as.matrix(order_strain), cluster_rows = T, clustering_distance_rows = "correlation", 
                     col= hmcols_FC,
                     breaks = bk_FC,cluster_cols = F, legend=TRUE, show_rownames=FALSE,scale ="row")
dev.off()


data_last <- rbind(
                   order_B6spe[rownames(order_B6spe)[p1_B6spe$tree_row$order],],
                   order_common[rownames(order_common)[p1_common$tree_row$order],],
                   order_Svspe[rownames(order_Svspe)[p1_Svspe$tree_row$order],],
                   order_strain[rownames(order_strain)[p1_strain$tree_row$order],]
)
                   
pdf("rawdata_heatmapHFDvsNC.pdf", width=4, height = 8) 
pheatmap(as.matrix(data_last), cluster_rows = F, clustering_distance_rows = "correlation", 
                      col= hmcols_FC,
                      gaps_row = c(45,130,199),
                      breaks = bk_FC,cluster_cols = F, legend=TRUE, show_rownames=FALSE,scale ="row")
dev.off()



