library(org.Hs.eg.db)
library(enrichplot)#GO,KEGG,GSEA
library(clusterProfiler)#GO,KEGG,GSEA
library(GOplot)#弦图，弦表图，系统聚类图
library(DOSE)
library(ggnewscale)
library(topGO)#绘制通路网络图
library(circlize)#绘制富集分析圈图
library(ComplexHeatmap)#绘制图例
GO_database <- 'org.Hs.eg.db'
KEGG_database <- 'hsa'

setwd("D:/20230315")
a <- read.table("miRanda_DF2kinds_final_res.txt",sep = "\t")
a1 <- a[,2]
b <- c(1,2,3)
for(i in 1:length(a1)){b[i] <- strsplit(a1[i],"\\|")[[1]][4]}
b2 <- unique(b)
targetgene <- b2
SUB <- bitr(targetgene,fromType = 'SYMBOL',toType = 'ENSEMBL',OrgDb = GO_database)
SUB2 <- targetgene[grep("ENSG",targetgene)]
SUB2L <- data.frame(SYMBOL=SUB2,ENSEMBL=SUB2)
SUBAL <- rbind(SUB,SUB2L)
SUBAL <- SUBAL[!duplicated(SUBAL$SYMBOL),]#17317
SUBAL <- SUBAL[!duplicated(SUBAL$ENSEMBL),]#17312
write.csv(SUBAL,"3Dgenetaget_to_ENSEMBL_2kinds.csv",row.names = F)


GO<-enrichGO( SUBAL$ENSEMBL,#GO富集分析
              OrgDb = GO_database,
              keyType = "ENSEMBL",#设定读取的gene ID类型
              ont = "ALL",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
              pvalueCutoff = 0.05,#设定p值阈值
              qvalueCutoff = 0.05,#设定q值阈值
              readable = T)
library(ggplot2)
barplot(GO, split="ONTOLOGY")+facet_grid(ONTOLOGY~., scale="free")#柱状图

dotplot(GO, split="ONTOLOGY",label_format=100,font.size=25)+facet_grid(ONTOLOGY~., scale="free")+theme(legend.text = element_text(size=20),legend.title=element_text(size=25))#点状图
ggsave(filename = "3Dgene_GOdotplot_2kinds.pdf",width=25,height = 15,units = "in",dpi=600)
enrichplot::cnetplot(GO,circular=FALSE,colorEdge = TRUE)


#KEGG
gene <- bitr(targetgene,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
gene2 <- bitr(targetgene,fromType = 'ENSEMBL',toType = 'ENTREZID',OrgDb = GO_database)
colnames(gene2)[1] <- "SYMBOL"
GENE <- rbind(gene,gene2)#16998
KEGG<-enrichKEGG(GENE$ENTREZID,#KEGG富集分析
                 organism = KEGG_database,
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 0.05)

barplot(KEGG,showCategory = 40,title = 'KEGG Pathway',label_format=100)#label_format=100可以避免标签重叠
ggsave(filename = "3Dgene_KEGGbarplot_2kinds.pdf",width=25,height = 15,units = "in",dpi=600)
dotplot(GO, split="ONTOLOGY",label_format=100)+facet_grid(ONTOLOGY~., scale="free")#点状图
dotplot(KEGG,showCategory = 40,label_format=100)
enrichplot::cnetplot(GO,circular=FALSE,colorEdge = TRUE)
enrichplot::cnetplot(KEGG,circular=FALSE,colorEdge = TRUE)

GO2 <- pairwise_termsim(GO)
KEGG2 <- pairwise_termsim(KEGG)
enrichplot::emapplot(GO2,showCategory = 50, color = "p.adjust", layout = "kk")#通路间关联网络图
enrichplot::emapplot(KEGG2,showCategory =50, color = "p.adjust", layout = "kk")
write.table(KEGG$ID, file = "KEGG_IDs.txt", #将所有KEGG富集到的通路写入本地文件查看
            append = FALSE, quote = TRUE, sep = " ",
            eol = "\n", na = "NA", dec = ".", row.names = TRUE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")
# browseKEGG(KEGG,"hsa05166")#选择其中的hsa05166通路进行展示



#25kinds
rm(list=ls())
GO_database <- 'org.Hs.eg.db'
KEGG_database <- 'hsa'
a <- read.table("miRanda_DF25kinds_final_res.txt",sep = "\t")
a1 <- a[,2]
b <- c(1,2,3)
for(i in 1:length(a1)){b[i] <- strsplit(a1[i],"\\|")[[1]][4]}
b2 <- unique(b)
targetgene <- b2
SUB <- bitr(targetgene,fromType = 'SYMBOL',toType = 'ENSEMBL',OrgDb = GO_database)
SUB2 <- targetgene[grep("ENSG",targetgene)]
SUB2L <- data.frame(SYMBOL=SUB2,ENSEMBL=SUB2)
SUBAL <- rbind(SUB,SUB2L)
SUBAL <- SUBAL[!duplicated(SUBAL$SYMBOL),]#17317
SUBAL <- SUBAL[!duplicated(SUBAL$ENSEMBL),]#17312
write.csv(SUBAL,"3Dgenetaget_to_ENSEMBL_25kinds.csv",row.names = F)


GO<-enrichGO( SUBAL$ENSEMBL,#GO富集分析
              OrgDb = GO_database,
              keyType = "ENSEMBL",#设定读取的gene ID类型
              ont = "ALL",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
              pvalueCutoff = 0.05,#设定p值阈值
              qvalueCutoff = 0.05,#设定q值阈值
              readable = T)
library(ggplot2)

dotplot(GO, split="ONTOLOGY",label_format=100,font.size=25)+facet_grid(ONTOLOGY~., scale="free")+theme(legend.text = element_text(size=20),legend.title=element_text(size=25))#点状图
ggsave(filename = "3Dgene_GOdotplot_25kinds.pdf",width=25,height = 15,units = "in",dpi=600)



#KEGG
gene <- bitr(targetgene,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
gene2 <- bitr(targetgene,fromType = 'ENSEMBL',toType = 'ENTREZID',OrgDb = GO_database)
colnames(gene2)[1] <- "SYMBOL"
GENE <- rbind(gene,gene2)#16985
KEGG<-enrichKEGG(GENE$ENTREZID,#KEGG富集分析
                 organism = KEGG_database,
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 0.05)

barplot(KEGG,showCategory = 40,title = 'KEGG Pathway',label_format=100)#label_format=100可以避免标签重叠
ggsave(filename = "3Dgene_KEGGbarplot_25kinds.pdf",width=25,height = 15,units = "in",dpi=600)
