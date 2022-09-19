library(DESeq2)
library(ggplot2)
library(EnhancedVolcano)
library(EDASeq)
library(ggrepel)
library(org.Hs.eg.db)
library(clusterProfiler)
library(ggpubr)
library(openxlsx)
setwd("~/TFM/r/mapping")


#### macrophage DATA

### controls
macrophage_controls<-NULL
for (f in c("M-Blank-I-3hr_ReadsPerGene.out.tab", "M-Blank-II-3hr_ReadsPerGene.out.tab", 
            "M-Blank-I-24hr_ReadsPerGene.out.tab", "M-Blank-II-24hr_ReadsPerGene.out.tab"))
{
  file<-read.table(f, header = TRUE, sep = "\t")
  macrophage_controls<-cbind(macrophage_controls,file[,2])
}

rownames(macrophage_controls)<-file[,1]
colnames(macrophage_controls)<-c("M-Blank-I-3hr", "M-Blank-II-3hr", "M-Blank-I-24hr", "M-Blank-II-24hr")
macrophage_controls<-macrophage_controls[-c(1,2,3),]


### infection

infection<-NULL
for (f in c("M-A-I-3hr_ReadsPerGene.out.tab", "M-B-I-3hr_ReadsPerGene.out.tab", "M-C-I-3hr_ReadsPerGene.out.tab", "M-D-I-3hr_ReadsPerGene.out.tab",
            "M-A-II-3hr_ReadsPerGene.out.tab", "M-B-II-3hr_ReadsPerGene.out.tab", "M-C-II-3hr_ReadsPerGene.out.tab", "M-D-II-3hr_ReadsPerGene.out.tab",
            "M-A-III-3hr_ReadsPerGene.out.tab", "M-B-III-3hr_ReadsPerGene.out.tab", "M-C-III-3hr_ReadsPerGene.out.tab", "M-D-III-3hr_ReadsPerGene.out.tab",
            "M-A-I-24hr_ReadsPerGene.out.tab", "M-B-I-24hr_ReadsPerGene.out.tab", "M-C-I-24hr_ReadsPerGene.out.tab", "M-D-I-24hr_ReadsPerGene.out.tab",
            "M-A-II-24hr_ReadsPerGene.out.tab", "M-B-II-24hr_ReadsPerGene.out.tab", "M-C-II-24hr_ReadsPerGene.out.tab", "M-D-II-24hr_ReadsPerGene.out.tab", 
            "M-A-III-24hr_ReadsPerGene.out.tab", "M-B-III-24hr_ReadsPerGene.out.tab", "M-C-III-24hr_ReadsPerGene.out.tab", "M-D-III-24hr_ReadsPerGene.out.tab"))
{
  file<-read.table(f, header = TRUE, sep = "\t")
  infection<-cbind(infection,file[,2])
}


rownames(infection)<-file[,1]
colnames(infection)<-c("M-A-I-3hr", "M-B-I-3hr", "M-C-I-3hr", "M-D-I-3hr",
                       "M-A-II-3hr", "M-B-II-3hr", "M-C-II-3hr", "M-D-II-3hr",
                       "M-A-III-3hr", "M-B-III-3hr", "M-C-III-3hr", "M-D-III-3hr",
                       "M-A-I-24hr", "M-B-I-24hr", "M-C-I-24hr", "M-D-I-24hr",
                       "M-A-II-24hr", "M-B-II-24hr", "M-C-II-24hr", "M-D-II-24hr", 
                       "M-A-III-24hr", "M-B-III-24hr", "M-C-III-24hr", "M-D-III-24hr")
infection<-infection[-c(1,2,3),]

infection_fungus <- infection[rownames(infection) %in% rownames(macrophage_controls),]


### all macrophage data

macrophage_data <- merge(infection_fungus,macrophage_controls, by="row.names")
rownames(macrophage_data)<- macrophage_data$Row.names
macrophage_data <- macrophage_data[,-c(1)]

### PCA plot for all data

Strain <- factor(c(rep(c("A","B","C", "D"),6),
                   rep(c("Blank"),4)), 
                 levels=c("A","B","C", "D", "Blank"))

Time <- factor(c("3hr", "3hr", "3hr", "3hr",
                 "3hr", "3hr", "3hr", "3hr",
                 "3hr", "3hr", "3hr", "3hr",
                 "24hr", "24hr", "24hr", "24hr",
                 "24hr", "24hr", "24hr", "24hr", 
                 "24hr", "24hr", "24hr", "24hr",
                 "3hr","3hr","24hr", "24hr"), 
               levels=c("3hr","24hr"))

Media <- factor(c(rep(c("infection"),24),
                  rep(c("initial"),4)),
                levels = c("infection","initial"))

colData_pca <- DataFrame(Strain=Strain,Time=Time,Media=Media)
colData_pca$Condition <- paste0(colData_pca$Time,"_",colData_pca$Media)


pre_dds_pca <- DESeqDataSetFromMatrix(macrophage_data, colData_pca, design = ~Time)#just a fake desing, does not affect the data normalization and transforamtion
dds_pca <- DESeq(pre_dds_pca)

### PCA
for_pca_vst<-vst(dds_pca)
for_pca_vst_counts<-assay(for_pca_vst)

pca <- prcomp(t(for_pca_vst_counts))
percentVar <- pca$sdev^2/sum(pca$sdev^2)
df_plotting<-data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2],colData_pca) 

ggplot(df_plotting, aes(PC1,PC2,label=rownames(df_plotting)))+
  geom_point(aes(colour=Media), size=5)+
  geom_text_repel(label=rownames(df_plotting), size = 3, segment.color = "black",arrow = arrow(length = unit(0, 'npc')))+
  xlab(paste0("PC1: ",round(percentVar[1]*100),"% variance")) +
  ylab(paste0("PC2: ",round(percentVar[2]*100),"% variance")) + 
  coord_fixed()+ theme_bw()+guides(colour = guide_legend(order = 1),size = guide_legend(order = 2))+
  ggsave("~/TFM/r/results/pcaplots/2M/pca.png", width = 8, height = 6, units = "in", dpi = 300)

################
#### DESeq2 ####
################


macrophage_data_for_deseq2 <- macrophage_data[,-c(25:28)]
group1g <- factor(x=c("A", "B", "C", "D",
                      "A", "B", "C", "D",
                      "A", "B", "C", "D",
                      "A", "B", "C", "D",
                      "A", "B", "C", "D", 
                      "A", "B", "C", "D"),
                  levels=c("A", "B", "C", "D"))
group2g <- factor(x=c("3hr", "3hr", "3hr", "3hr",
                      "3hr", "3hr", "3hr", "3hr",
                      "3hr", "3hr", "3hr", "3hr",
                      "24hr", "24hr", "24hr", "24hr",
                      "24hr", "24hr", "24hr", "24hr", 
                      "24hr", "24hr", "24hr", "24hr"),
                  levels=c("3hr","24hr"))
colDatag <- DataFrame(strain=group1g,time=group2g)
pre_ddsg <- DESeqDataSetFromMatrix(macrophage_data_for_deseq2, colDatag, design = ~time+strain)#variable of interest is the second one
ddsg <- DESeq(pre_ddsg)

B_vs_A <- results(ddsg, contrast = c("strain","B","A"))
C_vs_A <- results(ddsg, contrast = c("strain","C","A"))
D_vs_A <- results(ddsg, contrast = c("strain","D","A"))

C_vs_B <- results(ddsg, contrast = c("strain","C","B"))
D_vs_B <- results(ddsg, contrast = c("strain","D","B"))

D_vs_C <- results(ddsg, contrast = c("strain","D","C"))

###

get_up_down_reg <- function(res, name){
  res <- as.data.frame(res)
  up <- row.names(res[!is.na(res$padj) & res$padj<0.01 & res$log2FoldChange > 1,])
  down <- row.names(res[!is.na(res$padj) & res$padj<0.01 & res$log2FoldChange < -1,])
  
  assign(paste0(name,"_up"),up,envir = .GlobalEnv)
  assign(paste0(name,"_down"),down,envir = .GlobalEnv)
  
}

###

get_up_down_reg(B_vs_A,"B_vs_A")
get_up_down_reg(C_vs_A,"C_vs_A")
get_up_down_reg(D_vs_A,"D_vs_A")

length(B_vs_A_up)
length(B_vs_A_down)

length(C_vs_A_up)
length(C_vs_A_down)

length(D_vs_A_up)
length(D_vs_A_down)

###


B_vs_A_up <- row.names(B_vs_A[!is.na(B_vs_A$padj) & B_vs_A$padj<0.01 & B_vs_A$log2FoldChange > 1,])
B_vs_A_down <- row.names(B_vs_A[!is.na(B_vs_A$padj) & B_vs_A$padj<0.01 & B_vs_A$log2FoldChange < -1,])

C_vs_A_up <- row.names(C_vs_A[!is.na(C_vs_A$padj) & C_vs_A$padj<0.01 & C_vs_A$log2FoldChange > 1,])
C_vs_A_down <- row.names(C_vs_A[!is.na(C_vs_A$padj) & C_vs_A$padj<0.01 & C_vs_A$log2FoldChange < -1,])

D_vs_A_up <- row.names(D_vs_A[!is.na(D_vs_A$padj) & D_vs_A$padj<0.01 & D_vs_A$log2FoldChange > 1,])
D_vs_A_down <- row.names(D_vs_A[!is.na(D_vs_A$padj) & D_vs_A$padj<0.01 & D_vs_A$log2FoldChange < -1,])

#################################################################
####################   END   ####################################
#################################################################

library(clusterProfiler)
library(GO.db)
library(org.Hs.eg.db)

setwd("~/TFM/r")

### in the function, "file" is just a list of dif. expressed gene, "name" is just the name of the output file you want,
### and ontology is Biological Process (BP), Molecular Function (MF), Cellular Compartment (CC) 

############ HUMAN GO TERMS  #####################

human_goterms<-function(file,name,ontology){
  ego<-enrichGO(file, 
                OrgDb = org.Hs.eg.db,
                keyType = "ENSEMBL",
                pvalueCutoff = 0.05, 
                pAdjustMethod = "BH",
                ont = ontology)
  if (!is.null(ego)){
    p<-dotplot(ego,showCategory = 30)
    ggplot2::ggsave(sprintf("%s_%s.png",name,ontology),
                    units="in", width=10, heigh=7, dpi=600)
    write.table(ego,sprintf("%s_%s.txt",name,ontology), sep="\t")
    ego_df <- as.data.frame(ego)
    ego_df <- cbind.data.frame("Cluster"=name,ego_df)
    assign(paste0(substring(name,4),"_df"),as.data.frame(ego_df), envir = .GlobalEnv)
    
  }
}

###

B_vs_A_up <- B_vs_A_up[!(grepl("^LINC",B_vs_A_up) | grepl("^LOC",B_vs_A_up))]
B_vs_A_up <- sub("_1","",B_vs_A_up)
B_vs_A_up <- sub("_2","",B_vs_A_up)
mapIds(org.Hs.eg.db, B_vs_A_up, "REFSEQ", "SYMBOL")
library(EnsDb.Hsapiens.v86)
hsens=EnsDb.Hsapiens.v86
my.symbols <- B_vs_A_up
B_vs_A_up <- select(hsens,  
                keys = my.symbols, 
                columns = c("ENTREZID", "SYMBOL", "GENEID"), 
                keytype = "SYMBOL")

B_vs_A_down <- B_vs_A_down[!(grepl("^LINC",B_vs_A_down) | grepl("^LOC",B_vs_A_down))]
B_vs_A_down <- sub("_1","",B_vs_A_down)
B_vs_A_down <- sub("_2","",B_vs_A_down)
mapIds(org.Hs.eg.db, B_vs_A_down, "REFSEQ", "SYMBOL")
library(EnsDb.Hsapiens.v86)
hsens=EnsDb.Hsapiens.v86
my.symbols <- B_vs_A_down
B_vs_A_down <- select(hsens,  
                    keys = my.symbols, 
                    columns = c("ENTREZID", "SYMBOL", "GENEID"), 
                    keytype = "SYMBOL")

C_vs_A_up <- C_vs_A_up[!(grepl("^LINC",C_vs_A_up) | grepl("^LOC",C_vs_A_up))]
C_vs_A_up <- sub("_1","",C_vs_A_up)
C_vs_A_up <- sub("_2","",C_vs_A_up)
mapIds(org.Hs.eg.db, C_vs_A_up, "REFSEQ", "SYMBOL")
library(EnsDb.Hsapiens.v86)
hsens=EnsDb.Hsapiens.v86
my.symbols <- C_vs_A_up
C_vs_A_up <- select(hsens,  
                    keys = my.symbols, 
                    columns = c("ENTREZID", "SYMBOL", "GENEID"), 
                    keytype = "SYMBOL")

C_vs_A_down <- C_vs_A_down[!(grepl("^LINC",C_vs_A_down) | grepl("^LOC",C_vs_A_down))]
C_vs_A_down <- sub("_1","",C_vs_A_down)
C_vs_A_down <- sub("_2","",C_vs_A_down)
mapIds(org.Hs.eg.db, C_vs_A_down, "REFSEQ", "SYMBOL")
library(EnsDb.Hsapiens.v86)
hsens=EnsDb.Hsapiens.v86
my.symbols <- C_vs_A_down
C_vs_A_down <- select(hsens,  
                    keys = my.symbols, 
                    columns = c("ENTREZID", "SYMBOL", "GENEID"), 
                    keytype = "SYMBOL")

D_vs_A_up <- D_vs_A_up[!(grepl("^LINC",D_vs_A_up) | grepl("^LOC",D_vs_A_up))]
D_vs_A_up <- sub("_1","",D_vs_A_up)
D_vs_A_up <- sub("_2","",D_vs_A_up)
mapIds(org.Hs.eg.db, D_vs_A_up, "REFSEQ", "SYMBOL")
library(EnsDb.Hsapiens.v86)
hsens=EnsDb.Hsapiens.v86
my.symbols <- D_vs_A_up
D_vs_A_up <- select(hsens,  
                    keys = my.symbols, 
                    columns = c("ENTREZID", "SYMBOL", "GENEID"), 
                    keytype = "SYMBOL")

D_vs_A_down <- D_vs_A_down[!(grepl("^LINC",D_vs_A_down) | grepl("^LOC",D_vs_A_down))]
D_vs_A_down <- sub("_1","",D_vs_A_down)
D_vs_A_down <- sub("_2","",D_vs_A_down)
mapIds(org.Hs.eg.db, D_vs_A_down, "REFSEQ", "SYMBOL")
library(EnsDb.Hsapiens.v86)
hsens=EnsDb.Hsapiens.v86
my.symbols <- D_vs_A_down
D_vs_A_down <- select(hsens,  
                      keys = my.symbols, 
                      columns = c("ENTREZID", "SYMBOL", "GENEID"), 
                      keytype = "SYMBOL")

###

setwd("~/TFM/r/results")

human_goterms(B_vs_A_up$GENEID, "2M_B_vs_A_up", "BP")
human_goterms(B_vs_A_up$GENEID, "2M_B_vs_A_up", "MF")
human_goterms(B_vs_A_up$GENEID, "2M_B_vs_A_up", "CC")

human_goterms(C_vs_A_up$GENEID, "2M_C_vs_A_up", "BP")
human_goterms(C_vs_A_up$GENEID, "2M_C_vs_A_up", "MF")
human_goterms(C_vs_A_up$GENEID, "2M_C_vs_A_up", "CC")

human_goterms(D_vs_A_up$GENEID, "2M_D_vs_A_up", "BP")
human_goterms(D_vs_A_up$GENEID, "2M_D_vs_A_up", "MF")
human_goterms(D_vs_A_up$GENEID, "2M_D_vs_A_up", "CC")

###

human_goterms(B_vs_A_down$GENEID, "2M_B_vs_A_down", "BP")
human_goterms(B_vs_A_down$GENEID, "2M_B_vs_A_down", "MF")
human_goterms(B_vs_A_down$GENEID, "2M_B_vs_A_down", "CC")

human_goterms(C_vs_A_down$GENEID, "2M_C_vs_A_down", "BP")
human_goterms(C_vs_A_down$GENEID, "2M_C_vs_A_down", "MF")
human_goterms(C_vs_A_down$GENEID, "2M_C_vs_A_down", "CC")

human_goterms(D_vs_A_down$GENEID, "2M_D_vs_A_down", "BP")
human_goterms(D_vs_A_down$GENEID, "2M_D_vs_A_down", "MF")
human_goterms(D_vs_A_down$GENEID, "2M_D_vs_A_down", "CC")

###

B_vs_A_maplot <- ggmaplot(B_vs_A,
                          fdr = 0.01, fc = 2, size = 2,
                          palette = c("#B31B21", "#1465AC", "darkgray"),
                          genenames = as.vector(row.names(B_vs_A)),
                          legend = "top",
                          font.legend = c("plain", 15),
                          ggtheme = ggplot2::theme_bw(),
                          top = 0) + font("xlab", size = 10) +
  font("ylab", size = 10) +
  font("xy.text", size = 10)
ggsave("~/TFM/r/results/maplots/2M/B_vs_A_macrophage.png", plot = B_vs_A_maplot, width = 8, height = 4, units = "in", dpi = 300)

C_vs_A_maplot <- ggmaplot(C_vs_A,
                          fdr = 0.01, fc = 2, size = 2,
                          palette = c("#B31B21", "#1465AC", "darkgray"),
                          genenames = as.vector(row.names(C_vs_A)),
                          legend = "top",
                          font.legend = c("plain", 15),
                          ggtheme = ggplot2::theme_bw(),
                          top = 0) + font("xlab", size = 10) +
  font("ylab", size = 10) +
  font("xy.text", size = 10)
ggsave("~/TFM/r/results/maplots/2M/C_vs_A_macrophage.png", plot = C_vs_A_maplot, width = 8, height = 4, units = "in", dpi = 300)

D_vs_A_maplot <- ggmaplot(D_vs_A,
                          fdr = 0.01, fc = 2, size = 2,
                          palette = c("#B31B21", "#1465AC", "darkgray"),
                          genenames = as.vector(row.names(D_vs_A)),
                          legend = "top",
                          font.legend = c("plain", 15),
                          ggtheme = ggplot2::theme_bw(),
                          top = 0) + font("xlab", size = 10) +
  font("ylab", size = 10) +
  font("xy.text", size = 10)
ggsave("~/TFM/r/results/maplots/2M/D_vs_A_macrophage.png", plot = D_vs_A_maplot, width = 8, height = 4, units = "in", dpi = 300)

C_vs_B_maplot <- ggmaplot(C_vs_B,
                          fdr = 0.01, fc = 2, size = 2,
                          palette = c("#B31B21", "#1465AC", "darkgray"),
                          genenames = as.vector(row.names(C_vs_B)),
                          legend = "top",
                          font.legend = c("plain", 15),
                          ggtheme = ggplot2::theme_bw(),
                          top = 0) + font("xlab", size = 10) +
  font("ylab", size = 10) +
  font("xy.text", size = 10)
ggsave("~/TFM/r/results/maplots/2M/C_vs_B_macrophage.png", plot = C_vs_B_maplot, width = 8, height = 4, units = "in", dpi = 300)

D_vs_B_maplot <- ggmaplot(D_vs_B,
                          fdr = 0.01, fc = 2, size = 2,
                          palette = c("#B31B21", "#1465AC", "darkgray"),
                          genenames = as.vector(row.names(D_vs_B)),
                          legend = "top",
                          font.legend = c("plain", 15),
                          ggtheme = ggplot2::theme_bw(),
                          top = 0) + font("xlab", size = 10) +
  font("ylab", size = 10) +
  font("xy.text", size = 10)
ggsave("~/TFM/r/results/maplots/2M/D_vs_B_macrophage.png", plot = D_vs_B_maplot, width = 8, height = 4, units = "in", dpi = 300)

D_vs_C_maplot <- ggmaplot(D_vs_C,
                          fdr = 0.01, fc = 2, size = 2,
                          palette = c("#B31B21", "#1465AC", "darkgray"),
                          genenames = as.vector(row.names(D_vs_C)),
                          legend = "top",
                          font.legend = c("plain", 15),
                          ggtheme = ggplot2::theme_bw(),
                          top = 0) + font("xlab", size = 10) +
  font("ylab", size = 10) +
  font("xy.text", size = 10)
ggsave("~/TFM/r/results/maplots/2M/D_vs_C_macrophage.png", plot = D_vs_C_maplot, width = 8, height = 4, units = "in", dpi = 300)

###

all_up<-rbind.data.frame(B_vs_A_up_df,
                         C_vs_A_up_df,
                         D_vs_A_up_df)
rownames(all_up)<-NULL

all_up$Cluster <- ifelse(grepl("B", all_up$Cluster), "B vs A", 
                         ifelse(grepl("C",all_up$Cluster),"C vs A","D vs A"))


all_down <-rbind.data.frame(B_vs_A_down_df,
                            C_vs_A_down_df,
                            D_vs_A_down_df)
rownames(all_down)<-NULL

all_down$Cluster <- ifelse(grepl("B", all_down$Cluster), "B vs A", 
                           ifelse(grepl("C",all_down$Cluster),"C vs A","D vs A"))


## clusterprofiler exampledata
mydf <- data.frame(Entrez=c('1', '100', '1000', '100101467',
                            '100127206', '100128071'),
                   logFC = c(1.1, -0.5, 5, 2.5, -3, 3),
                   group = c('A', 'A', 'A', 'B', 'B', 'B'),
                   othergroup = c('good', 'good', 'bad', 'bad', 'good', 'bad'))

to_hack <- compareCluster(Entrez~group+othergroup, data=mydf,
                          fun='enrichGO', OrgDb='org.Hs.eg.db')

to_hack_up <- to_hack
to_hack_down <- to_hack
### hack to use for dotplot



to_hack_up@compareClusterResult <- all_up

### UP
plot_up<-dotplot(to_hack_up, showCategory=15)
plot_up+ggplot2::theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))+
  # ggsave("./results/project_4/up_GO_general.pdf",units="in", width=14, height=7, dpi=300)
  ggsave("~/TFM/r/results/GO/2M/up_GO_general.png",units="in", width=14, height=7, dpi=300)



## down

to_hack_down@compareClusterResult <- all_down
plot_down<-dotplot(to_hack_down, showCategory=15)
plot_down+ggplot2::theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))+
  # ggsave("./results/project_4/down_GO_general.pdf",units="in", width=10, height=7, dpi=300)
  ggsave("~/TFM/r/results/GO/2M/down_GO_general.png",units="in", width=10, height=7, dpi=300)




####


### Creating an excel file with DE data

wb <- createWorkbook("Project_2_excel")

header<-c("gene", "baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")
file_names<-c("B_vs_A",
              "C_vs_A",
              "D_vs_A")

descriptions<-c("fluconazole-resistant (FLZR) isolate derived from CBS138 vs CBS138",
                "MDR isolate derived from the same FLZR isolate (from isolate B), which carries R631G in HS1-FKS1 vs CBS138",
                "MDR derived from the same FLZR background (from isolate B), which carries S663P in HS1-FKS2 vs CBS138")

sheets_and_descriptions<-cbind.data.frame("Data_sheets"=file_names,"Description"=descriptions)


addWorksheet(wb,"Sheets_and_descriptions")
writeData(wb, "Sheets_and_descriptions", sheets_and_descriptions,keepNA = T)

for (file in file_names){
  assign("fung",as.data.frame(eval(parse(text=paste0(file)))))
  fung <- cbind("genes"=rownames(fung), data.frame(fung, row.names=NULL))
  addWorksheet(wb, file)
  writeData(wb, file, fung,keepNA = T)
}



#saveWorkbook(wb, "./results/project_4/Project_4_excel.xlsx", overwrite = TRUE)
saveWorkbook(wb, "~/TFM/r/results/excel/2M/Project_2_excel.xlsx", overwrite = TRUE)

