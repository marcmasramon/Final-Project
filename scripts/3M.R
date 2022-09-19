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
for (f in c("M-A-I-3hr_ReadsPerGene.out.tab", "M-E-I-3hr_ReadsPerGene.out.tab", "M-F-I-3hr_ReadsPerGene.out.tab", "M-G-I-3hr_ReadsPerGene.out.tab",
            "M-A-II-3hr_ReadsPerGene.out.tab", "M-E-II-3hr_ReadsPerGene.out.tab", "M-F-II-3hr_ReadsPerGene.out.tab", "M-G-II-3hr_ReadsPerGene.out.tab",
            "M-A-III-3hr_ReadsPerGene.out.tab", "M-E-III-3hr_ReadsPerGene.out.tab", "M-F-III-3hr_ReadsPerGene.out.tab", "M-G-III-3hr_ReadsPerGene.out.tab",
            "M-A-I-24hr_ReadsPerGene.out.tab", "M-E-I-24hr_ReadsPerGene.out.tab", "M-F-I-24hr_ReadsPerGene.out.tab", "M-G-I-24hr_ReadsPerGene.out.tab",
            "M-A-II-24hr_ReadsPerGene.out.tab", "M-E-II-24hr_ReadsPerGene.out.tab", "M-F-II-24hr_ReadsPerGene.out.tab", "M-G-II-24hr_ReadsPerGene.out.tab", 
            "M-A-III-24hr_ReadsPerGene.out.tab", "M-E-III-24hr_ReadsPerGene.out.tab", "M-F-III-24hr_ReadsPerGene.out.tab", "M-G-III-24hr_ReadsPerGene.out.tab"))
{
  file<-read.table(f, header = TRUE, sep = "\t")
  infection<-cbind(infection,file[,2])
}


rownames(infection)<-file[,1]
colnames(infection)<-c("M-A-I-3hr", "M-E-I-3hr", "M-F-I-3hr", "M-G-I-3hr",
                       "M-A-II-3hr", "M-E-II-3hr", "M-F-II-3hr", "M-G-II-3hr",
                       "M-A-III-3hr", "M-E-III-3hr", "M-F-III-3hr", "M-G-III-3hr",
                       "M-A-I-24hr", "M-E-I-24hr", "M-F-I-24hr", "M-G-I-24hr",
                       "M-A-II-24hr", "M-E-II-24hr", "M-F-II-24hr", "M-G-II-24hr", 
                       "M-A-III-24hr", "M-E-III-24hr", "M-F-III-24hr", "M-G-III-24hr")
infection<-infection[-c(1,2,3),]

infection_fungus <- infection[rownames(infection) %in% rownames(macrophage_controls),]


### all macrophage data

macrophage_data <- merge(infection_fungus,macrophage_controls, by="row.names")
rownames(macrophage_data)<- macrophage_data$Row.names
macrophage_data <- macrophage_data[,-c(1)]

### PCA plot for all data

Strain <- factor(c(rep(c("A","E","F", "G"),6),
                   rep(c("Blank"),4)), 
                 levels=c("A","E","F", "G", "Blank"))

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
  ggsave("~/TFM/r/results/pcaplots/3M/pca.png", width = 8, height = 6, units = "in", dpi = 300)


################
#### DESeq2 ####
################


macrophage_data_for_deseq2 <- macrophage_data[,-c(25:28)]
group1g <- factor(x=c("A", "E", "F", "G",
                      "A", "E", "F", "G",
                      "A", "E", "F", "G",
                      "A", "E", "F", "G",
                      "A", "E", "F", "G", 
                      "A", "E", "F", "G"),
                  levels=c("A", "E", "F", "G"))
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

E_vs_A <- results(ddsg, contrast = c("strain","E","A"))
F_vs_A <- results(ddsg, contrast = c("strain","F","A"))
G_vs_A <- results(ddsg, contrast = c("strain","G","A"))

F_vs_E <- results(ddsg, contrast = c("strain","F","E"))
G_vs_E <- results(ddsg, contrast = c("strain","G","E"))

G_vs_F <- results(ddsg, contrast = c("strain","G","F"))

###

get_up_down_reg <- function(res, name){
  res <- as.data.frame(res)
  up <- row.names(res[!is.na(res$padj) & res$padj<0.01 & res$log2FoldChange > 1,])
  down <- row.names(res[!is.na(res$padj) & res$padj<0.01 & res$log2FoldChange < -1,])
  
  assign(paste0(name,"_up"),up,envir = .GlobalEnv)
  assign(paste0(name,"_down"),down,envir = .GlobalEnv)
  
}

###

get_up_down_reg(E_vs_A,"E_vs_A")
get_up_down_reg(F_vs_A,"F_vs_A")
get_up_down_reg(F_vs_A,"G_vs_A")
get_up_down_reg(G_vs_E,"G_vs_E")
get_up_down_reg(G_vs_F,"G_vs_F")

length(E_vs_A_up)
length(E_vs_A_down)

length(F_vs_A_up)
length(F_vs_A_down)

length(G_vs_A_up)
length(G_vs_A_down)

length(G_vs_E_up)
length(G_vs_E_down)

length(G_vs_F_up)
length(G_vs_F_down)

###


E_vs_A_up <- row.names(E_vs_A[!is.na(E_vs_A$padj) & E_vs_A$padj<0.01 & E_vs_A$log2FoldChange > 1,])
E_vs_A_down <- row.names(E_vs_A[!is.na(E_vs_A$padj) & E_vs_A$padj<0.01 & E_vs_A$log2FoldChange < -1,])

F_vs_A_up <- row.names(F_vs_A[!is.na(F_vs_A$padj) & F_vs_A$padj<0.01 & F_vs_A$log2FoldChange > 1,])
F_vs_A_down <- row.names(F_vs_A[!is.na(F_vs_A$padj) & F_vs_A$padj<0.01 & F_vs_A$log2FoldChange < -1,])

G_vs_A_up <- row.names(G_vs_A[!is.na(G_vs_A$padj) & G_vs_A$padj<0.01 & G_vs_A$log2FoldChange > 1,])
G_vs_A_down <- row.names(G_vs_A[!is.na(G_vs_A$padj) & G_vs_A$padj<0.01 & G_vs_A$log2FoldChange < -1,])


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
    ggplot2::ggsave(sprintf("%s_%s_macrophage.png",name,ontology),
                    units="in", width=10, heigh=7, dpi=600)
    write.table(ego,sprintf("%s_%s.txt",name,ontology), sep="\t")
    ego_df <- as.data.frame(ego)
    ego_df <- cbind.data.frame("Cluster"=name,ego_df)
    assign(paste0(substring(name,4),"_df"),as.data.frame(ego_df), envir = .GlobalEnv)
    
  }
}

###

E_vs_A_up <- E_vs_A_up[!(grepl("^LINC",E_vs_A_up) | grepl("^LOC",E_vs_A_up))]
E_vs_A_up <- sub("_1","",E_vs_A_up)
E_vs_A_up <- sub("_2","",E_vs_A_up)
mapIds(org.Hs.eg.db, E_vs_A_up, "REFSEQ", "SYMBOL")
library(EnsDb.Hsapiens.v86)
hsens=EnsDb.Hsapiens.v86
my.symbols <- E_vs_A_up
E_vs_A_up <- select(hsens,  
                    keys = my.symbols, 
                    columns = c("ENTREZID", "SYMBOL", "GENEID"), 
                    keytype = "SYMBOL")

E_vs_A_down <- E_vs_A_down[!(grepl("^LINC",E_vs_A_down) | grepl("^LOC",E_vs_A_down))]
E_vs_A_down <- sub("_1","",E_vs_A_down)
E_vs_A_down <- sub("_2","",E_vs_A_down)
mapIds(org.Hs.eg.db, E_vs_A_down, "REFSEQ", "SYMBOL")
library(EnsDb.Hsapiens.v86)
hsens=EnsDb.Hsapiens.v86
my.symbols <- E_vs_A_down
E_vs_A_down <- select(hsens,  
                      keys = my.symbols, 
                      columns = c("ENTREZID", "SYMBOL", "GENEID"), 
                      keytype = "SYMBOL")

F_vs_A_up <- F_vs_A_up[!(grepl("^LINC",F_vs_A_up) | grepl("^LOC",F_vs_A_up))]
F_vs_A_up <- sub("_1","",F_vs_A_up)
F_vs_A_up <- sub("_2","",F_vs_A_up)
mapIds(org.Hs.eg.db, F_vs_A_up, "REFSEQ", "SYMBOL")
library(EnsDb.Hsapiens.v86)
hsens=EnsDb.Hsapiens.v86
my.symbols <- F_vs_A_up
F_vs_A_up <- select(hsens,  
                    keys = my.symbols, 
                    columns = c("ENTREZID", "SYMBOL", "GENEID"), 
                    keytype = "SYMBOL")

F_vs_A_down <- F_vs_A_down[!(grepl("^LINC",F_vs_A_down) | grepl("^LOC",F_vs_A_down))]
F_vs_A_down <- sub("_1","",F_vs_A_down)
F_vs_A_down <- sub("_2","",F_vs_A_down)
mapIds(org.Hs.eg.db, F_vs_A_down, "REFSEQ", "SYMBOL")
library(EnsDb.Hsapiens.v86)
hsens=EnsDb.Hsapiens.v86
my.symbols <- F_vs_A_down
F_vs_A_down <- select(hsens,  
                      keys = my.symbols, 
                      columns = c("ENTREZID", "SYMBOL", "GENEID"), 
                      keytype = "SYMBOL")

G_vs_A_up <- G_vs_A_up[!(grepl("^LINC",G_vs_A_up) | grepl("^LOC",G_vs_A_up))]
G_vs_A_up <- sub("_1","",G_vs_A_up)
G_vs_A_up <- sub("_2","",G_vs_A_up)
mapIds(org.Hs.eg.db, G_vs_A_up, "REFSEQ", "SYMBOL")
library(EnsDb.Hsapiens.v86)
hsens=EnsDb.Hsapiens.v86
my.symbols <- G_vs_A_up
G_vs_A_up <- select(hsens,  
                    keys = my.symbols, 
                    columns = c("ENTREZID", "SYMBOL", "GENEID"), 
                    keytype = "SYMBOL")

G_vs_A_down <- G_vs_A_down[!(grepl("^LINC",G_vs_A_down) | grepl("^LOC",G_vs_A_down))]
G_vs_A_down <- sub("_1","",G_vs_A_down)
G_vs_A_down <- sub("_2","",G_vs_A_down)
mapIds(org.Hs.eg.db, G_vs_A_down, "REFSEQ", "SYMBOL")
library(EnsDb.Hsapiens.v86)
hsens=EnsDb.Hsapiens.v86
my.symbols <- G_vs_A_down
G_vs_A_down <- select(hsens,  
                      keys = my.symbols, 
                      columns = c("ENTREZID", "SYMBOL", "GENEID"), 
                      keytype = "SYMBOL")

###

setwd("~/TFM/r/results")

human_goterms(E_vs_A_up$GENEID, "3M_E_vs_A_up", "BP")
human_goterms(E_vs_A_up$GENEID, "3M_E_vs_A_up", "MF")
human_goterms(E_vs_A_up$GENEID, "3M_E_vs_A_up", "CC")

human_goterms(F_vs_A_up$GENEID, "3M_F_vs_A_up", "BP")
human_goterms(F_vs_A_up$GENEID, "3M_F_vs_A_up", "MF")
human_goterms(F_vs_A_up$GENEID, "3M_F_vs_A_up", "CC")

human_goterms(G_vs_A_up$GENEID, "3M_G_vs_A_up", "BP")
human_goterms(G_vs_A_up$GENEID, "3M_G_vs_A_up", "MF")
human_goterms(G_vs_A_up$GENEID, "3M_G_vs_A_up", "CC")

###

human_goterms(E_vs_A_down$GENEID, "3M_E_vs_A_down", "BP")
human_goterms(E_vs_A_down$GENEID, "3M_E_vs_A_down", "MF")
human_goterms(E_vs_A_down$GENEID, "3M_E_vs_A_down", "CC")

human_goterms(F_vs_A_down$GENEID, "3M_F_vs_A_down", "BP")
human_goterms(F_vs_A_down$GENEID, "3M_F_vs_A_down", "MF")
human_goterms(F_vs_A_down$GENEID, "3M_F_vs_A_down", "CC")

human_goterms(G_vs_A_down$GENEID, "3M_G_vs_A_down", "BP")
human_goterms(G_vs_A_down$GENEID, "3M_G_vs_A_down", "MF")
human_goterms(G_vs_A_down$GENEID, "3M_G_vs_A_down", "CC")

###

E_vs_A_maplot <- ggmaplot(E_vs_A,
                             fdr = 0.01, fc = 2, size = 2,
                             palette = c("#B31B21", "#1465AC", "darkgray"),
                             genenames = as.vector(row.names(E_vs_A)),
                             legend = "top",
                             font.legend = c("plain", 15),
                             ggtheme = ggplot2::theme_bw(),
                             top = 0) + font("xlab", size = 10) +
  font("ylab", size = 10) +
  font("xy.text", size = 10)
ggsave("~/TFM/r/results/maplots/3M/E_vs_A_macrophage.png", plot = E_vs_A_maplot, width = 8, height = 4, units = "in", dpi = 300)

F_vs_A_maplot <- ggmaplot(F_vs_A,
                           fdr = 0.01, fc = 2, size = 2,
                           palette = c("#B31B21", "#1465AC", "darkgray"),
                           genenames = as.vector(row.names(F_vs_A)),
                           legend = "top",
                           font.legend = c("plain", 15),
                           ggtheme = ggplot2::theme_bw(),
                           top = 0) + font("xlab", size = 10) +
  font("ylab", size = 10) +
  font("xy.text", size = 10)
ggsave("~/TFM/r/results/maplots/3M/F_vs_A_macrophage.png", plot = F_vs_A_maplot, width = 8, height = 4, units = "in", dpi = 300)

G_vs_A_maplot <- ggmaplot(G_vs_A,
                              fdr = 0.01, fc = 2, size = 2,
                              palette = c("#B31B21", "#1465AC", "darkgray"),
                              genenames = as.vector(row.names(G_vs_A)),
                              legend = "top",
                              font.legend = c("plain", 15),
                              ggtheme = ggplot2::theme_bw(),
                              top = 0) + font("xlab", size = 10) +
  font("ylab", size = 10) +
  font("xy.text", size = 10)
ggsave("~/TFM/r/results/maplots/3M/G_vs_A_macrophage.png", plot = G_vs_A_maplot, width = 8, height = 4, units = "in", dpi = 300)

F_vs_E_maplot <- ggmaplot(F_vs_E,
                                   fdr = 0.01, fc = 2, size = 2,
                                   palette = c("#B31B21", "#1465AC", "darkgray"),
                                   genenames = as.vector(row.names(F_vs_E)),
                                   legend = "top",
                                   font.legend = c("plain", 15),
                                   ggtheme = ggplot2::theme_bw(),
                                   top = 0) + font("xlab", size = 10) +
  font("ylab", size = 10) +
  font("xy.text", size = 10)
ggsave("~/TFM/r/results/maplots/3M/F_vs_E_macrophage.png", plot = F_vs_E_maplot, width = 8, height = 4, units = "in", dpi = 300)

G_vs_E_maplot <- ggmaplot(G_vs_E,
                              fdr = 0.01, fc = 2, size = 2,
                              palette = c("#B31B21", "#1465AC", "darkgray"),
                              genenames = as.vector(row.names(G_vs_E)),
                              legend = "top",
                              font.legend = c("plain", 15),
                              ggtheme = ggplot2::theme_bw(),
                              top = 0) + font("xlab", size = 10) +
  font("ylab", size = 10) +
  font("xy.text", size = 10)
ggsave("~/TFM/r/results/maplots/3M/G_vs_E_macrophage.png", plot = G_vs_E_maplot, width = 8, height = 4, units = "in", dpi = 300)

G_vs_F_maplot <- ggmaplot(G_vs_F,
                                 fdr = 0.01, fc = 2, size = 2,
                                 palette = c("#B31B21", "#1465AC", "darkgray"),
                                 genenames = as.vector(row.names(G_vs_F)),
                                 legend = "top",
                                 font.legend = c("plain", 15),
                                 ggtheme = ggplot2::theme_bw(),
                                 top = 0) + font("xlab", size = 10) +
  font("ylab", size = 10) +
  font("xy.text", size = 10)
ggsave("~/TFM/r/results/maplots/3M/G_vs_F_macrophage.png", plot = G_vs_F_maplot, width = 8, height = 4, units = "in", dpi = 300)

###

all_up<-rbind.data.frame(E_vs_A_up_df,
                         F_vs_A_up_df,
                         G_vs_A_up_df)
rownames(all_up)<-NULL

all_up$Cluster <- ifelse(grepl("E", all_up$Cluster), "E vs A",
                         ifelse(grepl("F",all_up$Cluster),"F vs A","G vs A"))


all_down <-rbind.data.frame(E_vs_A_down_df,
                            F_vs_A_down_df,
                            G_vs_A_down_df)
rownames(all_down)<-NULL

all_down$Cluster <- ifelse(grepl("E", all_down$Cluster), "E vs A", 
                           ifelse(grepl("F",all_down$Cluster),"F vs A","G vs A"))


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
  ggsave("~/TFM/r/results/GO/3M/up_GO_general.png",units="in", width=14, height=7, dpi=300)



## down

to_hack_down@compareClusterResult <- all_down
plot_down<-dotplot(to_hack_down, showCategory=15)
plot_down+ggplot2::theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))+
  # ggsave("./results/project_4/down_GO_general.pdf",units="in", width=10, height=7, dpi=300)
  ggsave("~/TFM/r/results/GO/3M/down_GO_general.png",units="in", width=10, height=7, dpi=300)




####


### Creating an excel file with DE data

wb <- createWorkbook("Project_3_excel")

header<-c("gene", "baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")
file_names<-c("E_vs_A",
              "F_vs_A",
              "G_vs_A")


descriptions<-c("petite isolate derived from CBS138 vs CBS138",
                "clinical isolate that is not petite vs CBS138",
                "petite obtained from a patient vs CBS138")

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
saveWorkbook(wb, "~/TFM/r/results/excel/3M/Project_3_excel.xlsx", overwrite = TRUE)

