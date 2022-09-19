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

#### FUNGAL DATA

### controls
fungal_controls<-NULL
for (f in c("H2O2-1_ReadsPerGene.out.tab","H2O2-2_ReadsPerGene.out.tab", "H2O2-3_ReadsPerGene.out.tab", 
            "PH1_ReadsPerGene.out.tab", "PH2_ReadsPerGene.out.tab", "PH3_ReadsPerGene.out.tab", 
            "Spent-1_ReadsPerGene.out.tab", "Spent-2_ReadsPerGene.out.tab", "Spent-3-new_ReadsPerGene.out.tab",
            "Spent-H2O2-1-New_ReadsPerGene.out.tab", "Spent-H21O2-2_ReadsPerGene.out.tab", "Spent-H21O2-3_ReadsPerGene.out.tab",
            "Initial-A_ReadsPerGene.out.tab", "R-A-3hr_ReadsPerGene.out.tab"))
{
  file<-read.table(f, header = TRUE, sep = "\t")
  fungal_controls<-cbind(fungal_controls,file[,2])
}

rownames(fungal_controls)<-file[,1]
colnames(fungal_controls)<-c("H2O2-1","H2O2-2", "H2O2-3", 
                             "PH1", "PH2", "PH3", 
                             "Spent-1", "Spent-2", "Spent-3-new",
                             "Spent_H2O2-1-New", "Spent-H21O2-2", "Spent-H21O2-3",
                             "Initial-A", "R-A-3hr")
fungal_controls<-fungal_controls[-c(1,2,3),]


### infection

infection<-NULL
for (f in c("g-A-I-3hr_ReadsPerGene.out.tab",
            "g-A-II-3hr_ReadsPerGene.out.tab",
            "g-A-III-3hr_ReadsPerGene.out.tab"))
{
  file<-read.table(f, header = TRUE, sep = "\t")
  infection<-cbind(infection,file[,2])
}


rownames(infection)<-file[,1]
colnames(infection)<-c("g-A-I-3hr",
                       "g-A-II-3hr",
                       "g-A-III-3hr")
infection<-infection[-c(1,2,3),]

infection_fungus <- infection[rownames(infection) %in% rownames(fungal_controls),]


### all fungal data

fungal_data <- merge(infection_fungus,fungal_controls, by="row.names")
rownames(fungal_data)<- fungal_data$Row.names
fungal_data <- fungal_data[,-c(1)]

### PCA plot for all data

Media<- factor(c(rep(c("Macrophage"),3),
                 rep(c("H2O2"),3),
                 rep(c("PH"),3),
                 rep(c("Spent"),3),
                 rep(c("Spent_H2O2"),3),
                 rep(c("Initial"),1),
                 rep(c("RPMI_3h"),1)),
               levels=c("Macrophage","H2O2","PH", "Spent", "Spent_H2O2","Initial","RPMI_3h"))

colData_pca <- DataFrame(Media=Media)

pre_dds_pca <- DESeqDataSetFromMatrix(fungal_data, colData_pca, design = ~Media)#just a fake desing, does not affect the data normalization and transforamtion
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
  ggsave("~/TFM/r/results/pcaplots/4F/pca.png", width = 8, height = 6, units = "in", dpi = 300)

################
#### DESeq2 ####
################


fungal_data_for_deseq2 <- fungal_data[,-c(16:17)]
Media<- factor(c(rep(c("Macrophage"),3),
                 rep(c("H2O2"),3),
                 rep(c("PH"),3),
                 rep(c("Spent"),3),
                 rep(c("Spent_H2O2"),3)),
               levels=c("Macrophage","H2O2","PH", "Spent", "Spent_H2O2"))
colDatag <- DataFrame(Media=Media)
pre_ddsg <- DESeqDataSetFromMatrix(fungal_data_for_deseq2, colDatag, design = ~Media)#variable of interest is the second one
ddsg <- DESeq(pre_ddsg)

###

H2O2_vs_Macrophage <- results(ddsg, contrast = c("Media","H2O2","Macrophage"))
PH_vs_Macrophage <- results(ddsg, contrast = c("Media","PH","Macrophage"))
Spent_vs_Macrophage <- results(ddsg, contrast = c("Media","Spent","Macrophage"))
Spent_H2O2_vs_Macrophage <- results(ddsg, contrast = c("Media","Spent_H2O2","Macrophage"))

PH_vs_H2O2 <- results(ddsg, contrast = c("Media","PH","H2O2"))
Spent_vs_H2O2 <- results(ddsg, contrast = c("Media","Spent","H2O2"))
Spent_H2O2_vs_H2O2 <- results(ddsg, contrast = c("Media","Spent_H2O2","H2O2"))

Spent_vs_PH <- results(ddsg, contrast = c("Media","Spent","PH"))
Spent_H2O2_vs_PH <- results(ddsg, contrast = c("Media","Spent_H2O2","PH"))

Spent_H2O2_vs_Spent <- results(ddsg, contrast = c("Media","Spent_H2O2","Spent"))

###

get_up_down_reg <- function(res, name){
  res <- as.data.frame(res)
  up <- row.names(res[!is.na(res$padj) & res$padj<0.01 & res$log2FoldChange > 1,])
  down <- row.names(res[!is.na(res$padj) & res$padj<0.01 & res$log2FoldChange < -1,])
  
  assign(paste0(name,"_up"),up,envir = .GlobalEnv)
  assign(paste0(name,"_down"),down,envir = .GlobalEnv)
  
}

###

get_up_down_reg(H2O2_vs_Macrophage,"H2O2_vs_Macrophage")
get_up_down_reg(PH_vs_Macrophage,"PH_vs_Macrophage")
get_up_down_reg(Spent_vs_Macrophage,"Spent_vs_Macrophage")
get_up_down_reg(Spent_H2O2_vs_Macrophage,"Spent_H2O2_vs_Macrophage")

length(H2O2_vs_Macrophage_up)
length(H2O2_vs_Macrophage_down)

length(PH_vs_Macrophage_up)
length(PH_vs_Macrophage_down)

length(Spent_H2O2_vs_Macrophage_up)
length(Spent_H2O2_vs_Macrophage_down)

length(Spent_vs_Macrophage_up)
length(Spent_vs_Macrophage_down)

###

H2O2_vs_Macrophage_up <- row.names(H2O2_vs_Macrophage[!is.na(H2O2_vs_Macrophage$padj) & H2O2_vs_Macrophage$padj<0.01 & H2O2_vs_Macrophage$log2FoldChange > 1,])
H2O2_vs_Macrophage_down <- row.names(H2O2_vs_Macrophage[!is.na(H2O2_vs_Macrophage$padj) & H2O2_vs_Macrophage$padj<0.01 & H2O2_vs_Macrophage$log2FoldChange < -1,])

PH_vs_Macrophage_up <- row.names(PH_vs_Macrophage[!is.na(PH_vs_Macrophage$padj) & PH_vs_Macrophage$padj<0.01 & PH_vs_Macrophage$log2FoldChange > 1,])
PH_vs_Macrophage_down <- row.names(PH_vs_Macrophage[!is.na(PH_vs_Macrophage$padj) & PH_vs_Macrophage$padj<0.01 & PH_vs_Macrophage$log2FoldChange < -1,])

Spent_vs_Macrophage_up <- row.names(Spent_vs_Macrophage[!is.na(Spent_vs_Macrophage$padj) & Spent_vs_Macrophage$padj<0.01 & Spent_vs_Macrophage$log2FoldChange > 1,])
Spent_vs_Macrophage_down <- row.names(Spent_vs_Macrophage[!is.na(Spent_vs_Macrophage$padj) & Spent_vs_Macrophage$padj<0.01 & Spent_vs_Macrophage$log2FoldChange < -1,])

Spent_H2O2_vs_Macrophage_up <- row.names(Spent_H2O2_vs_Macrophage[!is.na(Spent_H2O2_vs_Macrophage$padj) & Spent_H2O2_vs_Macrophage$padj<0.01 & Spent_H2O2_vs_Macrophage$log2FoldChange > 1,])
Spent_H2O2_vs_Macrophage_down <- row.names(Spent_H2O2_vs_Macrophage[!is.na(Spent_H2O2_vs_Macrophage$padj) & Spent_H2O2_vs_Macrophage$padj<0.01 & Spent_H2O2_vs_Macrophage$log2FoldChange < -1,])

#################################################################
####################   END   ####################################
#################################################################

library(clusterProfiler)
library(GO.db)
library(org.Hs.eg.db)

setwd("~/TFM/r")

### in the function, "file" is just a list of dif. expressed gene, "name" is just the name of the output file you want,
### and ontology is Biological Process (BP), Molecular Function (MF), Cellular Compartment (CC) 


### FUNGAL GO TERMS
### Finding MF,BP and CC terms
all_go_tems<-as.data.frame(GOTERM)

all_go_tems<-all_go_tems[!(all_go_tems$go_id %in% c("GO:0003674","GO:0008150","GO:0005575")),]


MF_terms<-all_go_tems$go_id[all_go_tems$Ontology=="MF"]
BP_terms<-all_go_tems$go_id[all_go_tems$Ontology=="BP"]
CC_terms<-all_go_tems$go_id[all_go_tems$Ontology=="CC"]


### Association of GOID with description
goterms <- Term(GOTERM)
a<-as.data.frame(goterms)
go_names<-cbind(row.names(a),a)


### CGLAB GOIDS
cglab_go<-read.table("cglab_go.txt")


MF_universe<-cglab_go[cglab_go$V1 %in% MF_terms,]
BP_universe<-cglab_go[cglab_go$V1 %in% BP_terms,]
CC_universe<-cglab_go[cglab_go$V1 %in% CC_terms,]



fungal_goterms<-function(file,name,ont){
  ego<-enricher(file, pvalueCutoff = 0.05, 
                pAdjustMethod = "BH", universe = eval(parse(text=paste0("as.character(",ont,"_universe$V2)"))),minGSSize = 2, 
                maxGSSize = 10000, TERM2GENE = eval(parse(text=paste0(ont,"_universe"))),TERM2NAME = go_names)
  
  
  if (!is.null(ego)){
    p<-dotplot(ego,showCategory = 30)
    ggplot2::ggsave(sprintf("./results/4F/%s_%s.png",name,ont),
                    units="in", width=10, heigh=7, dpi=600)
    write.table(ego,sprintf("./results/4F/%s_%s.txt",name,ont), sep="\t")
    ego_df <- as.data.frame(ego)
    ego_df <- cbind.data.frame("Cluster"=name,ego_df)
    assign(paste0(substring(name,4),"_df"),as.data.frame(ego_df), envir = .GlobalEnv)
    
  }
  
}

###

setwd("~/TFM/r/results")

fungal_goterms(H2O2_vs_Macrophage_up, "4F_H2O2_vs_Macrophage_up", "BP")
fungal_goterms(H2O2_vs_Macrophage_up, "4F_H2O2_vs_Macrophage_up", "MF")
fungal_goterms(H2O2_vs_Macrophage_up, "4F_H2O2_vs_Macrophage_up", "CC")

fungal_goterms(PH_vs_Macrophage_up, "4F_PH_vs_Macrophage_up", "BP")
fungal_goterms(PH_vs_Macrophage_up, "4F_PH_vs_Macrophage_up", "MF")
fungal_goterms(PH_vs_Macrophage_up, "4F_PH_vs_Macrophage_up", "CC")

fungal_goterms(Spent_vs_Macrophage_up, "4F_Spent_vs_Macrophage_up", "BP")
fungal_goterms(Spent_vs_Macrophage_up, "4F_Spent_vs_Macrophage_up", "MF")
fungal_goterms(Spent_vs_Macrophage_up, "4F_Spent_vs_Macrophage_up", "CC")

fungal_goterms(Spent_H2O2_vs_Macrophage_up, "4F_Spent_H2O2_vs_Macrophage_up", "BP")
fungal_goterms(Spent_H2O2_vs_Macrophage_up, "4F_Spent_H2O2_vs_Macrophage_up", "MF")
fungal_goterms(Spent_H2O2_vs_Macrophage_up, "4F_Spent_H2O2_vs_Macrophage_up", "CC")

###

fungal_goterms(H2O2_vs_Macrophage_down, "4F_H2O2_vs_Macrophage_down", "BP")
fungal_goterms(H2O2_vs_Macrophage_down, "4F_H2O2_vs_Macrophage_down", "MF")
fungal_goterms(H2O2_vs_Macrophage_down, "4F_H2O2_vs_Macrophage_down", "CC")

fungal_goterms(PH_vs_Macrophage_down, "4F_PH_vs_Macrophage_down", "BP")
fungal_goterms(PH_vs_Macrophage_down, "4F_PH_vs_Macrophage_down", "MF")
fungal_goterms(PH_vs_Macrophage_down, "4F_PH_vs_Macrophage_down", "CC")

fungal_goterms(Spent_vs_Macrophage_down, "4F_Spent_vs_Macrophage_down", "BP")
fungal_goterms(Spent_vs_Macrophage_down, "4F_Spent_vs_Macrophage_down", "MF")
fungal_goterms(Spent_vs_Macrophage_down, "4F_Spent_vs_Macrophage_down", "CC")

fungal_goterms(Spent_H2O2_vs_Macrophage_down, "4F_Spent_H2O2_vs_Macrophage_down", "BP")
fungal_goterms(Spent_H2O2_vs_Macrophage_down, "4F_Spent_H2O2_vs_Macrophage_down", "MF")
fungal_goterms(Spent_H2O2_vs_Macrophage_down, "4F_Spent_H2O2_vs_Macrophage_down", "CC")

###

H2O2_vs_Macrophage_maplot <- ggmaplot(H2O2_vs_Macrophage,
              fdr = 0.01, fc = 2, size = 2,
              palette = c("#B31B21", "#1465AC", "darkgray"),
              genenames = as.vector(row.names(H2O2_vs_Macrophage)),
              legend = "top",
              font.legend = c("plain", 15),
              ggtheme = ggplot2::theme_bw(),
              top = 0) + font("xlab", size = 10) +
              font("ylab", size = 10) +
              font("xy.text", size = 10)
ggsave("~/TFM/r/results/maplots/4F/H2O2_vs_Macrophage_fungal.png", plot = H2O2_vs_Macrophage_maplot, width = 8, height = 4, units = "in", dpi = 300)
  
PH_vs_Macrophage_maplot <- ggmaplot(PH_vs_Macrophage,
                        fdr = 0.01, fc = 2, size = 2,
                        palette = c("#B31B21", "#1465AC", "darkgray"),
                        genenames = as.vector(row.names(PH_vs_Macrophage)),
                        legend = "top",
                        font.legend = c("plain", 15),
                        ggtheme = ggplot2::theme_bw(),
                        top = 0) + font("xlab", size = 10) +
  font("ylab", size = 10) +
  font("xy.text", size = 10)
ggsave("~/TFM/r/results/maplots/4F/PH_vs_Macrophage_fungal.png", plot = PH_vs_Macrophage_maplot, width = 8, height = 4, units = "in", dpi = 300)

Spent_vs_Macrophage_maplot <- ggmaplot(Spent_vs_Macrophage,
                        fdr = 0.01, fc = 2, size = 2,
                        palette = c("#B31B21", "#1465AC", "darkgray"),
                        genenames = as.vector(row.names(Spent_vs_Macrophage)),
                        legend = "top",
                        font.legend = c("plain", 15),
                        ggtheme = ggplot2::theme_bw(),
                        top = 0) + font("xlab", size = 10) +
  font("ylab", size = 10) +
  font("xy.text", size = 10)
ggsave("~/TFM/r/results/maplots/4F/Spent_vs_Macrophage_fungal.png", plot = Spent_vs_Macrophage_maplot, width = 8, height = 4, units = "in", dpi = 300)

Spent_H2O2_vs_Macrophage_maplot <- ggmaplot(Spent_H2O2_vs_Macrophage,
                        fdr = 0.01, fc = 2, size = 2,
                        palette = c("#B31B21", "#1465AC", "darkgray"),
                        genenames = as.vector(row.names(Spent_H2O2_vs_Macrophage)),
                        legend = "top",
                        font.legend = c("plain", 15),
                        ggtheme = ggplot2::theme_bw(),
                        top = 0) + font("xlab", size = 10) +
  font("ylab", size = 10) +
  font("xy.text", size = 10)
ggsave("~/TFM/r/results/maplots/4F/Spent_H2O2_vs_Macrophage_fungal.png", plot = Spent_H2O2_vs_Macrophage_maplot, width = 8, height = 4, units = "in", dpi = 300)

PH_vs_H2O2_maplot <- ggmaplot(PH_vs_H2O2,
                        fdr = 0.01, fc = 2, size = 2,
                        palette = c("#B31B21", "#1465AC", "darkgray"),
                        genenames = as.vector(row.names(PH_vs_H2O2)),
                        legend = "top",
                        font.legend = c("plain", 15),
                        ggtheme = ggplot2::theme_bw(),
                        top = 0) + font("xlab", size = 10) +
  font("ylab", size = 10) +
  font("xy.text", size = 10)
ggsave("~/TFM/r/results/maplots/4F/PH_vs_H2O2_fungal.png", plot = PH_vs_H2O2_maplot, width = 8, height = 4, units = "in", dpi = 300)

Spent_vs_H2O2_maplot <- ggmaplot(Spent_vs_H2O2,
                              fdr = 0.01, fc = 2, size = 2,
                              palette = c("#B31B21", "#1465AC", "darkgray"),
                              genenames = as.vector(row.names(Spent_vs_H2O2)),
                              legend = "top",
                              font.legend = c("plain", 15),
                              ggtheme = ggplot2::theme_bw(),
                              top = 0) + font("xlab", size = 10) +
  font("ylab", size = 10) +
  font("xy.text", size = 10)
ggsave("~/TFM/r/results/maplots/4F/Spent_vs_H2O2_fungal.png", plot = Spent_vs_H2O2_maplot, width = 8, height = 4, units = "in", dpi = 300)

Spent_H2O2_vs_H2O2_maplot <- ggmaplot(Spent_H2O2_vs_H2O2,
                              fdr = 0.01, fc = 2, size = 2,
                              palette = c("#B31B21", "#1465AC", "darkgray"),
                              genenames = as.vector(row.names(Spent_H2O2_vs_H2O2)),
                              legend = "top",
                              font.legend = c("plain", 15),
                              ggtheme = ggplot2::theme_bw(),
                              top = 0) + font("xlab", size = 10) +
  font("ylab", size = 10) +
  font("xy.text", size = 10)
ggsave("~/TFM/r/results/maplots/4F/Spent_H2O2_vs_H2O2_fungal.png", plot = Spent_H2O2_vs_H2O2_maplot, width = 8, height = 4, units = "in", dpi = 300)

Spent_vs_PH_maplot <- ggmaplot(Spent_vs_PH,
                                      fdr = 0.01, fc = 2, size = 2,
                                      palette = c("#B31B21", "#1465AC", "darkgray"),
                                      genenames = as.vector(row.names(Spent_vs_PH)),
                                      legend = "top",
                                      font.legend = c("plain", 15),
                                      ggtheme = ggplot2::theme_bw(),
                                      top = 0) + font("xlab", size = 10) +
  font("ylab", size = 10) +
  font("xy.text", size = 10)
ggsave("~/TFM/r/results/maplots/4F/Spent_vs_PH_fungal.png", plot = Spent_vs_PH_maplot, width = 8, height = 4, units = "in", dpi = 300)

Spent_H2O2_vs_PH_maplot <- ggmaplot(Spent_H2O2_vs_PH,
                                      fdr = 0.01, fc = 2, size = 2,
                                      palette = c("#B31B21", "#1465AC", "darkgray"),
                                      genenames = as.vector(row.names(Spent_H2O2_vs_PH)),
                                      legend = "top",
                                      font.legend = c("plain", 15),
                                      ggtheme = ggplot2::theme_bw(),
                                      top = 0) + font("xlab", size = 10) +
  font("ylab", size = 10) +
  font("xy.text", size = 10)
ggsave("~/TFM/r/results/maplots/4F/Spent_H2O2_vs_PH_fungal.png", plot = Spent_H2O2_vs_PH_maplot, width = 8, height = 4, units = "in", dpi = 300)

Spent_H2O2_vs_Spent_maplot <- ggmaplot(Spent_H2O2_vs_Spent,
                                    fdr = 0.01, fc = 2, size = 2,
                                    palette = c("#B31B21", "#1465AC", "darkgray"),
                                    genenames = as.vector(row.names(Spent_H2O2_vs_Spent)),
                                    legend = "top",
                                    font.legend = c("plain", 15),
                                    ggtheme = ggplot2::theme_bw(),
                                    top = 0) + font("xlab", size = 10) +
  font("ylab", size = 10) +
  font("xy.text", size = 10)
ggsave("~/TFM/r/results/maplots/4F/Spent_H2O2_vs_Spent_fungal.png", plot = Spent_H2O2_vs_Spent_maplot, width = 8, height = 4, units = "in", dpi = 300)

###

all_up<-rbind.data.frame(H2O2_vs_Macrophage_up_df,
                         PH_vs_Macrophage_up_df,
                         Spent_vs_Macrophage_up_df,
                         Spent_H2O2_vs_Macrophage_up_df)
rownames(all_up)<-NULL

all_up$Cluster <- ifelse(grepl("PH", all_up$Cluster), "PH vs Macrophage", 
                         ifelse(grepl("Spent_H2O2",all_up$Cluster),"Spent-H202 vs Macrophage", 
                                ifelse(grepl("Spent_",all_up$Cluster),"Spent vs Macrophage","H2O2 vs Macrophage")))


all_down <-rbind.data.frame(H2O2_vs_Macrophage_down_df,
                            PH_vs_Macrophage_down_df,
                            Spent_vs_Macrophage_down_df,
                            Spent_H2O2_vs_Macrophage_down_df)
rownames(all_down)<-NULL

all_down$Cluster <- ifelse(grepl("PH", all_down$Cluster), "PH vs Macrophage", 
                           ifelse(grepl("Spent_H2O2",all_down$Cluster),"Spent-H202 vs Macrophage", 
                                  ifelse(grepl("Spent_",all_down$Cluster),"Spent vs Macrophage","H2O2 vs Macrophage")))


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
  ggsave("~/TFM/r/results/GO/4F/up_GO_general.png",units="in", width=14, height=7, dpi=300)



## down

to_hack_down@compareClusterResult <- all_down
plot_down<-dotplot(to_hack_down, showCategory=15)
plot_down+ggplot2::theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))+
  # ggsave("./results/project_4/down_GO_general.pdf",units="in", width=10, height=7, dpi=300)
  ggsave("~/TFM/r/results/GO/4F/down_GO_general.png",units="in", width=10, height=7, dpi=300)




####


### Creating an excel file with DE data

wb <- createWorkbook("Project_4_excel")

header<-c("gene", "baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")
file_names<-c("H2O2_vs_Macrophage",
              "PH_vs_Macrophage",
              "Spent_vs_Macrophage",
              "Spent_H2O2_vs_Macrophage")



descriptions<-c("C. glabrata in H2O2 vs macrophage",
                "C. glabrata in PH vs macrophage",
                "C. glabrata in Spent vs macrophage",
                "C. glabrata in Spent_H2O2 vs macrophage")

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
saveWorkbook(wb, "~/TFM/r/results/excel/4F/Project_4_excel.xlsx", overwrite = TRUE)
