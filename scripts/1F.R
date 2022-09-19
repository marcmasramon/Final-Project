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
for (f in c("Initial-A_ReadsPerGene.out.tab","Initial-H_ReadsPerGene.out.tab","Initial-I_ReadsPerGene.out.tab","Initial-J_ReadsPerGene.out.tab", "Initial-K_ReadsPerGene.out.tab", 
            "R-A-3hr_ReadsPerGene.out.tab","R-H-3hr_ReadsPerGene.out.tab", "R-I-3hr_ReadsPerGene.out.tab", "R-J-3hr_ReadsPerGene.out.tab", "R-K-3hr_ReadsPerGene.out.tab", 
            "R-A-24hr_ReadsPerGene.out.tab","R-H-24hr_ReadsPerGene.out.tab", "R-I-24hr_ReadsPerGene.out.tab", "R-J-24hr_ReadsPerGene.out.tab", "R-K-24hr_ReadsPerGene.out.tab"))
{
  file<-read.table(f, header = TRUE, sep = "\t")
  fungal_controls<-cbind(fungal_controls,file[,2])
}

rownames(fungal_controls)<-file[,1]
colnames(fungal_controls)<-c("Initial-A","Initial-H","Initial-I","Initial-J", "Initial-K", 
                             "R-A-3hr","R-H-3hr", "R-I-3hr", "R-J-3hr", "R-K-3hr", 
                             "R-A-24hr","R-H-24hr", "R-I-24hr", "R-J-24hr", "R-K-24hr")
fungal_controls<-fungal_controls[-c(1,2,3),]


### infection

infection<-NULL
for (f in c("g-A-I-3hr_ReadsPerGene.out.tab","g-H-I-3hr_ReadsPerGene.out.tab", "g-I-I-3hr_ReadsPerGene.out.tab", "g-J-I-3hr_ReadsPerGene.out.tab", "g-K-I-3hr_ReadsPerGene.out.tab",
            "g-A-II-3hr_ReadsPerGene.out.tab","g-H-II-3hr_ReadsPerGene.out.tab", "g-I-II-3hr_ReadsPerGene.out.tab", "g-J-II-3hr_ReadsPerGene.out.tab", "g-K-II-3hr_ReadsPerGene.out.tab",
            "g-A-III-3hr_ReadsPerGene.out.tab","g-H-III-3hr_ReadsPerGene.out.tab", "g-I-III-3hr_ReadsPerGene.out.tab", "g-J-III-3hr_ReadsPerGene.out.tab", "g-K-III-3hr_ReadsPerGene.out.tab",
            "g-A-I-24hr_ReadsPerGene.out.tab","g-H-I-24hr_ReadsPerGene.out.tab", "g-I-I-24hr_ReadsPerGene.out.tab", "g-J-I-24hr_ReadsPerGene.out.tab", "g-K-I-24hr_ReadsPerGene.out.tab",
            "g-A-II-24hr_ReadsPerGene.out.tab","g-H-II-24hr_ReadsPerGene.out.tab", "g-I-II-24hr_ReadsPerGene.out.tab", "g-J-II-24hr_ReadsPerGene.out.tab", "g-K-II-24hr_ReadsPerGene.out.tab", 
            "g-A-III-24hr_ReadsPerGene.out.tab","g-H-III-24hr_ReadsPerGene.out.tab", "g-I-III-24hr_ReadsPerGene.out.tab", "g-J-III-24hr_ReadsPerGene.out.tab", "g-K-III-24hr_ReadsPerGene.out.tab")) 
 {
  file<-read.table(f, header = TRUE, sep = "\t")
  infection<-cbind(infection,file[,2])
}


rownames(infection)<-file[,1]
colnames(infection)<-c("g-A-I-3hr","g-H-I-3hr", "g-I-I-3hr", "g-J-I-3hr", "g-K-I-3hr",
                       "g-A-II-3hr","g-H-II-3hr", "g-I-II-3hr", "g-J-II-3hr", "g-K-II-3hr",
                       "g-A-III-3hr","g-H-III-3hr", "g-I-III-3hr", "g-J-III-3hr", "g-K-III-3hr",
                       "g-A-I-24hr","g-H-I-24hr", "g-I-I-24hr", "g-J-I-24hr", "g-K-I-24hr",
                       "g-A-II-24hr","g-H-II-24hr", "g-I-II-24hr", "g-J-II-24hr", "g-K-II-24hr", 
                       "g-A-III-24hr","g-H-III-24hr", "g-I-III-24hr", "g-J-III-24hr", "g-K-III-24hr")
infection<-infection[-c(1,2,3),]

infection_fungus <- infection[rownames(infection) %in% rownames(fungal_controls),]


### all fungal data

fungal_data <- merge(infection_fungus,fungal_controls, by="row.names")
rownames(fungal_data)<- fungal_data$Row.names
fungal_data <- fungal_data[,-c(1)]

### PCA plot for all data

Strain <- factor(rep(c("A","H","I","J","K"),9), 
                 levels=c("A","H","I","J","K"))

Time <- factor(c("3hr", "3hr", "3hr", "3hr","3hr",
                 "3hr", "3hr", "3hr", "3hr","3hr",
                 "3hr", "3hr", "3hr", "3hr","3hr",
                 "24hr", "24hr", "24hr", "24hr","24hr",
                 "24hr", "24hr", "24hr", "24hr", "24hr",
                 "24hr", "24hr", "24hr", "24hr","24hr",
                 "Initial","Initial","Initial", "Initial", "Initial",
                 "3hr", "3hr", "3hr", "3hr", "3hr",
                 "24hr", "24hr", "24hr", "24hr","24hr"), 
               levels=c("Initial","3hr","24hr"))

Media <- factor(c(rep(c("infection"),30),
                  rep(c("initial"),5),
                  rep(c("RPMI"),10)),
                levels = c("infection","initial","RPMI"))

colData_pca <- DataFrame(Strain=Strain,Time=Time,Media=Media)
colData_pca$Condition <- paste0(colData_pca$Time,"_",colData_pca$Media)


pre_dds_pca <- DESeqDataSetFromMatrix(fungal_data, colData_pca, design = ~Time)#just a fake desing, does not affect the data normalization and transforamtion
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
  ggsave("~/TFM/r/results/pcaplots/1F/pca.png", width = 8, height = 6, units = "in", dpi = 300)

################
#### DESeq2 ####
################


fungal_data_for_deseq2 <- fungal_data[,-c(31:45)]
group1g <- factor(x=c("A","H","I","J","K",
                      "A","H","I","J","K",
                      "A","H","I","J","K",
                      "A","H","I","J","K",
                      "A","H","I","J","K", 
                      "A","H","I","J","K"),
                  levels=c("A","H","I","J","K"))
group2g <- factor(x=c("3hr", "3hr", "3hr", "3hr","3hr",
                      "3hr", "3hr", "3hr", "3hr","3hr",
                      "3hr", "3hr", "3hr", "3hr","3hr",
                      "24hr", "24hr", "24hr", "24hr","24hr",
                      "24hr", "24hr", "24hr", "24hr", "24hr",
                      "24hr", "24hr", "24hr", "24hr","24hr"),
                  levels=c("3hr","24hr"))
colDatag <- DataFrame(strain=group1g,time=group2g)
pre_ddsg <- DESeqDataSetFromMatrix(fungal_data_for_deseq2, colDatag, design = ~time+strain)#variable of interest is the second one
ddsg <- DESeq(pre_ddsg)

H_vs_A <- results(ddsg, contrast = c("strain","H","A"))
I_vs_A <- results(ddsg, contrast = c("strain","I","A"))
J_vs_A <- results(ddsg, contrast = c("strain","J","A"))
K_vs_A <- results(ddsg, contrast = c("strain","K","A"))

I_vs_H <- results(ddsg, contrast = c("strain","I","H"))
J_vs_H <- results(ddsg, contrast = c("strain","J","H"))
K_vs_H <- results(ddsg, contrast = c("strain","K","H"))

J_vs_I <- results(ddsg, contrast = c("strain","J","I"))
K_vs_I <- results(ddsg, contrast = c("strain","K","I"))

K_vs_J <- results(ddsg, contrast = c("strain","K","J"))

###

get_up_down_reg <- function(res, name){
  res <- as.data.frame(res)
  up <- row.names(res[!is.na(res$padj) & res$padj<0.01 & res$log2FoldChange > 1,])
  down <- row.names(res[!is.na(res$padj) & res$padj<0.01 & res$log2FoldChange < -1,])
  
  assign(paste0(name,"_up"),up,envir = .GlobalEnv)
  assign(paste0(name,"_down"),down,envir = .GlobalEnv)
  
}

###

get_up_down_reg(H_vs_A,"H_vs_A")
get_up_down_reg(I_vs_A,"I_vs_A")
get_up_down_reg(J_vs_A,"J_vs_A")
get_up_down_reg(K_vs_A,"K_vs_A")

length(H_vs_A_up)
length(H_vs_A_down)

length(I_vs_A_up)
length(I_vs_A_down)

length(J_vs_A_up)
length(J_vs_A_down)

length(K_vs_A_up)
length(K_vs_A_down)

###

H_vs_A_up <- row.names(H_vs_A[!is.na(H_vs_A$padj) & H_vs_A$padj<0.01 & H_vs_A$log2FoldChange > 1,])
H_vs_A_down <- row.names(H_vs_A[!is.na(H_vs_A$padj) & H_vs_A$padj<0.01 & H_vs_A$log2FoldChange < -1,])

I_vs_A_up <- row.names(I_vs_A[!is.na(I_vs_A$padj) & I_vs_A$padj<0.01 & I_vs_A$log2FoldChange > 1,])
I_vs_A_down <- row.names(I_vs_A[!is.na(I_vs_A$padj) & I_vs_A$padj<0.01 & I_vs_A$log2FoldChange < -1,])

J_vs_A_up <- row.names(J_vs_A[!is.na(J_vs_A$padj) & J_vs_A$padj<0.01 & J_vs_A$log2FoldChange > 1,])
J_vs_A_down <- row.names(J_vs_A[!is.na(J_vs_A$padj) & J_vs_A$padj<0.01 & J_vs_A$log2FoldChange < -1,])

K_vs_A_up <- row.names(K_vs_A[!is.na(K_vs_A$padj) & K_vs_A$padj<0.01 & K_vs_A$log2FoldChange > 1,])
K_vs_A_down <- row.names(K_vs_A[!is.na(K_vs_A$padj) & K_vs_A$padj<0.01 & K_vs_A$log2FoldChange < -1,])

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
    ggplot2::ggsave(sprintf("./results/1F/%s_%s.png",name,ont),
                    units="in", width=10, heigh=7, dpi=600)
    write.table(ego,sprintf("./results/1F/%s_%s.txt",name,ont), sep="\t")
    ego_df <- as.data.frame(ego)
    ego_df <- cbind.data.frame("Cluster"=name,ego_df)
    assign(paste0(substring(name,4),"_df"),as.data.frame(ego_df), envir = .GlobalEnv)
    
  }
  
}

###

setwd("~/TFM/r/results")

fungal_goterms(H_vs_A_up, "1F_H_vs_A_up", "BP")
fungal_goterms(H_vs_A_up, "1F_H_vs_A_up", "MF")
fungal_goterms(H_vs_A_up, "1F_H_vs_A_up", "CC")

fungal_goterms(I_vs_A_up, "1F_I_vs_A_up", "BP")
fungal_goterms(I_vs_A_up, "1F_I_vs_A_up", "MF")
fungal_goterms(I_vs_A_up, "1F_I_vs_A_up", "CC")

fungal_goterms(J_vs_A_up, "1F_J_vs_A_up", "BP")
fungal_goterms(J_vs_A_up, "1F_J_vs_A_up", "MF")
fungal_goterms(J_vs_A_up, "1F_J_vs_A_up", "CC")

fungal_goterms(K_vs_A_up, "1F_K_vs_A_up", "BP")
fungal_goterms(K_vs_A_up, "1F_K_vs_A_up", "MF")
fungal_goterms(K_vs_A_up, "1F_K_vs_A_up", "CC")

###

fungal_goterms(H_vs_A_down, "1F_H_vs_A_down", "BP")
fungal_goterms(H_vs_A_down, "1F_H_vs_A_down", "MF")
fungal_goterms(H_vs_A_down, "1F_H_vs_A_down", "CC")

fungal_goterms(I_vs_A_down, "1F_I_vs_A_down", "BP")
fungal_goterms(I_vs_A_down, "1F_I_vs_A_down", "MF")
fungal_goterms(I_vs_A_down, "1F_I_vs_A_down", "CC")

fungal_goterms(J_vs_A_down, "1F_J_vs_A_down", "BP")
fungal_goterms(J_vs_A_down, "1F_J_vs_A_down", "MF")
fungal_goterms(J_vs_A_down, "1F_J_vs_A_down", "CC")

fungal_goterms(K_vs_A_down, "1F_K_vs_A_down", "BP")
fungal_goterms(K_vs_A_down, "1F_K_vs_A_down", "MF")
fungal_goterms(K_vs_A_down, "1F_K_vs_A_down", "CC")

###

H_vs_A_maplot <- ggmaplot(H_vs_A,
                             fdr = 0.01, fc = 2, size = 2,
                             palette = c("#B31B21", "#1465AC", "darkgray"),
                             genenames = as.vector(row.names(H_vs_A)),
                             legend = "top",
                             font.legend = c("plain", 15),
                             ggtheme = ggplot2::theme_bw(),
                             top = 0) + font("xlab", size = 10) +
  font("ylab", size = 10) +
  font("xy.text", size = 10)
ggsave("~/TFM/r/results/maplots/1F/H_vs_A_fungal.png", plot = H_vs_A_maplot, width = 8, height = 4, units = "in", dpi = 300)

I_vs_A_maplot <- ggmaplot(I_vs_A,
                           fdr = 0.01, fc = 2, size = 2,
                           palette = c("#B31B21", "#1465AC", "darkgray"),
                           genenames = as.vector(row.names(I_vs_A)),
                           legend = "top",
                           font.legend = c("plain", 15),
                           ggtheme = ggplot2::theme_bw(),
                           top = 0) + font("xlab", size = 10) +
  font("ylab", size = 10) +
  font("xy.text", size = 10)
ggsave("~/TFM/r/results/maplots/1F/I_vs_A_fungal.png", plot = I_vs_A_maplot, width = 8, height = 4, units = "in", dpi = 300)

J_vs_A_maplot <- ggmaplot(J_vs_A,
                              fdr = 0.01, fc = 2, size = 2,
                              palette = c("#B31B21", "#1465AC", "darkgray"),
                              genenames = as.vector(row.names(J_vs_A)),
                              legend = "top",
                              font.legend = c("plain", 15),
                              ggtheme = ggplot2::theme_bw(),
                              top = 0) + font("xlab", size = 10) +
  font("ylab", size = 10) +
  font("xy.text", size = 10)
ggsave("~/TFM/r/results/maplots/1F/J_vs_A_fungal.png", plot = J_vs_A_maplot, width = 8, height = 4, units = "in", dpi = 300)

K_vs_A_maplot <- ggmaplot(K_vs_A,
                                   fdr = 0.01, fc = 2, size = 2,
                                   palette = c("#B31B21", "#1465AC", "darkgray"),
                                   genenames = as.vector(row.names(K_vs_A)),
                                   legend = "top",
                                   font.legend = c("plain", 15),
                                   ggtheme = ggplot2::theme_bw(),
                                   top = 0) + font("xlab", size = 10) +
  font("ylab", size = 10) +
  font("xy.text", size = 10)
ggsave("~/TFM/r/results/maplots/1F/K_vs_A_fungal.png", plot = K_vs_A_maplot, width = 8, height = 4, units = "in", dpi = 300)

I_vs_H_maplot <- ggmaplot(I_vs_H,
                              fdr = 0.01, fc = 2, size = 2,
                              palette = c("#B31B21", "#1465AC", "darkgray"),
                              genenames = as.vector(row.names(I_vs_H)),
                              legend = "top",
                              font.legend = c("plain", 15),
                              ggtheme = ggplot2::theme_bw(),
                              top = 0) + font("xlab", size = 10) +
  font("ylab", size = 10) +
  font("xy.text", size = 10)
ggsave("~/TFM/r/results/maplots/1F/I_vs_H_fungal.png", plot = I_vs_H_maplot, width = 8, height = 4, units = "in", dpi = 300)

J_vs_H_maplot <- ggmaplot(J_vs_H,
                                 fdr = 0.01, fc = 2, size = 2,
                                 palette = c("#B31B21", "#1465AC", "darkgray"),
                                 genenames = as.vector(row.names(J_vs_H)),
                                 legend = "top",
                                 font.legend = c("plain", 15),
                                 ggtheme = ggplot2::theme_bw(),
                                 top = 0) + font("xlab", size = 10) +
  font("ylab", size = 10) +
  font("xy.text", size = 10)
ggsave("~/TFM/r/results/maplots/1F/J_vs_H_fungal.png", plot = J_vs_H_maplot, width = 8, height = 4, units = "in", dpi = 300)

K_vs_H_maplot <- ggmaplot(K_vs_H,
                                      fdr = 0.01, fc = 2, size = 2,
                                      palette = c("#B31B21", "#1465AC", "darkgray"),
                                      genenames = as.vector(row.names(K_vs_H)),
                                      legend = "top",
                                      font.legend = c("plain", 15),
                                      ggtheme = ggplot2::theme_bw(),
                                      top = 0) + font("xlab", size = 10) +
  font("ylab", size = 10) +
  font("xy.text", size = 10)
ggsave("~/TFM/r/results/maplots/1F/K_vs_H_fungal.png", plot = K_vs_H_maplot, width = 8, height = 4, units = "in", dpi = 300)

J_vs_I_maplot <- ggmaplot(J_vs_I,
                               fdr = 0.01, fc = 2, size = 2,
                               palette = c("#B31B21", "#1465AC", "darkgray"),
                               genenames = as.vector(row.names(J_vs_I)),
                               legend = "top",
                               font.legend = c("plain", 15),
                               ggtheme = ggplot2::theme_bw(),
                               top = 0) + font("xlab", size = 10) +
  font("ylab", size = 10) +
  font("xy.text", size = 10)
ggsave("~/TFM/r/results/maplots/1F/J_vs_I_fungal.png", plot = J_vs_I_maplot, width = 8, height = 4, units = "in", dpi = 300)

K_vs_I_maplot <- ggmaplot(K_vs_I,
                                    fdr = 0.01, fc = 2, size = 2,
                                    palette = c("#B31B21", "#1465AC", "darkgray"),
                                    genenames = as.vector(row.names(K_vs_I)),
                                    legend = "top",
                                    font.legend = c("plain", 15),
                                    ggtheme = ggplot2::theme_bw(),
                                    top = 0) + font("xlab", size = 10) +
  font("ylab", size = 10) +
  font("xy.text", size = 10)
ggsave("~/TFM/r/results/maplots/1F/K_vs_I_fungal.png", plot = K_vs_I_maplot, width = 8, height = 4, units = "in", dpi = 300)

K_vs_J_maplot <- ggmaplot(K_vs_J,
                                       fdr = 0.01, fc = 2, size = 2,
                                       palette = c("#B31B21", "#1465AC", "darkgray"),
                                       genenames = as.vector(row.names(K_vs_J)),
                                       legend = "top",
                                       font.legend = c("plain", 15),
                                       ggtheme = ggplot2::theme_bw(),
                                       top = 0) + font("xlab", size = 10) +
  font("ylab", size = 10) +
  font("xy.text", size = 10)
ggsave("~/TFM/r/results/maplots/1F/K_vs_J_fungal.png", plot = K_vs_J_maplot, width = 8, height = 4, units = "in", dpi = 300)

###

all_up<-rbind.data.frame(H_vs_A_up_df,
                         I_vs_A_up_df,
                         J_vs_A_up_df,
                         K_vs_A_up_df)
rownames(all_up)<-NULL

all_up$Cluster <- ifelse(grepl("H", all_up$Cluster), "H vs A", 
                         ifelse(grepl("I",all_up$Cluster),"I vs A", 
                                ifelse(grepl("J",all_up$Cluster),"J vs A","K vs A")))


all_down <-rbind.data.frame(H_vs_A_down_df,
                            I_vs_A_down_df,
                            J_vs_A_down_df,
                            K_vs_A_down_df)
rownames(all_down)<-NULL

all_down$Cluster <- ifelse(grepl("H", all_down$Cluster), "H vs A", 
                           ifelse(grepl("I",all_down$Cluster),"I vs A", 
                                  ifelse(grepl("J",all_down$Cluster),"J vs A","K vs A")))


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
  ggsave("~/TFM/r/results/GO/1F/up_GO_general.png",units="in", width=14, height=7, dpi=300)



## down

to_hack_down@compareClusterResult <- all_down
plot_down<-dotplot(to_hack_down, showCategory=15)
plot_down+ggplot2::theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))+
  # ggsave("./results/project_4/down_GO_general.pdf",units="in", width=10, height=7, dpi=300)
  ggsave("~/TFM/r/results/GO/1F/down_GO_general.png",units="in", width=10, height=7, dpi=300)




####


### Creating an excel file with DE data

wb <- createWorkbook("Project_1_excel")

header<-c("gene", "baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")
file_names<-c("H_vs_A",
              "I_vs_A",
              "J_vs_A",
              "K_vs_A")



descriptions<-c("deletant lacking FKS1 and derived from CBS138 vs CBS138",
                "deletant lacking FKS2 and derived from CBS138 vs CBS138",
                "echinocandin-resistant (ECR) isolate carrying F659del and derived from CBS138 vs CBS138",
                "ECR carrying F659V and derived fro CBS138 vs CBS138")

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
saveWorkbook(wb, "~/TFM/r/results/excel/1F/Project_1_excel.xlsx", overwrite = TRUE)
