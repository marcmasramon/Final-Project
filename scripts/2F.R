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
for (f in c("Initial-A_ReadsPerGene.out.tab","Initial-B_ReadsPerGene.out.tab","Initial-C_ReadsPerGene.out.tab", "Initial-D_ReadsPerGene.out.tab", 
            "R-A-3hr_ReadsPerGene.out.tab", "R-B-3hr_ReadsPerGene.out.tab", "R-C-3hr_ReadsPerGene.out.tab", "R-D-3hr_ReadsPerGene.out.tab", 
            "R-A-24hr_ReadsPerGene.out.tab", "R-B-24hr_ReadsPerGene.out.tab", "R-C-24hr_ReadsPerGene.out.tab", "R-D-24hr_ReadsPerGene.out.tab"))
{
  file<-read.table(f, header = TRUE, sep = "\t")
  fungal_controls<-cbind(fungal_controls,file[,2])
}

rownames(fungal_controls)<-file[,1]
colnames(fungal_controls)<-c("Initial-A","Initial-B","Initial-C", "Initial-D", 
                  "R-A-3hr", "R-B-3hr", "R-C-3hr", "R-D-3hr", 
                  "R-A-24hr", "R-B-24hr", "R-C-24hr", "R-D-24hr")
fungal_controls<-fungal_controls[-c(1,2,3),]


### infection

infection<-NULL
for (f in c("g-A-I-3hr_ReadsPerGene.out.tab", "g-B-I-3hr_ReadsPerGene.out.tab", "g-C-I-3hr_ReadsPerGene.out.tab", "g-D-I-3hr_ReadsPerGene.out.tab",
            "g-A-II-3hr_ReadsPerGene.out.tab", "g-B-II-3hr_ReadsPerGene.out.tab", "g-C-II-3hr_ReadsPerGene.out.tab", "g-D-II-3hr_ReadsPerGene.out.tab",
            "g-A-III-3hr_ReadsPerGene.out.tab", "g-B-III-3hr_ReadsPerGene.out.tab", "g-C-III-3hr_ReadsPerGene.out.tab", "g-D-III-3hr_ReadsPerGene.out.tab",
            "g-A-I-24hr_ReadsPerGene.out.tab", "g-B-I-24hr_ReadsPerGene.out.tab", "g-C-I-24hr_ReadsPerGene.out.tab", "g-D-I-24hr_ReadsPerGene.out.tab",
            "g-A-II-24hr_ReadsPerGene.out.tab", "g-B-II-24hr_ReadsPerGene.out.tab", "g-C-II-24hr_ReadsPerGene.out.tab", "g-D-II-24hr_ReadsPerGene.out.tab", 
            "g-A-III-24hr_ReadsPerGene.out.tab", "g-B-III-24hr_ReadsPerGene.out.tab", "g-C-III-24hr_ReadsPerGene.out.tab", "g-D-III-24hr_ReadsPerGene.out.tab"))
{
  file<-read.table(f, header = TRUE, sep = "\t")
  infection<-cbind(infection,file[,2])
}


rownames(infection)<-file[,1]
colnames(infection)<-c("g-A-I-3hr", "g-B-I-3hr", "g-C-I-3hr", "g-D-I-3hr",
                   "g-A-II-3hr", "g-B-II-3hr", "g-C-II-3hr", "g-D-II-3hr",
                   "g-A-III-3hr", "g-B-III-3hr", "g-C-III-3hr", "g-D-III-3hr",
                   "g-A-I-24hr", "g-B-I-24hr", "g-C-I-24hr", "g-D-I-24hr",
                   "g-A-II-24hr", "g-B-II-24hr", "g-C-II-24hr", "g-D-II-24hr", 
                   "g-A-III-24hr", "g-B-III-24hr", "g-C-III-24hr", "g-D-III-24hr")
infection<-infection[-c(1,2,3),]

infection_fungus <- infection[rownames(infection) %in% rownames(fungal_controls),]


### all fungal data

fungal_data <- merge(infection_fungus,fungal_controls, by="row.names")
rownames(fungal_data)<- fungal_data$Row.names
fungal_data <- fungal_data[,-c(1)]

### PCA plot for all data

Strain <- factor(rep(c("A","B","C","D"),9), 
                  levels=c("A","B","C", "D"))

Time <- factor(c("3hr", "3hr", "3hr", "3hr",
                      "3hr", "3hr", "3hr", "3hr",
                      "3hr", "3hr", "3hr", "3hr",
                      "24hr", "24hr", "24hr", "24hr",
                      "24hr", "24hr", "24hr", "24hr", 
                      "24hr", "24hr", "24hr", "24hr",
                  "Initial","Initial","Initial", "Initial", 
                  "3hr", "3hr", "3hr", "3hr", 
                  "24hr", "24hr", "24hr", "24hr"), 
                  levels=c("Initial","3hr","24hr"))

Media <- factor(c(rep(c("infection"),24),
                    rep(c("initial"),4),
                    rep(c("RPMI"),8)),
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
  ggsave("~/TFM/r/results/pcaplots/2F/pca.png", width = 8, height = 6, units = "in", dpi = 300)


################
#### DESeq2 ####
################


fungal_data_for_deseq2 <- fungal_data[,-c(25:36)]
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
pre_ddsg <- DESeqDataSetFromMatrix(fungal_data_for_deseq2, colDatag, design = ~time+strain)#variable of interest is the second one
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

############
### test


as.data.frame(B_vs_A["CAGL0A00451g",])

A_vs_B <- results(ddsg, contrast=c("strain","A","B"))
as.data.frame(A_vs_B["CAGL0A00451g",])


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
    ggplot2::ggsave(sprintf("./results/2F/%s_%s.png",name,ont),
                    units="in", width=10, heigh=7, dpi=600)
    write.table(ego,sprintf("./results/2F/%s_%s.txt",name,ont), sep="\t")
    ego_df <- as.data.frame(ego)
    ego_df <- cbind.data.frame("Cluster"=name,ego_df)
    assign(paste0(substring(name,4),"_df"),as.data.frame(ego_df), envir = .GlobalEnv)
    
  }
  
}

###

setwd("~/TFM/r/results")

fungal_goterms(B_vs_A_up, "2F_B_vs_A_up", "BP")
fungal_goterms(B_vs_A_up, "2F_B_vs_A_up", "MF")
fungal_goterms(B_vs_A_up, "2F_B_vs_A_up", "CC")

fungal_goterms(C_vs_A_up, "2F_C_vs_A_up", "BP")
fungal_goterms(C_vs_A_up, "2F_C_vs_A_up", "MF")
fungal_goterms(C_vs_A_up, "2F_C_vs_A_up", "CC")

fungal_goterms(D_vs_A_up, "2F_D_vs_A_up", "BP")
fungal_goterms(D_vs_A_up, "2F_D_vs_A_up", "MF")
fungal_goterms(D_vs_A_up, "2F_D_vs_A_up", "CC")

###

fungal_goterms(B_vs_A_down, "2F_B_vs_A_down", "BP")
fungal_goterms(B_vs_A_down, "2F_B_vs_A_down", "MF")
fungal_goterms(B_vs_A_down, "2F_B_vs_A_down", "CC")

fungal_goterms(C_vs_A_down, "2F_C_vs_A_down", "BP")
fungal_goterms(C_vs_A_down, "2F_C_vs_A_down", "MF")
fungal_goterms(C_vs_A_down, "2F_C_vs_A_down", "CC")

fungal_goterms(D_vs_A_down, "2F_D_vs_A_down", "BP")
fungal_goterms(D_vs_A_down, "2F_D_vs_A_down", "MF")
fungal_goterms(D_vs_A_down, "2F_D_vs_A_down", "CC")

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
ggsave("~/TFM/r/results/maplots/2F/B_vs_A_fungal.png", plot = B_vs_A_maplot, width = 8, height = 4, units = "in", dpi = 300)

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
ggsave("~/TFM/r/results/maplots/2F/C_vs_A_fungal.png", plot = C_vs_A_maplot, width = 8, height = 4, units = "in", dpi = 300)

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
ggsave("~/TFM/r/results/maplots/2F/D_vs_A_fungal.png", plot = D_vs_A_maplot, width = 8, height = 4, units = "in", dpi = 300)

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
ggsave("~/TFM/r/results/maplots/2F/C_vs_B_fungal.png", plot = C_vs_B_maplot, width = 8, height = 4, units = "in", dpi = 300)

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
ggsave("~/TFM/r/results/maplots/2F/D_vs_B_fungal.png", plot = D_vs_B_maplot, width = 8, height = 4, units = "in", dpi = 300)

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
ggsave("~/TFM/r/results/maplots/2F/D_vs_C_fungal.png", plot = D_vs_C_maplot, width = 8, height = 4, units = "in", dpi = 300)

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
  ggsave("~/TFM/r/results/GO/2F/up_GO_general.png",units="in", width=14, height=7, dpi=300)



## down

to_hack_down@compareClusterResult <- all_down
plot_down<-dotplot(to_hack_down, showCategory=15)
plot_down+ggplot2::theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))+
  # ggsave("./results/project_4/down_GO_general.pdf",units="in", width=10, height=7, dpi=300)
  ggsave("~/TFM/r/results/GO/2F/down_GO_general.png",units="in", width=10, height=7, dpi=300)




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
saveWorkbook(wb, "~/TFM/r/results/excel/2F/Project_2_excel.xlsx", overwrite = TRUE)
