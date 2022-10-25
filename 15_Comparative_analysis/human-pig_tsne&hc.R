BiocManager::install("Rtsne")
library(Rtsne)
library(ggplot2)
######tsne : reduce dim
load("/Users/baizhonghao/Downloads/human-pig GTEx/expressionSeurat.Rdata")
load("/Users/baizhonghao/Downloads/human-pig GTEx/DataMatrix20217.RData")
all_tissues<-unique(sort(experiment$Tissue_new))
all_tissues<-all_tissues[-8]
all_tissues
sample_sum<-NULL
for(i in 1:length(tissues)){
  sample<-experiment[which(experiment$Tissue_new == tissues[i]),]
  sample_sum<-rbind(sample_sum,sample)
}
experiment_tsne<-sample_sum
unique(sort(experiment_tsne$Tissue_new))
expression_inter<-expression_inter[match(experiment_tsne$SampleID,rownames(expression_inter)),]
expression_inter<-t(apply(expression_inter, MARGIN = 1, FUN = scale))
expression_tsne <- Rtsne(expression_inter,dims = 2, perplexity=150, theta=0.5, 
                         verbose=TRUE, max_iter = 1000,check_duplicates = FALSE,partial_pca = T,num_threads=50)

####plot clustering figure color by species
#tiff(file = "/Users/baizhonghao/Downloads/human-pig GTEx/TSNE_1-legend.tiff",##reqiured to change
#     res = 300, width = 2000, height = 1500,compression = "lzw")
par(mar=c(2,2,2,2))
info.norm<-as.data.frame(expression_tsne$Y)
colnames(info.norm)<-c("tsne1","tsne2")
colors<-c("blue","red")
#colors<-distinctColorPalette(length(unique(sort(experiment$species))))
names(colors)=c('Human',"Pig")
tsne1<-ggplot(info.norm, aes(x = tsne1, y = tsne2, colour = experiment_tsne$Species)) + 
  geom_point(alpha = 1,size=3) + theme_classic()+
  xlab("t-SNE-1")+ylab("t-SNE-2")+
  theme(panel.background=element_blank(),
        panel.grid.minor=element_blank(),panel.grid.major=element_blank(),
        legend.text = element_text(size=16),
        legend.title = element_text(size=18,face="bold"),
        #legend.position='none',
        axis.title.y=element_text(size=18,color="black"),
        axis.title.x=element_text(size=18,color="black"),
        axis.text.y=element_text(size=16,color="black"),
        axis.text.x=element_text(size=16,color="black",hjust=0.95,vjust=0.2))+
  scale_colour_manual(name="Species",values=colors)#+stat_ellipse(type="norm")
ggsave(tsne1,file='/Users/baizhonghao/Downloads/17tissues-tsne1.pdf',width=6,height=4,dpi=300)

###color by tissue

par(mar=c(2,2,5,5))
info.norm<-as.data.frame(expression_tsne$Y)
colnames(info.norm)<-c("tsne1","tsne2")

colors<-c("#CC66FF","#AAAAFF","#FF0000","#8EABD2","#FFD700","#33CCCC","#FDFDBF", "#8EA9DB","#8B0F55",
       "#E2EFDA","#7570B3","#AAEEFF","#99BB88","#FFDD99","#FF6600","#A6CEE3","#A6761D","#ed1299", "#09f9f5", "#246b93", "#cc8e12", "#d561dd", "#c93f00", "#ddd53e",
       "#4aef7b", "#e86502", "#9ed84e", "#39ba30", "#6ad157", "#8249aa", "#99db27", "#e07233", "#ff523f",
       "#ce2523", "#f7aa5d", "#cebb10", "#03827f", "#931635", "#373bbf", "#a1ce4c", "#ef3bb6", "#d66551",
       "#1a918f", "#ff66fc", "#2927c4", "#7149af")
names(colors)=c("Adipose","Artery","Blood","Colon","Frontal_cortex","Heart","Hypothalamus","Ileum",
              "Kidney","Liver","Lung","Muscle","Ovary", "Pituitary","Spleen","Testis","Uterus","Adrenal Gland","Bladder","Blastocyst","Blastomere","Cartilage",
              "Cervix Uteri","Duodenum","Embryo","Esophagus","Fallopian Tube","Fetal_thymus","Jejunum","Lymph_node","Macrophage",       
              "Milk","Morula","Nerve","Oocyte","Pancreas","pEPSC","Placenta","Prostate","Salivary Gland","Skin","Stomach","Synovial_membrane","Thyroid","Vagina","Breast" )

#colors<-distinctColorPalette(length(unique(sort(experiment$tissue))))
         
tsne2<-ggplot(info.norm, aes(x = tsne1, y = tsne2, colour =experiment_tsne$Tissue_new)) + 
  geom_point(alpha = 0.5,size=3) + theme_classic()+
  xlab("t-SNE-1")+ylab("t-SNE-2")+
  theme(panel.background=element_blank(),
        panel.grid.minor=element_blank(),panel.grid.major=element_blank(),
        legend.text = element_text(size=18),
        legend.title = element_text(size=1,face="bold"),
        legend.position='none',
        axis.title.y=element_text(size=18,color="black"),
        axis.title.x=element_text(size=18,color="black"),
        axis.text.y=element_text(size=16,color="black"),
        axis.text.x=element_text(size=16,color="black",hjust=0.95,vjust=0.2))+
  scale_colour_manual(name="Tissues",values=col)#+stat_ellipse(type="norm")
ggsave(tsne2,file='/Users/baizhonghao/Downloads/17tissues-tsne2.pdf',width=4,height=4,dpi=300)

####  expression median conservation
expression<-expression_inter[match(experiment_tsne$SampleID,rownames(expression_inter)),]
experiment_tsne$annotation<-paste(experiment_tsne$Tissue_new,experiment_tsne$Species,sep=" - ")
head(experiment_tsne)
length(table(experiment_tsne$annotation))
expression_hc<-expression[match(experiment_tsne$SampleID,rownames(expression)),]
median_tissue_tpm<-aggregate(expression_hc,list(as.factor(experiment_tsne$annotation)),median)
rownames(median_tissue_tpm)<-median_tissue_tpm[,1]
(median_tissue_tpm)[1:10,1:10]
corr_mat=cor(t(median_tissue_tpm[,-1]))
dim(corr_mat)
head(corr_mat[1:10,1:10])

#Annotation bar#

sample_info<-array(NA, dim=c(length(rownames(corr_mat)),2))
colnames(sample_info)<-c("Species","Tissues")
rownames(sample_info)<-rownames(corr_mat)
sample_info<-as.data.frame(sample_info)
split<-as.data.frame(strsplit(rownames(sample_info),split = " - "))
split<-as.data.frame(t(split))
sample_info$Species<-split$V2
sample_info$Tissues<-split$V1
ann_colors<-list(
  Species=c(Human="blue",Pig="red"),
  Tissues=c(Adipose="#CC66FF",Artery="#AAAAFF",Blood="#FF0000",Colon="#8EABD2",Frontal_cortex="#FFD700",Heart="#33CCCC",
            Hypothalamus="#FDFDBF",Ileum="#8EA9DB",Kidney="#8B0F55",Liver="#E2EFDA",Lung="#7570B3",Muscle="#AAEEFF",Ovary="#99BB88",
            Pituitary="#FFDD99",Spleen="#FF6600",Testis="#A6CEE3",Uterus="#A6761D",`Adrenal Gland`="#ed1299",Bladder="#09f9f5",
            Blastocyst="#246b93",Blastomere="#cc8e12",Cartilage="#d561dd",`Cervix Uteri`="#c93f00",Duodenum="#ddd53e",Embryo="#4aef7b",
            Esophagus="#e86502",`Fallopian Tube`="#9ed84e",Fetal_thymus="#39ba30",Jejunum="#6ad157",Lymph_node="#8249aa",Macrophage="#99db27",
            Milk="#e07233",Morula="#ff523f",Nerve="#ce2523",Oocyte="#f7aa5d",Pancreas="#cebb10",pEPSC="#03827f",Placenta="#931635",
            Prostate="#373bbf",`Salivary Gland`="#a1ce4c",Skin="#ef3bb6",Stomach="#d66551",Synovial_membrane="#1a918f",Thyroid="#ff66fc",Vagina="#2927c4",Breast="#7149af")
)
corr_mat[lower.tri(corr_mat)]=NA
#sample_info$Species<-factor(sample_info$Species,levels = c("Human","Pig"))
#sample_info$Tissues<-factor(sample_info$Tissues,levels = c("Adipose","Artery","Blood","Colon","Frontal_cortex","Heart","Hypothalamus","Ileum",
#                                                           "Kidney","Liver","Lung","Muscle","Ovary", "Pituitary","Spleen","Testis","Uterus","Adrenal Gland","Bladder","Blastocyst","Blastomere", "Cartilage",
#                         ,annotation_colors = ann_colors,                                  "Cervix Uteri","Duodenum","Embryo","Esophagus","Fallopian Tube","Fetal_thymus","Jejunum","Lymph_node","Macrophage",       
#                                                           "Milk","Morula","Nerve","Oocyte","Pancreas","pEPSC","Placenta","Prostate","Salivary Gland","Skin","Stomach","Synovial_membrane","Thyroid","Vagina" ))
library(pheatmap)
pheatmap(corr_mat,annotation_col = sample_info,annotation_colors = ann_colors,file="/Users/baizhonghao/Downloads/hc-median.pdf",width=20,height =20,cluster_rows = T,cluster_cols = T, show_rownames=F,show_colnames=F,res=300,fontsize =18)
pheatmap(corr_mat,annotation_col = sample_info,annotation_colors = ann_colors,file="/Users/baizhonghao/Downloads/hc-mediantri.pdf",width=20,height =20,cluster_rows = F,cluster_cols = F, show_rownames=F,show_colnames=F,res=300,fontsize =18,border_color = NA,na_col = NA)
expression<-expression_inter
save(experiment,expression,file="/Users/baizhonghao/Downloads/human-pig GTEx/DataMatrix20217.RData")
experiment<-experiment_tsne
expression<-expression_hc
