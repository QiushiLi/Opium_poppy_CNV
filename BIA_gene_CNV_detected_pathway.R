setwd("~/Desktop/CNV_heatmap")
library(devtools)
#install_github("jokergoo/ComplexHeatmap")
library(ComplexHeatmap)
options(stringsAsFactors = F)
library(pheatmap)
library(RColorBrewer)
library(vegan)
spe<-read.delim(file="BIA_gene_CNV_detected_pathway.txt",header=T,check.names=FALSE,sep="\t",fill=T,row.names = 1)
spe2<-spe

#dat.bray<-vegdist(spe2,method="bray",binary=FALSE,diag=FALSE,upper=FALSE,na.rm=FALSE)
dat.bray<-vegdist(spe2,method="bray")
write.table(spe2,"BIA_gene_CNV_detected_pathway.xls",sep="\t",eol="\n",quote=FALSE,col.names = T,row.names =T)
q1=0
q2=1
q3=2
q4=3
q5=15
library(RColorBrewer)
#display.brewer.all()
myColor<-c("#FFFFFF" ,"#cccccc", "#808080", colorRampPalette(c("#006600", "#ccff33","#ff6600") )(10))
brks<- c(-1,0,1,2,c(seq(3, 15, length.out=13)))

rowann<-read.table(file="rlg.txt",header=T,check.names=FALSE,sep="\t",fill=T,row.names=1)
rowann$cScaf<-factor(rowann$cScaf,levels = rowann$cScaf[!duplicated(rowann$cScaf)])
ann_color<-rep(brewer.pal(n=3,name = 'Set1')[1:2],5)
names(ann_color)<-rowann$cScaf[!duplicated(rowann$cScaf)]
anno_colors <- list( cScaf= ann_color)

pheatmap(spe2,filename = 'heatmap_output.pdf',
         show_colnames = T,breaks = brks,color=myColor,scale = "none", 
         silent=FALSE,fontsize=10,cellwidth = 36,cellheight = 20,
         treeheight_row=70,
         cluster_cols =F,cluster_rows =F,show_rownames = T,
         annotation_row=rowann,annotation_colors = anno_colors,annotation_names_row = TRUE,
         annotation_legend = F)

library(ComplexHeatmap)
df<-read.delim('pathway.txt',header = T)
mycol<-read.delim('color.txt',header = F)
mycol$V1<-paste("#",mycol$V1,sep = '')
pdf('legend.pdf',width = 12,height = 12,onefile = F)
gd = Legend(labels = df$Pathway, title = "Pathway", type = "lines",
            legend_gp = gpar(col =mycol$V1, lty = 1,lwd=2), grid_width = unit(1, "cm"),background = NA)

lgd_list_vertical = packLegend(gd )
grid.newpage()
grid.draw(lgd_list_vertical)
dev.off()
