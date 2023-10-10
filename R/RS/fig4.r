#set color

colors1 <- c("#8dd3c7", "#ffffb3", "#bebada", "#fb8072", "#80b1d3", "#fdb462", "#b3de69", "#fccde5", "#d9d9d9", "#bc80bd", "#ccebc5", "#ffed6f", "#1f78b4")

library(RColorBrewer)

brewer_colors <- brewer.pal(9, "Set3")
brewer_colors <- append(brewer_colors, brewer.pal(4, "Set1")[1:4])


## HDMB 饼图
## 主题
theme_zg <- function(..., bg='white'){
  require(grid)
  theme_classic(...) +
    theme(rect=element_rect(fill=bg),
          plot.margin=unit(rep(0.5,4), 'lines'),
          panel.background=element_rect(fill='transparent', color='black'),
          panel.border=element_rect(fill='transparent', color='transparent'),
          panel.grid=element_blank(),
          axis.title.x =element_text(size=10),
          axis.title.y=element_text(size=10),
          axis.title = element_text(color='black',vjust=0.1),
          legend.title=element_blank(),
          legend.key=element_rect(fill='transparent', color='transparent'))}
library(ggsci)
library(ggplot2)
pie <- read.table("./data/HDMB.txt", header=T, sep="\t",dec=".")
p <- ggplot(data = pie, mapping = aes(x = 'Content', y = Number,fill = Superclass)) + 
  geom_bar(stat = 'identity', position = 'stack') +
  scale_fill_manual(values =colors1 ) +
  theme_zg()
p

p + coord_polar(theta = 'y') + labs(x = '', y = '', title = '') + 
  theme(axis.text = element_blank()
  )
ggsave("./fig/hdmb.pdf",width = 15,height = 8,units = "cm")



 
##### 火山图

  rm(list = ls())
  theme_zg <- function(..., bg='white'){
    require(grid)
    theme_classic(...) +
      theme(rect=element_rect(fill=bg),
            plot.margin=unit(rep(0.5,4), 'lines'),
            panel.background=element_rect(fill='transparent', color='black'),
            panel.border=element_rect(fill='transparent', color='transparent'),
            panel.grid=element_blank(),
            axis.title.x =element_text(size=14),
            axis.title.y=element_text(size=14),
            axis.title = element_text(color='black',vjust=0.1),
            legend.title=element_blank(),
            legend.key=element_rect(fill='transparent', color='transparent'))}
  
  library(ggplot2)
  library(ggrepel)
  
  df <- read.delim("./data/vol-h.txt",header = T, sep = '\t',row.names = 1, stringsAsFactors = FALSE, check.names = FALSE)
  head(df)
  fcthreshold = 1
  pthreshold =0.05
  
  p.vol.h <- ggplot(data = df,
                    aes(x = FC,y = -log10(P_value),colour = Regulate,fill = Regulate))+
    scale_color_manual(values = c('green','grey','red'))+
    geom_point(alpha = 0.5,size = 1.0)+
    # 辅助线
    geom_vline(xintercept = c(-log2(fcthreshold),log2(fcthreshold)),lty = 4,col = "black",lwd = 0.8)+
    geom_hline(yintercept = -log10(pthreshold),lty = 4,col = "black",lwd = 0.8)+
    theme_zg()+
    labs(x = "log2(FoldChange)",y = "-log10(p-adjust)",title = "H+RS vs H")+
    # 坐标轴标题、标签和图例相关设置
    theme(axis.text = element_text(size = 11),axis.title = element_text(size = 10), # 坐标轴标签和标题
          plot.title = element_text(hjust = 0.5,size = 14,face = "bold"), # 标题
          legend.text = element_text(size = 11),legend.title = element_text(size = 10), # 图例标签和标题
          plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm")) # 图边距
  p.vol.h
  ggsave(p.vol.h,filename = "./fig/Volcano1.pdf",width = 10,height = 6,units = "cm")
  
  df <- read.delim("./data/vol-h.txt",header = T, sep = '\t',row.names = 1, stringsAsFactors = FALSE, check.names = FALSE)
  head(df)
  fcthreshold = 1
  pthreshold =0.05
  
  p.vol.ac <- ggplot(data = df,
                     aes(x = FC,y = -log10(P_value),colour = Regulate,fill = Regulate))+
    scale_color_manual(values = c('green','grey','red'))+
    geom_point(alpha = 0.5,size = 1.0)+
    # 辅助线
    geom_vline(xintercept = c(-log2(fcthreshold),log2(fcthreshold)),lty = 4,col = "black",lwd = 0.8)+
    geom_hline(yintercept = -log10(pthreshold),lty = 4,col = "black",lwd = 0.8)+
    theme_zg()+
    labs(x = "log2(FoldChange)",y = "-log10(p-adjust)",title = "AC+RS vs AC")+
    # 坐标轴标题、标签和图例相关设置
    theme(axis.text = element_text(size = 11),axis.title = element_text(size = 10), # 坐标轴标签和标题
          plot.title = element_text(hjust = 0.5,size = 14,face = "bold"), # 标题
          legend.text = element_text(size = 11),legend.title = element_text(size = 10), # 图例标签和标题
          plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm")) # 图边距
  p.vol.ac
  library(ggpubr)
  vol <- ggarrange(p.vol.h,p.vol.ac,ncol = 2,nrow = 1,common.legend = T,legend = "right")
  vol
  ggsave(vol,filename = "./fig/Volcano.pdf",width = 17,height = 6,units = "cm")

  
#### 代谢物热图
  abundance <- read.table("./data/heatmap.txt",header = TRUE,sep = "\t",quote = "")
  
  abundance <- as.data.frame(abundance)
  taxon <- gsub(".*(o__)","",abundance$Taxon)
  taxon <- gsub(".*(;)","",taxon)
  taxon <- gsub(".*(Other)","",taxon)
  taxon <- gsub("(\\[)","",taxon)
  taxon <- gsub("(\\])","",taxon)
  abundance$Taxon <- taxon
  abundance[abundance==""] <- NA
  abundance <- na.omit(abundance)
  abundance <- subset(abundance,Taxon != "Unclassified")
  abundance <- subset(abundance,Taxon != "norank")
  abundance <- subset(abundance,Taxon != "uncultured")
  abundance <- t(abundance)
  colnames(abundance) <- abundance[1,]
  abundance <- abundance[-1,]
  abundance <- as.matrix(abundance)
  f.abundance <- matrix(as.numeric(abundance),nrow = nrow(abundance))
  rownames(f.abundance) <- rownames(abundance)
  colnames(f.abundance) <- colnames(abundance)
  f.abundance1 <- t(f.abundance)
  sum <- apply(f.abundance1,1,sum) 
  f.abundance2 <- cbind(f.abundance1,sum)
  f.abundance3 <- as.data.frame(f.abundance2)
  f.abundance4 <- f.abundance3[order(f.abundance3[,"sum"],decreasing = T),]
  f.abundance5 <- subset(f.abundance4, select = -sum)
  f.abundance6 <- f.abundance5[1:25,]
  f.abundance7 <- t(f.abundance6)
  sample.name <- read.table("./data/map-group.txt",head=F,colClasses=c("character","character"),sep = "\t")
  sample.name <- sample.name[-1,]
  f.abundance8 <- as.data.frame(f.abundance7)
  f.abundance8$V1 <- rownames(f.abundance8)
  f.abundance9 <- merge(sample.name,f.abundance8,by = "V1",sort = F)
  rownames(f.abundance9) <- f.abundance9$V2
  f.abundance10 <- subset(f.abundance9, select = -V1)
  f.abundance11 <- f.abundance10[,-c(1,2)]
  f.abundance12 <- t(f.abundance11)

  library(pheatmap)

  
  
  library(dplyr)
  group<-unique(sample.name$V3)
  Gnum<-length(group)
  number <- sample.name %>% group_by(V3) %>% count(V3)
  breaks <- c()
  k <- 1
  for(i in 1:(Gnum-1)){
    breaks[k] <- sum(number[1:i,2])
    k <- k+1
  }
  
  annotation_col <- data.frame(Group=factor(rep(sample.name$V3)))
  rownames(annotation_col)
  head(annotation_col)
  rownames(annotation_col) <- colnames(f.abundance12)
  
  
  
  #bk <- c(seq(-6,-0.1,by=0.01),seq(0,8,by=0.01))
  #log2转换
  scale_test <- apply(f.abundance12,2,function(x)log2(x+1))
  #scale转换
  # scale_test <- apply(f.abundance12,2,scale)
  
  
  pdf("./fig/map.pdf",width = 7.2,height = 5)
  
  pheatmap(mat=scale_test,fontsize = 6,
           angle_row = 45,
           cellwidth = 6,#格子的宽度
           scale = "row",
           show_colnames = F,###是否显示列名
           show_rownames = TRUE,
           annotation_names_col = T,
           fontsize_row  = 9, # 字体大小
           cluster_rows = F,
          # cluster_cols = F,
           annotation_legend =T, ##显示图例
           annotation_col = annotation_col,##列的分组信息，如果想加颜色可以用annotation_colors
           gaps_col = breaks,
           color = colorRampPalette(c("blue","white","red"))(100))
  
dev.off()


##### KEGG富集分析
library (dplyr) 
library (ggplot2)  
library(tidyverse)
library(openxlsx)
#数值导入#
H = read.xlsx('./data/KEGG.xlsx',sheet= "Sheet1",sep=',')
head(H)
H <- H[1:20,]

KEGG <-  ggplot(H,aes(y=Description,x=Ratio_in_study))+
  geom_point(aes(size=Num,color=Padjust))+
  theme_zg()+
  scale_color_gradient(low = "red",high ="blue")+
  #  facet_wrap(~Term , scales="free" ,nrow = 3) +
  theme(plot.title = element_text(hjust=0.5)) +
  theme(axis.text.y = element_text(color = "black",size = 9.0)) + 
  labs(color=expression(Padjust,size="Count"), 
       x="Gene Ratio",y=" ",title="KEGG Enrichment")

KEGG
ggsave("./fig/KEGG.pdf",width = 15,height = 9,units = "cm")

