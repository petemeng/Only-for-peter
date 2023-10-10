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

##### 火山图
library(ggplot2)
library(ggrepel)

df <- read.delim("./data/vol-h.txt",header = T, sep = '\t',row.names = 1, stringsAsFactors = FALSE, check.names = FALSE)
head(df)
fcthreshold = 1
pthreshold =0.05

p.vol.h <- ggplot(data = df,
                  aes(x = Log2FC,y = -log10(Padjust),colour = Regulate,fill = Regulate))+
  scale_color_manual(values = c('green','grey','red'))+
  geom_point(alpha = 0.5,size = 1.0)+
  # 辅助线
  geom_vline(xintercept = c(-log2(fcthreshold),log2(fcthreshold)),lty = 4,col = "black",lwd = 0.8)+
  geom_hline(yintercept = -log10(pthreshold),lty = 4,col = "black",lwd = 0.8)+
  theme_zg()+
  labs(x = "log2(FoldChange)",y = "-log10(p-adjust)",title = "H+RS vs H")+
  # 坐标轴标题、标签和图例相关设置
  theme(axis.text = element_text(size = 8),axis.title = element_text(size = 7), # 坐标轴标签和标题
        plot.title = element_text(hjust = 0.5,size = 8,face = "bold"), # 标题
        legend.text = element_text(size = 8),legend.title = element_text(size = 7),# 图例标签和标题
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm")) # 图边距
p.vol.h
ggsave(p.vol.h,filename = "./fig/Volcano1.pdf",width = 7,height = 5,units = "cm")

df <- read.delim("./data/vol-ac.txt",header = T, sep = '\t',row.names = 1, stringsAsFactors = FALSE, check.names = FALSE)
head(df)
fcthreshold = 1
pthreshold =0.05

p.vol.ac <- ggplot(data = df,
                   aes(x = Log2FC,y = -log10(Padjust),colour = Regulate,fill = Regulate))+
  scale_color_manual(values = c('green','grey','red'))+
  geom_point(alpha = 0.5,size = 1.0)+
  # 辅助线
  geom_vline(xintercept = c(-log2(fcthreshold),log2(fcthreshold)),lty = 4,col = "black",lwd = 0.8)+
  geom_hline(yintercept = -log10(pthreshold),lty = 4,col = "black",lwd = 0.8)+
  theme_zg()+
  labs(x = "log2(FoldChange)",y = "-log10(p-adjust)",title = "AC+RS vs AC")+
  # 坐标轴标题、标签和图例相关设置
  theme(axis.text = element_text(size = 8),axis.title = element_text(size = 7), # 坐标轴标签和标题
        plot.title = element_text(hjust = 0.5,size = 8,face = "bold"), # 标题
        legend.text = element_text(size = 8),legend.title = element_text(size = 7),# 图例标签和标题
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm")) # 图边距
p.vol.ac

library(ggpubr)
vol <- ggarrange(p.vol.h,p.vol.ac,ncol = 2,nrow = 1,common.legend = T,legend = "right")
vol
ggsave(vol,filename = "./fig/Volcano.pdf",width = 16,height = 6,units = "cm")



#### GO富集
library (dplyr) 
library (ggplot2)  
library(tidyverse)
library(openxlsx)
#数值导入#
BP = read.xlsx('./data/BP-H.xlsx',sheet= "Sheet1",sep=',')

head(BP)

display_number = c(10, 10, 10)  ##这三个数字分别代表选取的BP、CC、MF的数量
go_result_BP = as.data.frame(BP)[1:display_number[1], ]

ego = go_result_BP
ego$order=factor(rev(as.integer(rownames(ego))),labels = rev(ego$Description))
head(ego)
GO_H <-  ggplot(ego,aes(y=order,x=Ratio_in_study))+
  geom_point(aes(size=Number,color=Padjust))+
  scale_size_continuous(range = c(1, 3))+
  scale_y_discrete(labels=function(x) str_wrap(x, width=30)) +
  theme_bw()+
  scale_color_gradient(low = "red",high ="blue")+
  #facet_wrap(~Term , scales="free" ,nrow = 3) +
  theme(plot.title = element_text(hjust=0.5,size = 8,face = "bold"),
        axis.text.x = element_text(color = "black",size = 7.0,face = "bold"),
        axis.text.y = element_text(color = "black",size = 8.0,face = "bold"),
        legend.title = element_text(size = 6),
        axis.title.x = element_text(color = "black",size = 8.0,face = "bold"),
        axis.title.y = element_text(color = "black",size = 8.0,face = "bold"),
        legend.text = element_text(size = 6),
       # legend.title = element_text(size = 6),
        legend.key.size = unit(0.2, "cm")
        ) + 
  labs(color=expression(Padjust,size="Count"), 
       x="Gene Ratio",y="GO term",title="GO Enrichment(H+RS)")

GO_H

#数值导入#
BP = read.xlsx('./data/BP-AC.xlsx',sheet= "Sheet1",sep=',')
head(BP)

display_number = c(10, 10, 10)  ##这三个数字分别代表选取的BP、CC、MF的数量
go_result_BP = as.data.frame(BP)[1:display_number[1], ]
#go_result_CC = as.data.frame(CC)[1:display_number[2], ]
#go_result_MF = as.data.frame(MF)[1:display_number[3], ]

#ego = rbind(go_result_BP,go_result_CC,go_result_MF) 
ego =go_result_BP
ego$order=factor(rev(as.integer(rownames(ego))),labels = rev(ego$Description))
head(ego)
GO_AC <- ggplot(ego, aes(y = order, x = Ratio_in_study)) +
  geom_point(aes(size = Number, color = Padjust)) +
  scale_size_continuous(range = c(1, 3))+
  scale_y_discrete(labels=function(x) str_wrap(x, width=30)) +
  theme_bw() +
  scale_color_gradient(low = "red", high = "blue") +
#  facet_wrap(~Term, scales = "free", nrow = 3) +
  theme(plot.title = element_text(hjust=0.5,size = 8,face = "bold"),
        axis.text.x = element_text(color = "black",size = 7.0,face = "bold"),
        axis.text.y = element_text(color = "black",size = 8.0,face = "bold"),
        legend.title = element_text(size = 6),
        axis.title.x = element_text(color = "black",size = 8.0,face = "bold"),
        axis.title.y = element_text(color = "black",size = 8.0,face = "bold"),
        legend.text = element_text(size = 6),
    #    legend.title = element_text(size = 6),
        legend.key.size = unit(0.2, "cm")
  ) + 
  labs(color = expression(Padjust, size = "Count"),
       x = "Gene Ratio", y = "GO term", title = "GO Enrichment(AC+RS)")


GO_AC
library(ggpubr)
GO <- ggarrange(GO_H,GO_AC,ncol = 2,nrow = 1,common.legend = T,legend = "top")
GO
ggsave(filename = "./fig/GOenrich.pdf",width = 17.2,height = 8,units = "cm")




#### KEGG
library (dplyr) 
library (ggplot2)  
library(tidyverse)
library(openxlsx)
#数值导入#
H = read.xlsx('./data/KEGG-H.xlsx',sheet= "Sheet1",sep=',')
head(H)
H <- H[1:15,]

GO_H <-  ggplot(H,aes(y=Description,x=Ratio_in_study))+
  geom_point(aes(size=Num,color=Padjust))+
  scale_y_discrete(labels=function(x) str_wrap(x, width=30)) +
  scale_size_continuous(range = c(1, 5))+
  theme_zg()+
  scale_color_gradient(low = "red",high ="blue")+
  labs(color=expression(Padjust,size="Count"), 
       x="Gene Ratio",y="",title="KEGG Enrichment(H+RS)")+
  #  facet_wrap(~Term , scales="free" ,nrow = 3) +
  theme(plot.title = element_text(hjust=0.5,size = 8,face = "bold"),
        axis.text.x = element_text(color = "black",size = 7.0,face = "bold"),
        axis.text.y = element_text(color = "black",size = 7.0,face = "bold"),
        axis.title.x = element_text(color = "black",size = 8.0,face = "bold"),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.2, "cm"))  
 

GO_H
ggsave(filename = "./fig/KEGGenrich.pdf",width = 10,height = 10,units = "cm")



AC = read.xlsx('./data/KEGG-AC.xlsx',sheet= "Sheet1",sep=',')
head(AC)
AC <- AC[1:20,]

GO_AC <-  ggplot(AC,aes(y=Description,x=Ratio_in_study))+
  geom_point(aes(size=Num,color=Padjust))+
  theme_zg()+
  scale_color_gradient(low = "red",high ="blue")+
  #  facet_wrap(~Term , scales="free" ,nrow = 3) +
  theme(plot.title = element_text(hjust=0.5)) +
  theme(axis.text.y = element_text(color = "black",size = 14.0)) + 
  labs(color=expression(Padjust,size="Count"), 
       x="Gene Ratio",y="GO term",title="KEGG Enrichment(AC+RS)")

GO_AC
GO <- ggarrange(GO_H,GO_AC,labels = c("A","B"),ncol = 1,nrow = 2)
GO
ggsave(filename = "./fig/KEGGenrich.pdf",width = 12,height = 14)


##### TF分析
library(ggsci)
library(ggplot2)
tf <- read.table("./data/TF.txt",header = T)
ggplot(tf, aes(TFfamily, Gene_number, fill = TFfamily)) +
  geom_bar(stat = "identity",width = 0.5) +
  theme_zg() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        axis.title.x = element_text(color = "black",size = 8.0,face = "bold"), 
        axis.title.y = element_text(color = "black",size = 8.0,face = "bold"),
        axis.text.y = element_text(size = 7),
        legend.position = "none") +
  scale_fill_npg()


ggsave("./fig/TF.pdf",width = 8,height = 8,units = "cm")

