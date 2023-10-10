library(ComplexHeatmap)
library(circlize)

#####代谢物聚类分析#####
# 导入数据
data <- read.table("./data/m.txt", header = TRUE, row.names = 1)

# 计算相关矩阵
cor_matrix <- cor(data)

# 定义聚类方法
hc_row = hclust(dist(cor_matrix), method = "ward.D2")
hc_col = hclust(dist(t(cor_matrix)), method = "ward.D2")

# 选择适当的截断数目以确定簇的数量（在这个例子中我们设置为3）
clusters <- cutree(hc_row, k = 3)

# 查看簇的分配情况
print(clusters)

# 定义颜色映射
#color_mapping = colorRamp2(c(0, 1), c("white", "red"))
color_mapping = colorRampPalette(c("blue","white","red"))(100)

# 绘制热图
pdf("./fig/meta_map.pdf",width = 3.5,height = 2.5)
Heatmap(cor_matrix,
                  cluster_rows = hc_row,       # 将聚类结果应用于行
                  cluster_columns = hc_col,    # 将聚类结果应用于列
                  col = color_mapping,         # 应用颜色映射
                  name = "correlation",        # 图例标题
                  show_row_dend = TRUE,        # 显示行聚类树状图
                  show_column_dend = TRUE,     # 显示列聚类树状图
                  show_row_names = F,       # 显示行名称
                  show_column_names = F, # 显示列名称
        )
 
                  




# 绘制热图
dev.off()



#####转录组聚类分析#####
# 导入数据
data <- read.table("./data/g.txt", header = TRUE, row.names = 1)

# 计算相关矩阵
cor_matrix <- cor(data)

# 定义聚类方法
hc_row = hclust(dist(cor_matrix), method = "ward.D2")
hc_col = hclust(dist(t(cor_matrix)), method = "ward.D2")

# 选择适当的截断数目以确定簇的数量（在这个例子中我们设置为3）
clusters <- cutree(hc_row, k = 3)

# 查看簇的分配情况
print(clusters)


# 定义颜色映射
#color_mapping = colorRamp2(c(0,0.5, 1), c("blue","white", "red"))
color_mapping = colorRampPalette(c("blue","white","red"))(100)

# 绘制热图
pdf("./fig/rna_map.pdf",width = 3.5,height = 2.5)
Heatmap(cor_matrix,
                  cluster_rows = hc_row,       # 将聚类结果应用于行
                  cluster_columns = hc_col,    # 将聚类结果应用于列
                  col = color_mapping,         # 应用颜色映射
                  name = "correlation",        # 图例标题
                  show_row_dend = TRUE,        # 显示行聚类树状图
                  show_column_dend = TRUE,     # 显示列聚类树状图
                  show_row_names = F,       # 显示行名称
                  show_column_names = F)    # 显示列名称

# 绘制热图
dev.off()



######代谢组PCA #####
library(vegan)
library(ape)
library(ggplot2)
library(grid)
library(dplyr)
library(multcomp)
library(patchwork)
data <- read.table("./data/m.txt", header = TRUE, row.names = 1)
data <- t(data)
groups <- read.table("./data/mgroup.txt", header = TRUE)
colnames(groups) <- c("V1","V2")

length=length(unique(as.character(groups$V1)))
times1=length%/%8
res1=length%%8
times2=length%/%5
res2=length%%5
col1=rep(1:8,times1)
col=c(col1,1:res1)
pich1=rep(c(21:24),times2)
pich=c(pich1,15:(15+res2))
cbbPalette <- c("#B2182B","#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7","#CC6666","#9999CC","#66CC99","#99999",                "#ADD1E5")

# 对数据进行标准化
data_scaled <- scale(data)
# 计算主成分分析
pca_result <- prcomp(data_scaled)
pca_result <- prcomp(data)
PC1 = pca_result$x[,1]
PC2 = pca_result$x[,2]
plotdata <- data.frame(rownames(pca_result$x),PC1,PC2,groups$V2)
colnames(plotdata) <-c("sample","PC1","PC2","Group")

# 计算所占方差比例
explained_var_ratio <- pca_result$sdev^2 / sum(pca_result$sdev^2)

# 提取前两个主成分的贡献值
PC1_contribution <- explained_var_ratio[1]
PC2_contribution <- explained_var_ratio[2]
# 输出结果
cat("PC1 contribution:", round(PC1_contribution * 100, 2), "%\n")
cat("PC2 contribution:", round(PC2_contribution * 100, 2), "%\n")

pc1 <- round(PC1_contribution * 100, 2)
pc2 <- round(PC2_contribution * 100, 2)

plotdata$Group <- factor(plotdata$Group,levels = c("AC","AC+RS","H","H+RS"))


yf <- plotdata
yd1 <- yf %>% group_by(Group) %>% summarise(Max = max(PC1))
yd2 <- yf %>% group_by(Group) %>% summarise(Max = max(PC2))
yd1$Max <- yd1$Max + max(yd1$Max)*0.1
yd2$Max <- yd2$Max + max(yd2$Max)*0.1
fit1 <- aov(PC1~Group,data = plotdata)
tuk1<-glht(fit1,linfct=mcp(Group="Tukey"))
res1 <- cld(tuk1,alpah=0.05)
fit2 <- aov(PC2~Group,data = plotdata)
tuk2<-glht(fit2,linfct=mcp(Group="Tukey"))
res2 <- cld(tuk2,alpah=0.05)
test <- data.frame(PC1 = res1$mcletters$Letters,PC2 = res2$mcletters$Letters,
                   yd1 = yd1$Max,yd2 = yd2$Max,Group = yd1$Group)
test$Group <- factor(test$Group,levels = c("AC","AC+RS","H","H+RS"))


p1 <- ggplot(plotdata,aes(Group,PC1)) +
  geom_boxplot(aes(fill = Group)) +
  geom_text(data = test,aes(x = Group,y = yd1,label = PC1),
            size = 3,color = "black",fontface = "bold") +
  coord_flip() +
  scale_fill_manual(values=cbbPalette) +
  theme_bw()+
  theme(axis.ticks.length = unit(0.4,"lines"), 
        axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"), 
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_text(colour='black',size=10,face = "bold"),
        axis.text.x=element_blank(),
        legend.position = "none")
p1
p3 <- ggplot(plotdata,aes(Group,PC2)) +
  geom_boxplot(aes(fill = Group)) +
  geom_text(data = test,aes(x = Group,y = yd2,label = PC2),
            size = 3,color = "black",fontface = "bold") +
  scale_fill_manual(values=cbbPalette) +
  theme_bw()+
  theme(axis.ticks.length = unit(0.4,"lines"), 
        axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"), 
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_text(colour='black',size=10,angle = 45,
                                 vjust = 1,hjust = 1,face = "bold"),
        axis.text.y=element_blank(),
        legend.position = "none")
p3

p2<-ggplot(plotdata, aes(PC1, PC2)) +
  geom_point(aes(fill=Group),size=2,pch = 21)+
  scale_fill_manual(values=cbbPalette,name = "Group")+
  xlab(paste("PC1 ( ",pc1,"%"," )",sep="")) + 
  ylab(paste("PC2 ( ",pc2,"%"," )",sep=""))+
  xlim(ggplot_build(p1)$layout$panel_scales_y[[1]]$range$range) +
  ylim(ggplot_build(p3)$layout$panel_scales_y[[1]]$range$range) +
  theme(text=element_text(size=10))+
  geom_vline(aes(xintercept = 0),linetype="dotted")+
  geom_hline(aes(yintercept = 0),linetype="dotted")+
  theme(panel.background = element_rect(fill='white', colour='black'),
        panel.grid=element_blank(), 
        axis.title = element_text(color='black',size=10,face = "bold"),
        axis.ticks.length = unit(0.4,"lines"), axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"), 
        axis.title.x=element_text(colour='black',face = "bold", size=10,vjust = 7),
        axis.title.y=element_text(colour='black',face = "bold",size=10,vjust = -7),
        axis.text=element_text(colour='black',size=10,face = "bold"),
       # legend.title=element_text(size = 10,face = "bold"),
        legend.text=element_text(size=10,face = "bold"),
        legend.title = element_blank(),
        legend.key=element_blank(),legend.position = c(0.78,0.75),
        legend.background = element_rect(colour = "black"),
        legend.key.height=unit(0.5,"cm")) +
  guides(fill = guide_legend(ncol = 1))
p2

otu.adonis=adonis2(data~V2,data = groups,distance = "bray")

p4 <- ggplot(plotdata, aes(PC1, PC2)) +
  geom_text(aes(x = -0.5,y = 0.6,label = paste("PERMANOVA:\ndf = ",otu.adonis$Df[1],  "\nR2 = ",round(otu.adonis$R2[1],4),  "\np-value = ",otu.adonis$`Pr(>F)`[1],sep = "")),
            size = 2) +
  theme_bw() +
  xlab("") + ylab("") +
  theme(panel.grid=element_blank(), 
        axis.title = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank())
p4

p5 <- p1 + p4 + p2 + p3 + 
  plot_layout(heights = c(1,4),widths = c(4,1),ncol = 2,nrow = 2)
p5
ggsave("./fig/meta_pca.pdf",width = 17,height = 10,units = "cm")

######转录组PCA #####
library(vegan)
library(ape)
library(ggplot2)
library(grid)
library(dplyr)
library(multcomp)
library(patchwork)
data <- read.table("./data/g2.txt", header = TRUE, row.names = 1)
data <- t(data)
groups <- read.table("./data/ggroup.txt", header = TRUE)
colnames(groups) <- c("V1","V2")

length=length(unique(as.character(groups$V1)))
times1=length%/%8
res1=length%%8
times2=length%/%5
res2=length%%5
col1=rep(1:8,times1)
col=c(col1,1:res1)
pich1=rep(c(21:24),times2)
pich=c(pich1,15:(15+res2))
cbbPalette <- c("#B2182B","#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7","#CC6666","#9999CC","#66CC99","#99999",                "#ADD1E5")

# 对数据进行标准化
data_scaled <- scale(data)
# 计算主成分分析
pca_result <- prcomp(data_scaled)
pca_result <- prcomp(data)
PC1 = pca_result$x[,1]
PC2 = pca_result$x[,2]
plotdata <- data.frame(rownames(pca_result$x),PC1,PC2,groups$V2)
colnames(plotdata) <-c("sample","PC1","PC2","Group")

# 计算所占方差比例
explained_var_ratio <- pca_result$sdev^2 / sum(pca_result$sdev^2)

# 提取前两个主成分的贡献值
PC1_contribution <- explained_var_ratio[1]
PC2_contribution <- explained_var_ratio[2]
# 输出结果
cat("PC1 contribution:", round(PC1_contribution * 100, 2), "%\n")
cat("PC2 contribution:", round(PC2_contribution * 100, 2), "%\n")

pc1 <- round(PC1_contribution * 100, 2)
pc2 <- round(PC2_contribution * 100, 2)

plotdata$Group <- factor(plotdata$Group,levels = c("AC","AC+RS","H","H+RS"))


yf <- plotdata
yd1 <- yf %>% group_by(Group) %>% summarise(Max = max(PC1))
yd2 <- yf %>% group_by(Group) %>% summarise(Max = max(PC2))
yd1$Max <- yd1$Max + max(yd1$Max)*0.1
yd2$Max <- yd2$Max + max(yd2$Max)*0.1
fit1 <- aov(PC1~Group,data = plotdata)
tuk1<-glht(fit1,linfct=mcp(Group="Tukey"))
res1 <- cld(tuk1,alpah=0.05)
fit2 <- aov(PC2~Group,data = plotdata)
tuk2<-glht(fit2,linfct=mcp(Group="Tukey"))
res2 <- cld(tuk2,alpah=0.05)
test <- data.frame(PC1 = res1$mcletters$Letters,PC2 = res2$mcletters$Letters,
                   yd1 = yd1$Max,yd2 = yd2$Max,Group = yd1$Group)
test$Group <- factor(test$Group,levels = c("AC","AC+RS","H","H+RS"))


p1 <- ggplot(plotdata,aes(Group,PC1)) +
  geom_boxplot(aes(fill = Group)) +
  geom_text(data = test,aes(x = Group,y = yd1,label = PC1),
            size = 3,color = "black",fontface = "bold") +
  coord_flip() +
  scale_fill_manual(values=cbbPalette) +
  theme_bw()+
  theme(axis.ticks.length = unit(0.4,"lines"), 
        axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"), 
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_text(colour='black',size=10,face = "bold"),
        axis.text.x=element_blank(),
        legend.position = "none")
p1
p3 <- ggplot(plotdata,aes(Group,PC2)) +
  geom_boxplot(aes(fill = Group)) +
  geom_text(data = test,aes(x = Group,y = yd2,label = PC2),
            size = 3,color = "black",fontface = "bold") +
  scale_fill_manual(values=cbbPalette) +
  theme_bw()+
  theme(axis.ticks.length = unit(0.4,"lines"), 
        axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"), 
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_text(colour='black',size=10,angle = 45,
                                 vjust = 1,hjust = 1,face = "bold"),
        axis.text.y=element_blank(),
        legend.position = "none")
p3

p2<-ggplot(plotdata, aes(PC1, PC2)) +
  geom_point(aes(fill=Group),size=2,pch = 21)+
  scale_fill_manual(values=cbbPalette,name = "Group")+
  xlab(paste("PC1 ( ",pc1,"%"," )",sep="")) + 
  ylab(paste("PC2 ( ",pc2,"%"," )",sep=""))+
  xlim(ggplot_build(p1)$layout$panel_scales_y[[1]]$range$range) +
  ylim(ggplot_build(p3)$layout$panel_scales_y[[1]]$range$range) +
  theme(text=element_text(size=10))+
  geom_vline(aes(xintercept = 0),linetype="dotted")+
  geom_hline(aes(yintercept = 0),linetype="dotted")+
  theme(panel.background = element_rect(fill='white', colour='black'),
        panel.grid=element_blank(), 
        axis.title = element_text(color='black',size=10,face = "bold"),
        axis.ticks.length = unit(0.4,"lines"), axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"), 
        axis.title.x=element_text(colour='black',face = "bold", size=10,vjust = 7),
        axis.title.y=element_text(colour='black',face = "bold",size=10,vjust = -2),
        axis.text=element_text(colour='black',size=10,face = "bold"),
        #legend.title=element_text(size = 10,face = "bold"),
        legend.text=element_text(size=10,face = "bold"),
        legend.title = element_blank(),
        legend.key=element_blank(),legend.position = c(0.78,0.75),
        legend.background = element_rect(colour = "black"),
        legend.key.height=unit(0.5,"cm")) +
  guides(fill = guide_legend(ncol = 1))
p2

otu.adonis=adonis2(data~V2,data = groups,distance = "bray")

p4 <- ggplot(plotdata, aes(PC1, PC2)) +
  geom_text(aes(x = -0.5,y = 0.6,label = paste("PERMANOVA:\ndf = ",otu.adonis$Df[1],  "\nR2 = ",round(otu.adonis$R2[1],4),  "\np-value = ",otu.adonis$`Pr(>F)`[1],sep = "")),
            size = 2) +
  theme_bw() +
  xlab("") + ylab("") +
  theme(panel.grid=element_blank(), 
        axis.title = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank())
p4

p5 <- p1 + p4 + p2 + p3 + 
  plot_layout(heights = c(1,4),widths = c(4,1),ncol = 2,nrow = 2)
p5
ggsave("./fig/rna_pca.pdf",width = 17,height = 10,units = "cm")
