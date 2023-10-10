## alpha多样性
library(tidyverse)
library(vegan)
library(agricolae)
library(ggplot2)
library(LorMe)
library(ggpubr)
library(ape)
## 主题
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
## 自定义函数
alpha <- function(x, tree = NULL, base = exp(1)) {
  est <- estimateR(x)
  Richness <- est[1, ]
  Chao1 <- est[2, ]
  ACE <- est[4, ]
  Shannon <- diversity(x, index = 'shannon', base = base)
  Simpson <- diversity(x, index = 'simpson')    #Gini-Simpson 指数
  Pielou <- Shannon / log(Richness, base)
  goods_coverage <- 1 - rowSums(x == 1) / rowSums(x)
  
  result <- data.frame(Richness, Shannon, Simpson, Pielou, Chao1, ACE, goods_coverage)
  if (!is.null(tree)) {
    PD_whole_tree <- pd(x, tree, include.root = FALSE)[1]
    names(PD_whole_tree) <- 'PD_whole_tree'
    result <- cbind(result, PD_whole_tree)
  }
  result
}

otu <- read.delim('./data/selected_data.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
groups <- read.table("./data/selected_group2.txt", header=T, sep="\t",dec=".")
alpha_all <- alpha(otu, base = 2)
group <- arrange(groups,sample)
samplealpha<-cbind(group,alpha_all)

myaov <- aov(Richness~group,data = samplealpha)
aov_result <- LSD.test(myaov,"group",p.adj = "none")
stat = aov_result$groups
stat

richness <- ggplot(samplealpha,aes(group,Richness))+
  geom_boxplot(aes(group=group,color=factor(group)))+
  geom_point(size=2,aes(color=factor(group))) +
  theme_zg()+
  theme(plot.title = element_text( size = 15, hjust = 0.5))+
  xlab("")
richness

ggsave("./fig/Richness.pdf",width = 5,height = 3)


###PCOA
otu <- read.delim('./data/selected_data.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
groups <- read.table("./data/selected_group2.txt", header=T, sep="\t",dec=".")
groups <- arrange(groups,sample)

pcoa_otu <- otu

## PCOA分析
pcoa_otu <- vegdist(pcoa_otu,method = "bray")
pcoa <- pcoa(pcoa_otu,correction = "none",rn = NULL)
PC1 = pcoa$vectors[,1]
PC2 = pcoa$vectors[,2]
pc1 <-floor(pcoa$values$Relative_eig[1]*100)
pc2 <-floor(pcoa$values$Relative_eig[2]*100)
plotdata <- data.frame(rownames(pcoa$vectors),PC1,PC2,groups$group)
colnames(plotdata) <- c("sample","PC1","PC2","group")

## adonis分析
adonis_result <- adonis2(pcoa_otu~group,groups,permutations = 999)
adonis_result

## 作图
PCOA <- ggplot(plotdata,aes(PC1,PC2)) +
  geom_point(aes(fill=group),size=5,pch=21)+
  xlab(paste("PCO1 ( ",pc1,"%"," )",sep="")) + 
  ylab(paste("PCO2 ( ",pc2,"%"," )",sep="")) +
  theme(text=element_text(size=14))+
  geom_vline(aes(xintercept = 0),linetype="dotted")+
  geom_hline(aes(yintercept = 0),linetype="dotted")+
  #  geom_text(aes(label=sample),size=4,vjust=-1)+
  stat_ellipse(aes(color = group),level = 0.95,geom = "polygon",linetype=2,show.legend = FALSE,alpha = 0.1)+
  theme(panel.background = element_rect(fill='white', colour='black'),
        panel.grid=element_blank(), 
        axis.title = element_text(color='black',size=14),
        axis.ticks.length = unit(0.4,"lines"), axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"), 
        axis.title.x=element_text(colour='black', size=14,vjust = 2),
        axis.title.y=element_text(colour='black', size=14,vjust = -1),
        axis.text=element_text(colour='black',size=14),
        legend.title=element_blank(),
        legend.text=element_text(size=12),
        legend.key=element_blank(),legend.position = c(0.87,0.20),
        legend.key.height=unit(0.2,"cm"))
PCOA
ggsave("./fig/PCOA.pdf",width = 5,height = 3)



## 门水平物种组成分析
Palette <- c("#B2182B","#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7","#CC6666","#9999CC","#66CC99","#999999","#ADD1E5")

otu <- read.delim('./data/selected_data.txt', sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
otu <- data.frame(t(otu))
otu2 <- data.frame(ASV = row.names(otu), otu)

taxonmy <- read.table('./data/taxonomy.tsv',header = T,sep = "\t")
taxonmy <- taxonmy %>% 
  separate('Taxon',into = c("kingdom","phylum","class","order","fami1y","genus"),
           sep = ";",convert = T)
colnames(taxonmy)[1] <- "ASV"
colnames(otu)[1] <- "ASV"

abundance <- left_join(otu2,taxonmy,by="ASV")
rownames(abundance) <- abundance[,1]
abundance <- abundance[,-1]
phylum.sum <- aggregate(abundance[1:32],by=list(abundance$phylum),FUN = sum)
abundance.1 <- phylum.sum
sample=read.table("./data/selected_group2.txt",header = TRUE,sep = "\t")
sample <- arrange(sample,sample)
abundance.2=as.data.frame(t(abundance.1[,-1]))
colnames(abundance.2)=abundance.1[,1]
abundance.3=cbind(sample,abundance.2)
abundance.4 = abundance.3 
#abundance.4 <- abundance.4[!rownames(abundance.4) %in% c("B39","B42","B45","B47") , ]
abundance.5=aggregate(abundance.4,by=list(abundance.4$group),mean)
abundance.5=abundance.5[,-c(2:5)]
colnames(abundance.5)[1]=c("group")
abundance.6=t(abundance.5[,-c(1)])
colnames(abundance.6)<-abundance.5[,1]
abundance.1=as.data.frame(abundance.6)

format_percent <- function(input){
  out = input
  for (x in 1:length(input)){
    out[x] = input[x]/sum(input)
  }
  return(out)
}

abundance1=as.data.frame(apply(abundance.1,1,as.numeric))
abundance2=as.data.frame(apply(abundance1,1,format_percent))
colSums(abundance2)
abundance=cbind(rownames(abundance.1),abundance2)
colnames(abundance)[1] <- "Taxon"

taxon <- gsub(".*(p__)","",abundance$Taxon)

abundance$Taxon <- taxon

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

f.abundance9 <- f.abundance5[1:10,]
sum <- apply(f.abundance9,2,sum) 
Others <- 1-sum
f.abundance10 <- rbind(f.abundance9,Others)
rownames(f.abundance10)[11]<-"Others"
sample.name <- abundance.5[,1]
f.abundance11 <- as.data.frame(t(f.abundance10))
f.abundance12 <-cbind(sample.name,f.abundance11)

library(reshape2)
taxon <- melt(f.abundance12,id.vars = c("sample.name"),variable.name = "variable",value.name = "value")
colnames(taxon) <- c("group","Taxon","value")

library(ggplot2)
library(ggalluvial)
library(ggThemeAssist)
taxon$value=as.numeric(taxon$value)

phylum <- ggplot(data = taxon,aes(x = group, y = value, alluvium = Taxon, stratum = Taxon))+
  #geom_alluvium(aes(fill = Taxon),alpha = .5,width = 0.6) +  ##桑基图
  geom_stratum(aes(fill = Taxon),width = 0.6)+
  #分面
  #facet_wrap(Class~.,scales = "free_x",ncol = 2)+
  ylab(label = "Relative Abundance") + 
#  xlab(label = "Group")+
  xlab("")+
  scale_fill_manual(values = Palette)+
  theme_bw()+ theme(panel.grid=element_blank()) + 
  theme(panel.border = element_blank()) +
  theme(panel.background=element_rect(fill='transparent', 
                                      color='black'),plot.margin = unit(c(3,5,1,1),"mm"))+
  theme(axis.text.x=element_text(colour="black",size=12,face = "bold",angle = -315,vjust = 0.5,hjust = 0)) + 
  theme(axis.text.y=element_text(colour = "black",size = 18)) + 
  theme(axis.line = element_line(colour = "black"))+ 
  theme(axis.title.y = element_text(size = 24,face = "bold",
                                    margin = unit(c(0,1,0,1),"lines")))+
  #scale_y_continuous(limits = c(0,1),expand = c(0,0))+
  theme(legend.text = element_text(colour = "black",size = 15)) + 
  theme(legend.title = element_text(size = 18,colour = "black"))+
  theme(text = element_text(family = "Times"))
phylum
ggsave("./fig/Phylum.pdf",width = 6,height = 6)  

## 网络分析
library(edgeR)#用于基于edgeR的组间差异OTU识别
library(indicspecies)#用于基于indicator species的组间差异OTU识别
library(igraph)#用于共现性网络的绘制
library(Hmisc)#进行共现性网络构建前的OTU两两相关性计算
library(sciplot)#用于基于模块丰度的绘图
library(reshape2)#用于长宽数据转换
library(ggpmisc)#用于其它观测指标与组间差异模块丰度关系分析

#####输入OTU表格数据 #####
otu_its <- read.table("./data/network/JJ.txt",row.names=1,sep="\t",header=T, blank.lines.skip=F,check.names=F)
otu_its <- as.matrix(otu_its)
head(otu_its)

#####输入OTU物种分类数据 #####
tax_its <- read.table("./data/network/d2.txt",row.names=1, sep="\t", header=F,stringsAsFactors=F,quote="")
colnames(tax_its) <- c("Kingdom","Phylum","Class","Order", "Family", "Genus", "Species")
tax_its[tax_its==""] <- "unassigned"
head(tax_its)

#####输入metadata数据 #####
design_its <- read.table("./data/network/d3-JJ.txt", header=T, row.names=1, na.strings="NA")
#design_its$tillage <- factor(design_its $tillage,c("CPT","PT","ZT"))
#design_its$stage <- factor(design_its $stage ,c("T","J","F"))
#design_its$tilstag <- factor(design_its $tilstag,c("CPTT","CPTJ","CPTF","PTT","PTJ", "PTF","ZTT","ZTJ","ZTF"))
head(design_its)



#####以其中的T时期为例进行network的绘制#####
tsamples <- rownames(design_its)[which(design_its$stage == "JJ")]
tsamples
otu_its_t <- otu_its[,tsamples]
#去除一些低丰度OTU，在至少i个样本中具有至少j个序列
otu_its_t <- otu_its_t[which(rowSums(otu_its_t >= 1) >= 3),]
design_its_t <- droplevels(design_its[tsamples,])
tax_t_its <- tax_its[rownames(otu_its_t),]

#####寻找基于indicator species的组间显著差异OTU并保存结果#####
edgeR_its_t <- DGEList(counts=otu_its_t,
                       group=design_its_t$Group,
                       genes=tax_t_its)
#CPM标准化
otu_norm_its_t <- cpm(edgeR_its_t, normalized.lib.sizes=T, log=F)
# 准备数据：T时期下组间的indicator species
indic_t_its <- as.data.frame(t(otu_norm_its_t))
indic_t_groups_its <- design_its_t$Group
# 设置随机数种子，保证置换检验可重复
set.seed(8046)
# 鉴定各组指示种
indicatorsp_t_its <- multipatt(indic_t_its,indic_t_groups_its,func = "r.g",control=how(nperm=1000))
indic_t_df_its <- indicatorsp_t_its$sign


#按照阈值p.value < 0.05筛选各组显著的指示OTU
HRS_t_its <- as.matrix(indic_t_df_its[which(indic_t_df_its$s.HRS == 1 & indic_t_df_its$p.value < 0.05),])
ACRS_t_its <- as.matrix(indic_t_df_its[which(indic_t_df_its$s.ACRS == 1 & indic_t_df_its$p.value < 0.05),])
#合并
t_r_values_its <- rbind(HRS_t_its,ACRS_t_its)
# 组名修正，删除多余的"s."
# colnames(t_r_values_its)[1:3] <- c("CPT","PT","ZT")
colnames(t_r_values_its)[1:2] <- gsub("s.","",colnames(t_r_values_its)[1:2])


######绘制能够体现组件丰度差异OTU/模块的共现性网络图co-occurence network #####
#取基于indicator species分析和edgeR分析得到的显著组间差异OTU的交集
indic_edge_its_t <- rownames(t_r_values_its)
#基于TMM标准化的OTU表格进行OTU的两两Spearman相关计算
t_its_otu_cor <- rcorr(t(otu_norm_its_t),type=c("spearman"))
## 邻接矩阵转化为边表
CorrDF <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    from = rownames(cormat)[col(cormat)[ut]],
    to = rownames(cormat)[row(cormat)[ut]],
    cor =(cormat)[ut],
    p = pmat[ut]
  )
}
t_its_cor_df <- CorrDF(t_its_otu_cor$r,t_its_otu_cor$P)
# p值校正
t_its_cor_df$padj <- p.adjust(t_its_cor_df$p, method = "none") #method可选c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")
#取Spearman's rho > 0.7且p-value < 0.001的关系作为入选共现性网络co-occurence network的边
t_its_cor_df_padj <- t_its_cor_df[which(t_its_cor_df$cor > 0.80),]
t_its_cor_df_padj <- t_its_cor_df_padj[which(t_its_cor_df_padj$padj < 0.001),]

#生成node属性表
# 边两列合并为结点
nodeattrib_t_its <- data.frame(node = union(t_its_cor_df_padj$from,t_its_cor_df_padj$to))
# 显著的添加标签
nodeattrib_t_its$indicgroup <- 0
for (i in as.character(nodeattrib_t_its$node))
{
  if (i %in% indic_edge_its_t == TRUE)
  {nodeattrib_t_its[nodeattrib_t_its$node==i,"indicgroup"] <- paste(colnames(t_r_values_its)[which(t_r_values_its[i,1:2]==1)],collapse = "_")}
  else
  { nodeattrib_t_its[nodeattrib_t_its$node==i,"indicgroup"] <- "NA"}
}
#将OTU，即节点分类信息添加到node属性表
nodeattrib_t_its <- cbind(nodeattrib_t_its,tax_its[as.character(nodeattrib_t_its$node),])
#用igraph绘制共现性网络图co-occurence network
t_net_its <- graph_from_data_frame(t_its_cor_df_padj,direct=F, vertices = nodeattrib_t_its)
## 网络中的节点相对丰度
t_ra_its <- apply(otu_norm_its_t,1,mean)
t_ra_its <- t_ra_its[V(t_net_its)$name]
#将上述显著组间差异的OTU着色，这里有CPT,ZT,PT三种处理，因此理论上有6种可能的差异丰度分布情况
cs <- c("HRS","ACRS")
unique(V(t_net_its)$indicgroup)
V(t_net_its)$color <- V(t_net_its)$indicgroup
V(t_net_its)$color[!V(t_net_its)$color %in% cs] <- "gray30"
V(t_net_its)$color[V(t_net_its)$color == "HRS"] <- "#F8766D"
V(t_net_its)$color[V(t_net_its)$color == "ACRS"] <- "#619CFF"
V(t_net_its)$frame.color <- V(t_net_its)$color
#上述着色信息映射到node属性表中
tits_nodes <- rownames(nodeattrib_t_its[nodeattrib_t_its$indicgroup %in% cs,])
#设置节点形状
V(t_net_its)$shape <- "circle"
##设置节点大小，非显著组间差异OTU设置为"3"，显著的为"6"
V(t_net_its)$size <- V(t_net_its)$name
V(t_net_its)$size[!V(t_net_its)$size %in% tits_nodes] <- 3
V(t_net_its)$size[V(t_net_its)$size %in% tits_nodes] <- 6
tits_nodesizes <- as.numeric(V(t_net_its)$size)
#向量化各类型的显著组间差异OTU以便后续计算
HRS_nodes_t <- rownames(nodeattrib_t_its[nodeattrib_t_its$indicgroup=="HRS",])
ACRS_nodes_t <- rownames(nodeattrib_t_its[nodeattrib_t_its$indicgroup=="ACRS",])
cs_nodes_t_all <- c(HRS_nodes_t,ACRS_nodes_t)
#将网络中的节点/OTU进行聚类，这里采用fast greedy法
cfg_t <- cluster_fast_greedy(as.undirected(t_net_its))
#查看包含OTU数量最多的10个模块，以进行后续的着色
t_modules <- sort(table(membership(cfg_t)),decr=T)
t_modules_10 <- t_modules[1:10]
sm10_plot <- t_modules_10
names(sm10_plot) <- as.factor(1:10)
#寻找包含OTU数量最多的10个模块中具有显著组间差异OTU的模块
t_modules_cs <- table(factor(membership(cfg_t)[cs_nodes_t_all],levels=names(t_modules)))
t_modules_cs_10 <- t_modules_cs[names(t_modules_10)]
smcs10_plot <- t_modules_cs_10
names(smcs10_plot) <- as.factor(1:10)
#将OTU数量最多的10个模块中的OTU向量化
t_modules_points <- membership(cfg_t)[membership(cfg_t) %in% names(t_modules_10)]
t_points <- NULL
for(i in t_modules_points){
  tx <- which(names(t_modules_10)==i)
  t_points <- c(t_points, tx)
}
names(t_points) <- names(t_modules_points)
#按照组间差异OTU的类型着色这些OTU
t_all_cols <- sort(t_points)
t_all_cols[!names(t_all_cols) %in% cs] <- "gray30"
t_all_cols[names(t_all_cols) %in% HRS_nodes_t] <- "#F8766D"
t_all_cols[names(t_all_cols) %in% ACRS_nodes_t] <- "#619CFF"
#设置节点样式，1为空心圆代表普通节点，16为实心圆代表组间差异OTU
t_all_pch <- sort(t_points)
t_all_pch[names(t_all_pch) %in% rownames(otu_norm_its_t)] <- 1
t_all_pch[names(t_all_pch) %in% intersect(rownames(otu_norm_its_t),cs_nodes_t_all)] <- 16
#设置节点缩放倍数，1为普通节点1倍，2为组间差异OTU2倍
t_all_cex <- sort(t_points)
t_all_cex[!names(t_all_cex) %in% cs_nodes_t_all] <- 1
t_all_cex[names(t_all_cex) %in% cs_nodes_t_all] <- 2
#哪些模块包含有组间差异OTU
t_mods_list_cs <- list()
for (i in names(t_modules_cs_10)){
  x1 <- names(membership(cfg_t)[membership(cfg_t)==i])
  x2 <- x1[x1 %in% cs_nodes_t_all]
  t_mods_list_cs[[i]] <- as.numeric(V(t_net_its)[x2])
}
t_mods_list_cs
#设定layout出图，这里选择Fruchterman & Reingold
set.seed(8051)
coords_t_its <- layout_(t_net_its,with_fr(niter=9999, grid="nogrid"))
#每次运行耗费时间，可以把文件储存
#write.table(coords_t_its,paste0(output,"coords_t_its.txt"),sep="\t",row.names=F,col.names=F,quote=F)
#coords_t_its <- as.matrix(read.table("coords_t_its.txt"))
#dimnames(coords_t_its) <- NULL
#出图
pdf(paste0("./fig/network.pdf"),width=5,height=3)
par(mfrow=c(1,1), mar=c(0,0,0,0))
t_cols <- c("antiquewhite","antiquewhite2","antiquewhite4")
m_cols <- c("#F8766D","#619CFF")
plot(t_net_its,vertex.label=NA,vertex.size=tits_nodesizes, layout=coords_t_its,
     mark.groups=list(t_mods_list_cs$`6`,t_mods_list_cs$`4`,t_mods_list_cs$`1`),
     mark.col=t_cols, mark.border=t_cols)
legend("right",legend=c("Module 3", "Module 2", "Module 1"),col=t_cols,
       bty="n",fill=t_cols,border=t_cols)
legend("bottomright",legend = c("H","AC"),fill = m_cols)
dev.off()

pdf(paste0("./fig/module_abundance.pdf"),width=5,height=2)
par(mfrow=c(1,4), mar=c(0.5,3.5,2,0))#根据模块数调整mfrow
CS_cols <- c("#619CFF","#F8766D")
names(CS_cols) <- c("AC","H")
# T module 1
bargraph.CI(design_its_t$Group, colSums(otu_norm_its_t[cfg_t[[1]],])/1000,
            las=2, ylab="", cex.lab=.5, cex.axis=.7, cex.names=.7,
            err.width=.025, main="T M1", col=CS_cols, border=F)
mtext("cumulative 
      relative abundance", side = 2, line = 1, cex = 1.2)

# T module 2
bargraph.CI(design_its_t$Group, colSums(otu_norm_its_t[cfg_t[[4]],])/1000,
            las=2, ylab="", cex.lab=.5, cex.axis=.7, cex.names=.7,
            err.width=.025, main="T M2", col=CS_cols, border=F)

# T module 3
bargraph.CI(design_its_t$Group, colSums(otu_norm_its_t[cfg_t[[6]],])/1000,
            las=2, ylab="", cex.lab=.5, cex.axis=.7, cex.names=.7,
            err.width=.025, main="T M3", col=CS_cols, border=F)


plot.new()
par(mar=c(0.5,0,2,0))
legend("left", bty="n", cex=1, #x.intersp=0.1, y.intersp=1,
       legend=names(CS_cols),
       fill=CS_cols,
       border=CS_cols , xpd=T)
dev.off()  



cfg_t[[1]]
colSums(otu_norm_its_t[cfg_t[[1]],])/1000
write.csv(cfg_t[[1]],"jj_module1_otu.csv")

m1 <- read.table("./data/network/jj-m1.txt", header=T, sep="\t",dec=".")
write.table(table(m1$Genus),"./data/network/jj-m1-pie.txt",quote = F,row.names = F)


library(ggsci)
pie <- read.table("./data/network/jj-m1-pie.txt", header=T, sep="\t",dec=".")
p <- ggplot(data = pie, mapping = aes(x = 'Content', y = value,fill = Genus)) + 
  geom_bar(stat = 'identity', position = 'stack') +
  scale_fill_npg(alpha = 0.8)+
  theme_zg()
p

p + coord_polar(theta = 'y') + labs(x = '', y = '', title = '') + 
  theme(axis.text = element_blank()
  )


ggsave("./fig/pie_M1.pdf",width = 5,height = 3)



####### 平板抑菌效果
library(ggplot2)
library(tidyverse)
library(ggsci)
data <- read.table("./data/m1.txt",header = T,sep = "\t")
data <- arrange(data,M)
ggplot(data) +
  geom_bar(aes(M,mm,fill=M),stat = 'identity',width = 0.8) +
  geom_errorbar( aes(x=M, ymin=mm-sd, ymax=mm+sd), width=0.4, alpha=0.9, size=1.3) +
  theme_zg() +
  theme(legend.position = "null") + ##去掉图例
  scale_fill_npg(alpha = 0.8)+
  xlab("Antagonists")+
  ylab("Inhibitory zone(mm)")
ggsave("./fig/M1.pdf",width = 3,height = 3)

### 发病率
library(ggpubr)
library(ggplot2)
data <- read.table("./data/发病率.txt",header = T,sep = "\t")
Dis <-  ggplot(data) +
  geom_bar(aes(Strain,Disease,fill = Strain),stat = 'identity',width = 0.6) +
  geom_errorbar( aes(x=Strain, ymin=Disease-sd, ymax=Disease+sd), width=0.3, alpha=0.9, size=1.3) +
  theme_zg() +
  scale_fill_npg(alpha = 0.8)+
  theme(legend.position = "null") + ##去掉图例
  xlab("Antagonists")+
  ylab("Disease incidences(%)")
Dis
ggsave("./fig/M1_disease.pdf",width = 3,height = 3)






