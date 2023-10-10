#########B

library(ggplot2)
##########致病性图##############
data <- read.table("./data/virulue-jj.txt",header=T,sep="\t")
ggplot(data,aes(Time,value,group = Treat,color = Treat))+
  geom_line()+
  geom_point(size = 0.5)+
  theme_classic()+
  theme(legend.position = c(0.20,0.75)) +
  guides(fill=guide_legend(title=NULL)) +
  geom_errorbar(aes(ymin=value-0.15*sd,
                    ymax=value+0.15*sd),
                width=0.2)+
  theme(panel.background = element_blank(),
        axis.line.x = element_line(),
        axis.line.y = element_line(),
        axis.title.x = element_text(size = 10,face = "bold"),  # 控制X轴标题字体的大小
        axis.text.x = element_text(size = 8),  # 控制X轴刻度标签字体的大小
        axis.title.y = element_text(size = 10,face = "bold"),  # 控制Y轴标题字体的大小
        axis.text.y = element_text(size = 8),  # 控制Y轴刻度标签字体的大小
        legend.text = element_text(size = 8))+
  scale_x_continuous(breaks = seq(0, 20, by = 4))+
  xlab("Days post inoculation(dpi)")+
  ylab("Disease index")
ggsave("./fig/virulue.pdf",width = 8.5,height = 5,units = "cm")



library(ggplot2)
library(ggpubr)

## 自然土壤共植实验
data <- read.table("./data/virulue1.txt",header=T,sep="\t")

A <- ggplot(data,aes(Time,value,group = Treat,color = Treat))+
  geom_line()+
  geom_point(size = 0.5)+
  theme_classic()+
  # xlim(0,12)+
  theme(legend.position = c(0.20,0.75)) +
  guides(fill=guide_legend(title=NULL)) +
  geom_errorbar(aes(ymin=value-0.15*sd,
                    ymax=value+0.15*sd),
                width=0.2)+
  theme(panel.background = element_blank(),
        axis.line.x = element_line(),
        axis.line.y = element_line(),
        axis.title.x = element_text(size = 10,face = "bold"),  # 控制X轴标题字体的大小
        axis.text.x = element_text(size = 8),  # 控制X轴刻度标签字体的大小
        axis.title.y = element_text(size = 10,face = "bold"),  # 控制Y轴标题字体的大小
        axis.text.y = element_text(size = 8),  # 控制Y轴刻度标签字体的大小
        legend.text = element_text(size = 8))+
  scale_x_continuous(breaks = seq(0, 30, by = 4))+
  xlab("Days post inoculation(dpi)")+
  ylab("Disease index")
A
ggsave("./fig/fig1.pdf",width = 8.5,height = 5,units = "cm")

## 基质土壤共植实验
data2 <- read.table("./data/virulue2.txt",header=T,sep="\t")

B<- ggplot(data2,aes(Time,value,group = Treat,color = Treat))+
  geom_line()+
  geom_point(size = 0.5)+
  theme_classic()+
  # xlim(0,12)+
  theme(legend.position = c(0.20,0.75)) +
  guides(fill=guide_legend(title=NULL)) +
  geom_errorbar(aes(ymin=value-0.15*sd,
                    ymax=value+0.15*sd),
                width=0.2)+
  theme(panel.background = element_blank(),
        axis.line.x = element_line(),
        axis.line.y = element_line(),
        axis.title.x = element_text(size = 10,face = "bold"),  # 控制X轴标题字体的大小
        axis.text.x = element_text(size = 8),  # 控制X轴刻度标签字体的大小
        axis.title.y = element_text(size = 10,face = "bold"),  # 控制Y轴标题字体的大小
        axis.text.y = element_text(size = 8),  # 控制Y轴刻度标签字体的大小
        legend.text = element_text(size = 8))+
  scale_x_continuous(breaks = seq(0, 30, by = 4))+
  xlab("Days post inoculation(dpi)")+
  ylab("Disease index")

ggsave("./fig/fig2.pdf",width = 8.5,height = 5,units = "cm")

