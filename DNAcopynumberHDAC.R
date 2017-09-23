library(ggplot2)
library(gcookbook)
library(reshape)
library(dplyr)
data <- read.csv("D:/bioinformatics/Guochen/HD.csv", sep = ",", header=T,stringsAsFactors = F)
names(data)=c("Gene_ID","Well_ID","HDAC","Ctrl")
data$HDAC[data$HDAC%in% 0]<-1 
data$Ctrl[data$Ctrl%in% 0]<-1

data$Gene_ID=gsub("|HGLibB_","",data$Gene_ID)
data$Gene_ID=substr(data$Gene_ID,1,nchar(data$Gene_ID)-6)
data=data[1:57028,]
df=mutate(data,Score=log2(data$HDAC/data$Ctrl))
df2=df[,c(1,2,5)]
write.csv(df2,"D:/bioinformatics/Guochen/HDAC.csv", row.names = F)
dflog=read.csv("D:/bioinformatics/Guochen/outputhd.csv", sep = ",", header=T,stringsAsFactors = F)
str(dflog)
df6=aggregate(Score ~ Gene_ID, data=dflog,
              function(x){c(mean(x))})
df7=aggregate(LogP ~ Gene_ID, data=dflog,
              function(x){c(mean(x))})
df8=data.frame(df6, LogP=df7[,2])
df8$LogP=-df8$LogP
df8$Score=-df8$Score
df8$LogP[df8$LogP<0]=0

table(df8$LogP)

df8$LogP[df8$LogP<0 & df8$LogP==0]=NA
df8$LogP[df8$LogP==0.018]=NA
df8$LogP[df8$LogP==0.139]=NA
df8$LogP[df8$LogP==0.55]=NA
names(df8)=c("Gene_ID","Log2","Relative_Log_RSA")

newx=subset(df8,df8$Log2>4.5& df8$Relative_Log_RSA>2)

ggplot(df8, aes(x=Log2,y=Relative_Log_RSA))+geom_point() +annotate("text",x=newx$Log2,y=newx$Relative_Log_RSA+0.1,parse = T,label = newx$Gene_ID,size=3, family="serif",fontface="italic", colour="blue")

data.sel <- newx
write.csv(data.sel,"D:/bioinformatics/Guochen/HDACsubsetnew.csv", row.names = F)
x <- df8$Log2
y <- df8$Relative_Log_RSA

plot(x,y, pch=20, xlim=c(-10,10), ylim=c(0,5), xlab = "Fold change (log2 scaled)", ylab = "-Log RSA")
par(new = T)
pdf("guchenfigure.pdf")
x <- data.sel$Log2
y <- data.sel$Relative_Log_RSA
plot(x,y, pch=20, xlim=c(-10,10), ylim=c(0,5), xlab = "Fold change (log2 scaled)", ylab = "-Log RSA", col="red",cex=1.5)

x <- data.sel$Log2
y <- data.sel$Relative_Log_RSA
points(x,y,col = "blue", pch=21)


dev.off()

