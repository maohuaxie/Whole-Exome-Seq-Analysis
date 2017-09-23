library(ggplot2)
library(gcookbook)
library(reshape)

data <- read.csv("D:/bioinformatics/maohuaxie/ATRi.csv", sep = ",", header=T,stringsAsFactors = F)
data$DMSO..Untreated.[data$DMSO..Untreated. %in% 0]<-1 
data$VE822..500.nM.[data$VE822..500.nM. %in% 0]<-1

df=data[,c(1,3,4)]
names(df)=c("gene","ctrl","vE822")

df$gene=gsub("|HGLibB_","",df$gene)
df$gene=substr(df$gene,1,nchar(df$gene)-6)

dfctrl=rollapply(df$ctrl,width=3,FUN=mean,fill=NA)
dfve822=rollapply(df$vE822,width=3,FUN=mean,fill=NA)
data=data[1:57028,]
df2=df[1:57028,]
logchange=log(df2$vE822/df2$ctrl,2)
logchange1=log(df2$ctrl/df2$vE822,2)
df3=cbind(gene=df2$gene,logchange)
df4=cbind(Gene_ID=df2$gene,Well_ID=data$sequence,Score=logchange)
write.csv(df4,"D:/bioinformatics/Guochen/ve.csv", row.names = F)
df5=cbind(Gene_ID=df2$gene,Well_ID=data$sequence,Score=logchange1)
write.csv(df5,"D:/bioinformatics/Guochen/ve2.csv", row.names = F)




dflog=read.csv("D:/bioinformatics/maohuaxie/outputve2.csv", sep = ",", header=T,stringsAsFactors = F)
str(dflog)

df6=aggregate(Score ~ Gene_ID, data=dflog,
              function(x){c(mean(x))})
df7=aggregate(LogP ~ Gene_ID, data=dflog,
              function(x){c(mean(x))})
df8=data.frame(df6, LogP=df7[,2])
df8$LogP=-df8$LogP

df8$LogP[df8$LogP<0]=0

ggplot(df8, aes(x=Score,y=LogP))+geom_point()



dflog1=read.csv("D:/bioinformatics/maohuaxie/outputve1.csv", sep = ",", header=T,stringsAsFactors = F)

df9=aggregate(Score ~ Gene_ID, data=dflog1,
              function(x){c(mean(x))})
df10=aggregate(LogP ~ Gene_ID, data=dflog,
              function(x){c(mean(x))})
df11=data.frame(df9, LogP=df10[,2])
df11$LogP=-df11$LogP
df11$LogP[df11$LogP<0]=0


table(df11$LogP)

df11$LogP[df11$LogP<0 & df11$LogP==0]=NA
df11$LogP[df11$LogP==0.024]=NA
df11$LogP[df11$LogP==0.167]=NA
df11$LogP[df11$LogP==0.617]=NA
names(df11)=c("Gene_ID","Log2","Relative_Log_RSA")
newx=subset(df11,df11$Log2>4.5& df11$Relative_Log_RSA>2)
ggplot(df11, aes(x=Log2,y=Relative_Log_RSA))+geom_point() +annotate("text",x=newx$Log2,y=newx$Relative_Log_RSA+0.1,parse = T,label = newx$Gene_ID,size=3, family="serif",fontface="italic", colour="blue")

pl<-ggplot(df11,aes(x=Log2,y=Relative_Log_RSA))
pl2<- pl+ geom_point() + xlab("Log2 ") + ylab("-Log RSA") 
print(pl2)


data.sel <- newx
write.csv(data.sel,"D:/bioinformatics/maohuaxie/VEsubset.csv", row.names = F)
x <- df11$Log2
y <- df11$Relative_Log_RSA

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

