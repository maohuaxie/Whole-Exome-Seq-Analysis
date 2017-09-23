setwd("~/Project/20160413_maohua/")
data <- read.table("./data.csv", sep = ",", stringsAsFactors = F)
class(data[,4]) <- "numeric"
class(data[,6]) <- "numeric"

gene <- c("Il1rn", "Ccl3", "Ccl4", "Ccl12", "Cxcl12")
data.sel <- data[data[,10] %in% gene,]

x <- data[,4]
y <- -log(data[,6])/log(10)

plot(x,y, pch=20, xlim=c(-10,10), ylim=c(0,280), xlab = "", ylab = "")
par(new = T)

pdf("figure.pdf")
a = -10
x <- data[abs(data[,4]) > 1 & data[,6] < 10^a,4]
y <- -log(data[abs(data[,4]) > 1 & data[,6] < 10^a,6])/log(10)
plot(x,y, pch=20, xlim=c(-10,10), ylim=c(0,280), xlab = "Fold change (log2 scaled)", ylab = "-log(p.value)", col="red")

x <- data.sel[,4]
y <- -log(data.sel[,6])/log(10)
points(x,y,col = "blue", pch=21)

new.x <- c(-6, 3, 5, 6, 8)
new.y <- c(75, 120, 100, 80, 60)
arrows(new.x, new.y, x, y, col="orange", lwd=2)

text(new.x, new.y + 5, label=gene)
dev.off()
