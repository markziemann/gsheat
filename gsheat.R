library(reshape2)

tmp<-read.table("3col.txt",header=F)
x<-as.matrix(acast(tmp, V1~V2, value.var="V3"))
xt<-t(x)
xts<-scale(xt)
y<-t(xts)

##Fix cols and rows with zeros
y = y[rowSums(!is.na(y))!=0, colSums(!is.na(y))!=0]

covar <- sd(x)/rowMeans(x)
z <- x[ which( covar > 10),]
zt<-t(z)
zts<-scale(zt)
z<-t(zts)

z = z[rowSums(!is.na(z))!=0, colSums(!is.na(z))!=0]

#x<-as.matrix(x)
attach(as.data.frame(x))

mc=(mc1.rnk+mc2.rnk+mc3.rnk)/3
x<-cbind(x,mc=mc)

md=(md1.rnk+md2.rnk+md3.rnk)/3
x<-cbind(x,md=md)

mm=(mm1.rnk+mm2.rnk+mm3.rnk)/3
x<-cbind(x,mm=mm)

mp=(mp1.rnk+mp2.rnk+mp3.rnk)/3
x<-cbind(x,mp=mp)

ms=(ms1.rnk+ms2.rnk+ms3.rnk)/3
x<-cbind(x,ms=ms)

oi=(oi1.rnk+oi2.rnk+oi3.rnk)/3
x<-cbind(x,oi=oi)

rn=(rn1.rnk+rn2.rnk+rn3.rnk)/3
x<-cbind(x,rn=rn)

xx<-x[,grep("rnk",colnames(x),inver=T)]
x<-x[,grep("rnk",colnames(x))]

covar <- sd(xx)/rowMeans(xx)
z <- xx[ which( covar > 10),]
zt<-t(z)
zts<-scale(zt)
z<-t(zts)

z = z[rowSums(!is.na(z))!=0, colSums(!is.na(z))!=0]

library(limma)
contrasts<-read.table("contrast.mx",header=T,row.names=1)
design <- model.matrix(~ 0+factor(c(1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6,6,7,7,7)))
colnames(design) <- c("mc", "md", "mm", "mp","ms","oi","rn")

fit <- lmFit(x, design)
fit2 <- contrasts.fit(fit, contrasts)
fit2 <- eBayes(fit2)
pw<-topTable(fit2, number=100000, coef=ncol(contrasts), adjust="BH")
sig.pw<-pw[which (pw$adj.P.Val<0.05),]

sig.pw<-head(n=100,sig.pw)

sig.pw2<-as.data.frame(merge(x,sig.pw,by="row.names"))
sig.pw2<-sig.pw2[,1:(ncol(x)+1)]
row.names(sig.pw2)=sig.pw2$Row.names
sig.pw2$Row.names=NULL
rownames(sig.pw2)=gsub("REACTOME_","",rownames(sig.pw2))
rownames(sig.pw2)=gsub("_"," ",rownames(sig.pw2))
colnames(sig.pw2)=gsub(".rnk","",colnames(sig.pw2))

pdf("gsheat_100.pdf",width=7, height=8)
heatmap(as.matrix(sig.pw2), margins = c(4,20), cexRow=.4,main="most significant gene sets")
dev.off()


sig.pw<-head(n=50,sig.pw)
sig.pw2<-as.data.frame(merge(x,sig.pw,by="row.names"))
sig.pw2<-sig.pw2[,1:(ncol(x)+1)]
row.names(sig.pw2)=sig.pw2$Row.names
sig.pw2$Row.names=NULL
rownames(sig.pw2)=gsub("REACTOME_","",rownames(sig.pw2))
rownames(sig.pw2)=gsub("_"," ",rownames(sig.pw2))
colnames(sig.pw2)=gsub(".rnk","",colnames(sig.pw2))

pdf("gsheat_50.pdf",width=7, height=8)
heatmap(as.matrix(sig.pw2), margins = c(4,20), cexRow=.6,main="most significant gene sets")
dev.off()

