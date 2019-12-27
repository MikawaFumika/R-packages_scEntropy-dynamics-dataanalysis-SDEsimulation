setwd("/Users/yeyusong/Desktop/在投/entropydynamics/shuju")
library(readxl)
liberaydata <- read_excel('GSE122380_Supplementary_Data_Table_S1.xlsx', sheet=1, na='NA')
rawdata <- read.table("GSE122380_raw_counts.txt",fill=TRUE,header=F)
template <- read.table("mart_export.txt",fill=TRUE,header=F,sep=',')


#
setwd("/Users/yeyusong/Desktop/在投/entropydynamics/shuju/GSE122662_RAW")
GSM2836267_D0.matrix.mtx
reprogramming1 <- read.table("GSM2836267_D0.matrix.mtx",fill=TRUE,header=F)
gene1 <- read.table("GSM2836267_D0.genes.tsv",fill=TRUE,header=F)

library(Matrix)
matrix_dir = "/Users/yeyusong/Desktop/在投/entropydynamics/shuju/GSE122662_RAW/"
barcode.path <- paste0(matrix_dir, "GSM2836267_D0.barcodes.tsv.gz")
features.path <- paste0(matrix_dir, "features.tsv.gz")
matrix.path <- paste0(matrix_dir, "GSM2836267_D0.matrix.mtx.gz")
mat <- readMM(file = matrix.path)
feature.names = read.delim(features.path, 
                           header = FALSE,
                           stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path, 
                           header = FALSE,
                           stringsAsFactors = FALSE)
colnames(mat) = barcode.names$V1
rownames(mat) = feature.names$V1



#exon length template no intersect
geneid <- rawdata[,1]
geneid <- (as.character(geneid))
geneid <- geneid[-1]
exontemplate <- matrix(geneid)
exonlength <- vector()
for(i in 1:16319)
{
generaw=which(template[,1]==geneid[i])
genelengthtotal=as.numeric(as.character(template[generaw,6]))-as.numeric(as.character((template[generaw,5])))
exonlength[i]=sum(genelengthtotal)
}
exontemplate=data.frame(exontemplate,exonlength1)

#exon length template1 intersect
geneid <- rawdata[,1]
geneid <- (as.character(geneid))
geneid <- geneid[-1]
exontemplate1 <- matrix(geneid)
exonlength1 <- vector()
for(i in 1:16319)
{
  generaw=which(template[,1]==geneid[i])
  start=as.numeric(as.character((template[generaw,5])))
  end=as.numeric(as.character(template[generaw,6]))
  genelengthtotal=seq(start[1],end[1],1)
  if(length(generaw)>1)
  {
  for(j in 2:(length(generaw)))
    genelengthtotal=union(genelengthtotal,seq(start[j],end[j],1))
  }
  exonlength1[i]=length(genelengthtotal)
  print(i)
}
exontemplate1=data.frame(exontemplate1,exonlength1)
write.table(exontemplate1,"exontemplate1",quote = FALSE,row.names = FALSE,col.names = FALSE,sep="\t")

setwd("/Users/yeyusong/Desktop/在投/entropydynamics/shuju")
template <- read.table("exontemplate1",fill=TRUE,header=F,sep='\t')
rawdata <- read.table("GSE122380_raw_counts.txt",fill=TRUE,header=F)


# process tpm
exonlength1=template[,2]
rpkmdata=as.numeric(as.character(rawdata[,2]))[-1]/exonlength1
rpkmdata=rpkmdata/sum(rpkmdata)*1000000
rpkmdata=log((rpkmdata+1),2)
for(j in 3:298)
{  
  rpkm=as.numeric(as.character(rawdata[,j]))[-1]/exonlength1
  rpkm=rpkm/sum(rpkm)*1000000
  rpkm=log((rpkm+1),2)
  rpkmdata=data.frame(rpkmdata,rpkm)
  print(j)
}


#cauculate the distance
require(mixtools)
#define initial w
w=vector()
for(i in 1:297)
  w[i]=1/297

#xunhuan 
for(loop in 1:10)
{

#define reference set
referenceset=as.numeric(as.character(rpkmdata[,1]))*w[1]
for(i in 2:297)
{
  referenceset=referenceset+as.numeric(as.character(rpkmdata[,i]))*w[i]
}

hist(referenceset)

#define special reference set
i=0
cells=vector()
cellaaa=0
for(cell in 1:298)
{
  if(unlist((strsplit(as.character(rawdata[1,cell]),split="_")))[2]==i)
  {cellaaa=cellaaa+1
  cells[cellaaa]=cell-1}
}  
w=vector()
for(i in 1:297)
  w[i]=0
for(i in cells)
  w[i]=1
w=w/length(cells)


#distance
distancestat=vector()
for(i in 1:297)
{
y=as.numeric(as.character(rpkmdata[,i]))-referenceset
disy=hist(y,xlim=c(-10,10),breaks=1000)
counts=disy$counts
counts=counts/sum(counts)
distance=0
for(j in 1:length(counts))
{
  if(counts[j]!=0)
    distance <- distance + counts[j]*log(counts[j])
}
distancestat[i]=-distance
print(i)
}
hist(distancestat,breaks=10)

#em reckon
em <- normalmixEM(distancestat)
print(em$lambda)
print(em$mu)
print(em$sigma)
plot(em, whichplots = 2,breaks=30)


for(i in 1:297)
{w[i]=dnorm(distancestat[i],(em$mu)[1],(em$sigma)[1])}
w=w/sum(w)
}

#define days
i=0
cells=vector()
cellaaa=0
for(cell in 1:298)
{
  if(unlist((strsplit(as.character(rawdata[1,cell]),split="_")))[2]==i)
  {cellaaa=cellaaa+1
  cells[cellaaa]=cell-1}
}  
Dayscell=list(cells)

for(i in 1:15)
{
  cells=vector()
  cellaaa=0
  for(cell in 1:298)
  {
    if(unlist((strsplit(as.character(rawdata[1,cell]),split="_")))[2]==i)
    {cellaaa=cellaaa+1
    cells[cellaaa]=cell-1}
  }
  cells=list(cells)
  Dayscell=c(Dayscell,cells)
}

#plot distance all
par(mar=c(3,3,1,1))
colors <- colorRampPalette(c(rgb(0/255,144/255,255/255), rgb(238/255,44/255,44/255)))(16)
for(i in 0:15)
{
  cells=unlist(Dayscell[i+1])
  plot((runif(length(distancestat[cells]))/2+0.75+i),distancestat[cells],xlim=c(0,17),ylim=c(4.5,6),pch=19,cex=0.5,col=colors[i+1],xlab="", ylab="", main="",xaxt="n",yaxt="n",mgp=c(0,0.5,0),tck=0.03)
  par(new=TRUE)
}
colors <- colorRampPalette(c("pink", "red"))(8)
for(i in 0:7)
{
  cells=unlist(cellb[i+1])
  plot((runif(length(distancestat[cells]))/2+0.75+i),distancestat[cells],xlim=c(0,17),ylim=c(4.5,6),pch=19,cex=0.5,col=colors[i+1],xlab="", ylab="", main="",xaxt="n",yaxt="n",mgp=c(0,0.5,0),tck=0.03)
  par(new=TRUE)
}
axis(side=1,las=1,at=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16),mgp=c(0,0.5,0),tck=0.03,las=1,labels=c("0","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15"),cex.axis=0.8)
axis(side=2,las=1,at=c(4.5,5,5.5,6),mgp=c(0,0.5,0),tck=0.03,las=1,labels=c("4.5","5","5.5","6"),cex.axis=1)
mtext("Days",side=1,line=1.5,cex=1)
mtext("scEntropy",side=2,line=1.7,cex=1)

#jun zhi均值fang cha
meandistancea=vector()
vardistancea=vector()
for(i in 1:16)
{  
  meandistancea[i]=mean(distancestat[unlist(Dayscell[i])])
  vardistancea[i]=var(distancestat[unlist(Dayscell[i])])
}

par(mar=c(3,3,1,1))
colors <- colorRampPalette(c(rgb(0/255,144/255,255/255), rgb(238/255,44/255,44/255)))(16)
plot(seq(0,15),meandistancea,type="p",pch=19,tck=0.03,cex=0.8,las=1,xlab="",col=colors[rank(1:16)], ylab="", main="",xlim=c(0,15),ylim=c(4.5,6),xaxt="n",yaxt="n",bty="l")
mtext("Days",side=1,line=0.5,cex=1)
mtext("Average scEntropy",side=2,line=1.5,cex=1)
axis(side=1,las=1,at=c(0,5,10,15),mgp=c(0,0.5,0),tck=0.03,las=1,labels=c("","","",""),cex.axis=1)
axis(side=2,las=1,at=c(4.5,5,5.5,6),mgp=c(0,0.5,0),tck=0.03,las=1,labels=c("4.5","5","5.5","6"),cex.axis=1)

par(mar=c(3,3.5,1,1))
colors <- colorRampPalette(c(rgb(0/255,144/255,255/255), rgb(238/255,44/255,44/255)))(16)
plot(seq(0,15),vardistancea,type="p",pch=19,tck=0.03,cex=0.8,las=1,xlab="",col=colors[rank(1:16)], ylab="", main="",xlim=c(0,15),ylim=c(0,0.12),xaxt="n",yaxt="n",bty="l")
mtext("Days",side=1,line=0.5,cex=1)
mtext("Variance of scEntropy",side=2,line=2,cex=1)
axis(side=1,las=1,at=c(0,5,10,15),mgp=c(0,0.5,0),tck=0.03,las=1,labels=c("","","",""),cex.axis=1)
axis(side=2,las=1,at=c(0,0.05,0.1),mgp=c(0,0.5,0),tck=0.03,las=1,labels=c(0,0.05,0.1),cex.axis=1)

#zhengtifenbu
hist(distancestat,xlim=c(4.5,6),col='gray',tck=0.03,las=1,xlab="", ylab="",xaxt="n",bty="l",main="",breaks=seq(4.5,6,0.1))
axis(side=1,las=1,at=c(4.5,5,5.5,6),mgp=c(0,0.5,0),tck=0,las=1,labels=c(4.5,5,5.5,6),cex.axis=1)
mtext("sc-Entropy",side=1,line=1.5,cex=1)


#force atlas
library(networkD3)
data(MisLinks)
data(MisNodes)
forceNetwork(Links = MisLinks, Nodes = MisNodes, Source = "source",
             Target = "target", Value = "value", NodeID = "name",
             Group = "group", opacity = 0.4, zoom = TRUE)


matrixforce=matrix(nrow=297,ncol=297)
for(i in 1:297)
{
  for (j in 1:297)
  {
    deltarpkm=rpkmdata[,i]-rpkmdata[,j]
    deltarpkm=hist(deltarpkm,xlim=c(-10,10),breaks=200)
    deltarpkm=deltarpkm$counts/(sum(deltarpkm$counts))
    matrixforce[i,j]=sum(deltarpkm*(log(deltarpkm+1)))
  }
  print(i)
}

axismatrix=matrix(nrow=297,ncol=2)

for(i in 1:297)
{
  for(j in 1:2)
  {
axismatrix[i,j]=runif(1)
}
}
#iteration
k=1
kr=0.1

#try
#disx=c(1)
#disy=c(4)
#distancestat=c(rep(disx,times=100),rep(disy,times=197))

iteration=30000
  
for(iter in 1:iteration)
{ 
  samplex=sample(seq(1,297,1),1)
  sampley=sample(seq(1,297,1),1)
  if(samplex!=sampley)
  { 
    currentdistance=((axismatrix[samplex,2]-axismatrix[sampley,2])^2+(axismatrix[samplex,1]-axismatrix[sampley,1])^2)^0.5
    #attraction=abs(1*currentdistance)/((matrixforce[samplex,sampley])^2)
    #repulsion=1/currentdistance^2
    repulsion=((distancestat[samplex]-distancestat[sampley])^2)/currentdistance^2
    attraction=abs(1*currentdistance)/((distancestat[samplex]-distancestat[sampley])^2)

    
    length=(repulsion-attraction)*kr
    if(length>0.5)
    {
      length=0.5
    }
    if(length<(-0.5))
    {
      length=-0.5
    }

    xlength=length*(axismatrix[samplex,1]-axismatrix[sampley,1])
    ylength=length*(axismatrix[samplex,2]-axismatrix[sampley,2])

    axismatrix[samplex,1]=axismatrix[samplex,1]+xlength
    axismatrix[samplex,2]=axismatrix[samplex,1]+ylength
  }
}
plot(axismatrix[,1],axismatrix[,2],cex=0.2,xlim=c(min(axismatrix[,1]),max(axismatrix[,1])),ylim=c(min(axismatrix[,2]),max(axismatrix[,2])))



colors <- colorRampPalette(c(rgb(0/255,144/255,255/255), rgb(238/255,44/255,44/255)))(16)
par(mar=c(1,1,1,1))
for(i in seq(0,15,1))
{
  i=i+1
  cells=unlist(Dayscell[i+1])
  plot(axismatrix[cells,1],axismatrix[cells,2],cex=0.2,xlim=c(0.2,0.8),ylim=c(0.2,0.8),col=colors[i+1],xlab="",ylab="",xaxt='n',yaxt='n',bty="l")
  par(new=TRUE)
  #Sys.sleep(1)
}
#pingjun guixian
for(i in seq(0,14,1))
{
  i=i+1
  cells=unlist(Dayscell[i+1])
  cells1=unlist(Dayscell[i+2])
  plot(c(mean(axismatrix[cells,1]),mean(axismatrix[cells1,1])),c(mean(axismatrix[cells,2]),mean(axismatrix[cells1,2])),type='l',lwd=3.0,xlim=c(0.2,0.8),ylim=c(0.2,0.8),col='darkgrey',xlab="",ylab="",xaxt='n',yaxt='n',bty="l")
  par(new=TRUE)
  #Sys.sleep(1)
}
#moniguixian
for(i in seq(2,15,1))
{
  plot(c(mean(axismatrix[scentropysimutotal[,i],1]),mean(axismatrix[scentropysimutotal[,i+1],1])),c(mean(axismatrix[scentropysimutotal[,i],2]),mean(axismatrix[scentropysimutotal[,i+1],2])),type='l',lwd=2.0,xlim=c(0.2,0.8),ylim=c(0.2,0.8),col='darkgrey',xlab="",ylab="",bty='l',xaxt='n',yaxt='n')
  par(new=TRUE)
  #Sys.sleep(1)
}


i=0
cells=unlist(Dayscell[i+1])
plot(axismatrix[cells,1],axismatrix[cells,2],cex=0.2,xlim=c(-30,30),ylim=c(-30,30),col=colors[i+1])
par(new=T)


#plot distance distribution
meandistance=vector()
vardistance=vector()
par(mfrow=c(4,4),mar=c(2,2,1,1))
colors <- colorRampPalette(c("blue", "red"))(16)
distance_density_change=vector()
for(i in 1:15)
  distance_density_change[i]=0
distance_density_change=data.frame(distance_density_change)

distance_mid_change=vector()
for(i in 1:15)
  distance_mid_change[i]=0
distance_mid_change=data.frame(distance_mid_change)

for(i in 0:15)
{
cells=vector()
cellaaa=0
for(cell in 1:298)
 {
  if(unlist((strsplit(as.character(rawdata[1,cell]),split="_")))[2]==i)
    {cellaaa=cellaaa+1
    cells[cellaaa]=cell-1}
 }  
#plot(distancestat[cells],runif(length(distancestat[cells])),xlim=c(3.5,5.5),ylim=c(0,1),breaks=10,pch=19,col=colors[i+1],xlab="", ylab="", main="")


meandistance[i+1]=mean(distancestat[cells])
vardistance[i+1]=var(distancestat[cells])


distancedistri=hist(distancestat[cells],xlim=c(4.5,6),col=colors[i+1],ylim=c(0,15),tck=0.03,las=1,xlab="", ylab="", main=(i),xaxt="n",yaxt="n",bty="l",breaks=seq(4.5,6,0.1))
axis(side=1,las=1,at=c(4.5,5,5.5,6),mgp=c(0,0.5,0),tck=0.03,las=1,labels=c("4.5","5","5.5","6"),cex.axis=1)
axis(side=2,las=1,at=c(0,5,10,15),mgp=c(0,0.5,0),tck=0.03,las=1,labels=c("0","5","10","15"),cex.axis=1)
distance_density_change=data.frame(distance_density_change,distancedistri$counts/sum(distancedistri$counts))
distance_mid_change=data.frame(distance_mid_change,distancedistri$mids)
}

par(mar=c(2,2,1,1))
hist(distancestat,xlim=c(4.5,6),ylim=c(0,80),tck=0.03,las=1,col="grey",xlab="sc-Entropy", ylab="",main="Distribution of 297 entropy",xaxt="n",yaxt="n",bty="l")
axis(side=1,las=1,at=c(4.5,5,5.5,6),mgp=c(0,0.5,0),tck=0.03,las=1,labels=c("4.5","5","5.5","6"),cex.axis=1)
axis(side=2,las=1,at=c(0,40,80),mgp=c(0,0.5,0),tck=0.03,las=1,labels=c("0","40","80"),cex.axis=1)

#fen lei分类 mean var

meandistancea=vector()
vardistancea=vector()
meandistanceb=vector()
vardistanceb=vector()
for(i in 1:8)
{  
  meandistanceb[i]=mean(distancestat[unlist(cellb[i])])
  vardistanceb[i]=var(distancestat[unlist(cellb[i])])
}
for(i in 1:16)
{  
  meandistancea[i]=mean(distancestat[unlist(cella[i])])
  vardistancea[i]=var(distancestat[unlist(cella[i])])
}
vardistanceb[8]=0

par(mar=c(3,3,1,1))
colors <- colorRampPalette(c("blue", "red"))(16)
plot(seq(0,15),meandistance,type="p",pch=19,tck=0.03,cex=0.8,las=1,xlab="",col=colors[rank(1:16)], ylab="", main="",xlim=c(0,15),ylim=c(4.5,6),xaxt="n",yaxt="n",bty="l")
mtext("Days",side=1,line=1.5,cex=1)
mtext("Average sc-Entropy",side=2,line=1.5,cex=1)
axis(side=1,las=1,at=c(0,5,10,15),mgp=c(0,0.5,0),tck=0.03,las=1,labels=c("0","5","10","15"),cex.axis=1)
axis(side=2,las=1,at=c(4.5,5,5.5,6),mgp=c(0,0.5,0),tck=0.03,las=1,labels=c("4.5","5","5.5","6"),cex.axis=1)

#fen l分类hua分类话
par(mar=c(3,3,1,1))
colors <- colorRampPalette(c("black", "blue"))(16)
plot(seq(0,15),meandistancea,type="p",pch=19,tck=0.03,cex=0.8,las=1,xlab="",col=colors[rank(1:16)], ylab="", main="",xlim=c(0,15),ylim=c(4.5,6),xaxt="n",yaxt="n",bty="l")
par(new=TRUE)
colors <- colorRampPalette(c("pink", "red"))(8)
plot(seq(0,7),meandistanceb,type="p",pch=19,tck=0.03,cex=0.8,las=1,xlab="",col=colors[rank(1:8)], ylab="", main="",xlim=c(0,15),ylim=c(4.5,6),xaxt="n",yaxt="n",bty="l")
mtext("Days",side=1,line=1.5,cex=1)
mtext("Average sc-Entropy",side=2,line=1.5,cex=1)
axis(side=1,las=1,at=c(0,5,10,15),mgp=c(0,0.5,0),tck=0.03,las=1,labels=c("0","5","10","15"),cex.axis=1)
axis(side=2,las=1,at=c(4.5,5,5.5,6),mgp=c(0,0.5,0),tck=0.03,las=1,labels=c("4.5","5","5.5","6"),cex.axis=1)


par(mar=c(3,3,1,1))
colors <- colorRampPalette(c("blue", "red"))(16)
plot(seq(0,15),vardistance,type="p",pch=19,tck=0.03,cex=0.8,las=1,xlab="",col=colors[rank(1:16)], ylab="", main="",xlim=c(0,15),ylim=c(0,0.1),xaxt="n",yaxt="n",bty="l")
mtext("Days",side=1,line=1.5,cex=1)
mtext("Var",side=2,line=2.0,cex=1)
axis(side=1,las=1,at=c(0,5,10,15),mgp=c(0,0.5,0),tck=0.03,las=1,labels=c("0","5","10","15"),cex.axis=1)
axis(side=2,las=1,at=c(0,0.05,0.1),mgp=c(0,0.5,0),tck=0.03,las=1,labels=c("0","0.05","0.1"),cex.axis=1)

#fen l分类hua分类话
par(mar=c(3,3,1,1))
colors <- colorRampPalette(c("black", "blue"))(16)
plot(seq(0,15),vardistancea,type="p",pch=19,tck=0.03,cex=0.8,las=1,xlab="",col=colors[rank(1:16)], ylab="", main="",xlim=c(0,15),ylim=c(0,0.1),xaxt="n",yaxt="n",bty="l")
par(new=TRUE)
colors <- colorRampPalette(c("pink", "red"))(8)
plot(seq(0,7),vardistanceb,type="p",pch=19,tck=0.03,cex=0.8,las=1,xlab="",col=colors[rank(1:8)], ylab="", main="",xlim=c(0,15),ylim=c(0,0.1),xaxt="n",yaxt="n",bty="l")
mtext("Days",side=1,line=1.5,cex=1)
mtext("Var",side=2,line=2.0,cex=1)
axis(side=1,las=1,at=c(0,5,10,15),mgp=c(0,0.5,0),tck=0.03,las=1,labels=c("0","5","10","15"),cex.axis=1)
axis(side=2,las=1,at=c(0,0.05,0.1),mgp=c(0,0.5,0),tck=0.03,las=1,labels=c("0","0.05","0.1"),cex.axis=1)




plot(seq(0,15),vardistance)

hist(w)

#corre with gene expression
orderofdistance=order(distancestat)
sortofdistance=sort(distancestat)
rpkmdatasort=rpkmdata[,which(distancestat==min(distancestat))]
for(i in 2:297)
  {
  rpkmdatasort1=rpkmdata[,orderofdistance[i]]
  rpkmdatasort=data.frame(rpkmdatasort,rpkmdatasort1)
}
pearson=vector()
for(i in 1:16319)
{
persontest=cor.test(sortofdistance,as.numeric(as.character(rpkmdatasort[i,])),method="pearson")
pearson[i]=persontest$estimate
print(i)
}
hist(pearson,col='grey',xlim=c(-1,1))
yuzhi=0.73
length(pearson[(abs(pearson))>yuzhi])
which(pearson>yuzhi | pearson<(-yuzhi))
pearson[which(pearson>yuzhi | pearson<(-yuzhi))]
template[which(pearson>yuzhi | pearson<(-yuzhi)),1]

#guanyu jiyinde bao yong biomanager
library("AnnotationDbi")
library("org.Hs.eg.db")

GENEID <- seq(1,11,1)
ENSEMBL <- c('ENSG00000008516','ENSG00000075975','ENSG00000088305','ENSG00000099940','ENSG00000105131','ENSG00000130182','ENSG00000164362','ENSG00000174099','ENSG00000177981','ENSG00000183765','ENSG00000198598')
df <- data.frame(ENSEMBL,GENEID)
df$symbol <- mapIds(org.Hs.eg.db,keys=ENSEMBL,column="SYMBOL",keytype="ENSEMBL",multiVals="first")




density(pearson,col="grey")
rawdata[(162+1),1]

ccc=c(194 , 1206 , 1604 , 1937 , 2686 , 5535 ,10079 ,11931 ,12457 ,13216 ,14692)
pearsonzhishu=c(-0.7309178,0.7338711,-0.7330128 ,0.7300553,-0.7310438,-0.7349551,-0.7335726 ,0.7315132, 0.7538376, -0.7321125, -0.7550320)
for(i in 1:10)
print(rawdata[ccc[i],1])
#194  1206  1604  1937  2686  5535 10079 11931 12457 13216 14692

colors <- colorRampPalette(c("blue", "red"))(297)
#plot cor plot
a=c(cell1a,cell2a,cell3a,cell4a,cell5a,cell6a,cell7a,cell8a,cell9a,cell10a,cell11a,cell12a,cell13a,cell14a,cell15a,cell16a)
b=c(cell1b,cell2b,cell3b,cell4b,cell5b,cell6b,cell7b,cell8b)

colors <- colorRampPalette(c("pink", "red"))(26)

par(mfrow=c(3,3),mar=c(2,2,2,1))
for(i in 1:9)
{
per=ccc[i]
plot(sortofdistance[a],rpkmdatasort[per,][a],xlim=c(4.5,6),ylim=range(rpkmdatasort[per,][a]),pch=19,cex=0.4,col='black',xlab="", ylab="", main=pearsonzhishu[i],mgp=c(0,0.5,0),tck=0.03)
par(new=T)
plot(sortofdistance[b],rpkmdatasort[per,][b],xlim=c(4.5,6),ylim=range(rpkmdatasort[per,][a]),pch=19,cex=0.4,col='red',xlab="", ylab="", main="",mgp=c(0,0.5,0),tck=0.03)
}

#clusting corre with gene expression 有错（

orderofdistancea=order(distancestat[c(cell1a,cell2a,cell3a,cell4a,cell5a,cell6a,cell7a,cell8a,cell9a,cell10a,cell11a,cell12a,cell13a,cell14a,cell15a,cell16a)])
sortofdistancea=sort(distancestat[c(cell1a,cell2a,cell3a,cell4a,cell5a,cell6a,cell7a,cell8a,cell9a,cell10a,cell11a,cell12a,cell13a,cell14a,cell15a,cell16a)])
rpkmdatasorta=rpkmdata
dui ying d对应的bu shi对应的不是rpkm
[,cell1a[which(distancestat[c(cell1a,cell2a,cell3a,cell4a,cell5a,cell6a,cell7a,cell8a,cell9a,cell10a,cell11a,cell12a,cell13a,cell14a,cell15a,cell16a)]==min(distancestat[c(cell1a,cell2a,cell3a,cell4a,cell5a,cell6a,cell7a,cell8a,cell9a,cell10a,cell11a,cell12a,cell13a,cell14a,cell15a,cell16a)]))]]
for(i in 2:271)
{
  rpkmdatasort1=rpkmdata[,orderofdistancea[i]]
  rpkmdatasorta=data.frame(rpkmdatasorta,rpkmdatasort1)
}
pearsona=vector()
for(i in 1:16319)
{
  persontest=cor.test(sortofdistancea,as.numeric(as.character(rpkmdatasorta[i,])),method="pearson")
  pearsona[i]=persontest$estimate
  print(i)
}
hist(pearsona)
yuzhi=0.39
length(pearsona[(abs(pearsona))>yuzhi])
which(pearsona>yuzhi | pearsona<(-yuzhi))
pearsona[which(pearsona>yuzhi | pearsona<(-yuzhi))]

orderofdistanceb=order(distancestat[c(cell1b,cell2b,cell3b,cell4b,cell5b,cell6b,cell7b,cell8b)])
sortofdistanceb=sort(distancestat[c(cell1b,cell2b,cell3b,cell4b,cell5b,cell6b,cell7b,cell8b)])
rpkmdatasortb=rpkmdata[,cell1b[which(distancestat[cell1b]==min(distancestat[cell1b]))]]
for(i in 2:26)
{
  rpkmdatasort1=rpkmdata[,orderofdistanceb[i]]
  rpkmdatasortb=data.frame(rpkmdatasortb,rpkmdatasort1)
}
pearsonb=vector()
for(i in 1:16319)
{
  persontest=cor.test(sortofdistanceb,as.numeric(as.character(rpkmdatasortb[i,])),method="pearson")
  pearsonb[i]=persontest$estimate
  print(i)
}
hist(pearsonb)
yuzhi=0.61
length(pearsonb[(abs(pearsonb))>yuzhi])
which(pearsonb>yuzhi | pearsonb<(-yuzhi))
pearson[which(pearsonb>yuzhi | pearsonb<(-yuzhi))]


#plot rpkm distri
par(mfrow=c(4,4),mar=c(1,1,1,1))
hist(rpkmdata[,1],xlim=c(-50,0),breaks=100,ylim=c(0,1000))
hist(rpkmdata[,9],xlim=c(-50,0),breaks=100,ylim=c(0,1000))
hist(rpkmdata[,10],xlim=c(-50,0),breaks=100,ylim=c(0,1000))
hist(rpkmdata[,11],xlim=c(-50,0),breaks=100,ylim=c(0,1000))
hist(rpkmdata[,12],xlim=c(-50,0),breaks=100,ylim=c(0,1000))
hist(rpkmdata[,13],xlim=c(-50,0),breaks=100,ylim=c(0,1000))
hist(rpkmdata[,14],xlim=c(-50,0),breaks=100,ylim=c(0,1000))
hist(rpkmdata[,15],xlim=c(-50,0),breaks=100,ylim=c(0,1000))
hist(rpkmdata[,16],xlim=c(-50,0),breaks=100,ylim=c(0,1000))
hist(rpkmdata[,2],xlim=c(-50,0),breaks=100,ylim=c(0,1000))
hist(rpkmdata[,3],xlim=c(-50,0),breaks=100,ylim=c(0,1000))
hist(rpkmdata[,4],xlim=c(-50,0),breaks=100,ylim=c(0,1000))
hist(rpkmdata[,5],xlim=c(-50,0),breaks=100,ylim=c(0,1000))
hist(rpkmdata[,6],xlim=c(-50,0),breaks=100,ylim=c(0,1000))
hist(rpkmdata[,7],xlim=c(-50,0),breaks=100,ylim=c(0,1000))


cellgeneexpression=rawdata[,8]
cellgeneexpression=as.numeric(as.character(cellgeneexpression))
cellgeneexpression=cellgeneexpression[-1]
hist(cellgeneexpression,xlim=c(0,5000),breaks=500)

n <- 1000
mean_s <- c(-1, 7)
y <- sample(c("head", "tail"), size = n, replace = TRUE, prob = c(0.4, 0.6))
x <- rnorm(n = 1000, mean = mean_s[1])
tails <- y %in% c("tail")
x[tails] <- rnorm(sum(tails), mean = mean_s[2],sd=2)

require(mixtools)
em <- normalmixEM(x)
print(em$lambda)
print(em$mu)
print(em$sigma)
plot(em, whichplots = 2)

hist(x)
par(new=TRUE)
plot(em, whichplots = 2,xlim=c(-10,10))


source("http://bioconductor.org/biocLite.R")
biocLite('biomaRt')
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()
install.packages('DT')
BiocManager::install("edgeR")
library(DT)
library(BiocManager)
library("edgeR")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("SNPlocs.Hsapiens.dbSNP142.GRCh37")



#nihe


-(x-x1)*(3-1/(1+160/step))


B=0.05
k=0.1
xinitial1=5.7
xinitial2=4.7

simufpe=vector()
for(i in 1:151)
  simufpe[i]=0
simufpe=data.frame(simufpe)

par(mfrow=c(4,4),mar=c(2,2,2,1))
for(time1 in 1:16)
{ 
  k=0.1*(2-1/(1+10/time1))
  
  mean=(xinitial1-5.6)*exp(-k*(time1))
  sd= (B^2/k/2*(1-exp(-2*k*(time1))))^0.5
  simufpehist1=dnorm(seq(4.5,6,0.01)-5.6,mean,sd)
  simufpehist1=simufpehist1/sum(simufpehist1)*10
  
  mean=(xinitial2-5.6)*exp(-k*(time1))
  simufpehist2=dnorm(seq(4.5,6,0.01)-5.6,mean,sd)
  simufpehist2=simufpehist2/sum(simufpehist2)*10
  
  simufpehist=0.25*simufpehist1+0.75*simufpehist2
  simufpe=data.frame(simufpe,simufpehist)
}

par(mfrow=c(4,4),mar=c(2,2,2,1))
for(i in 1:16)
 { 
  plot(seq(4.5,6,0.01),simufpe[,i+1],type='l',xlim=c(4.5,6),ylim=c(0,0.6),col='red',tck=0.03,mgp=c(0,0.5,0),las=1,xlab="", ylab="", main="")
  par(new=T)
  for(j in 1:15)
  {
    polygon(c(4.4+0.1*j+0.01,4.5+0.1*j-0.01,4.5+0.1*j-0.01,4.4+0.1*j+0.01),c(distance_density_change[j,i+1],distance_density_change[j,i+1],0,0), density = NULL, border = F, col = rgb(150/255,150/255,150/255),tck=0.03,mgp=c(0,0.5,0),las=1,xlab="", ylab="", main="")
    par(new=T)
  }
  plot(seq(4.5,6,0.01),simufpe[,i+1],type='l',xlim=c(4.5,6),ylim=c(0,0.6),col='red',tck=0.03,mgp=c(0,0.5,0),las=1,xlab="", ylab="", main="")
}





#numerical simulation

k1=0.15
B1=0.05
xinitial1=-0.9
xcenter=5.6
time1=1

k2=0.4
B2=0.1
xinitial2=-0.3
xcenter2=5.8

a1=0.5
a2=2




ynumeric=vector()

if((rnorm(1,0,1)>0)){
ynumeric[1]=5.5
color='red'
holdvalue=0
}else 
{ynumeric[1]=4.7
color='black'
holdvalue=0}


steps=10
sddeltat=(1/steps)^0.5

for(step in 1:(steps*16))
{
  mt=a1/(1+step/a2/steps)
  deltax=(-k1*(ynumeric[step]-xcenter))*(1-mt)+(-k2*(ynumeric[step]-xcenter2))*(mt)
  #deltaguassian=(1-mt)*B1*rnorm(1,mean=0,sd=sddeltat) +(mt)*B2*rnorm(1,mean=0,sd=sddeltat)
  deltaguassian=rnorm(1,mean=0,sd=sddeltat*((1-mt)*B1+(mt)*B2))
  ynumeric[step+1]=ynumeric[step]+deltax*1/steps+deltaguassian
  print(step/steps)
  if(runif(1)<holdvalue)
    break
}
xnumeric=seq(0,length(ynumeric)-1,1)
plot(xnumeric/steps,ynumeric,type='l',col=color,xlim=c(0,16),ylim=c(4.5,6),tck=0.03,mgp=c(0,0.5,0),las=1,xlab="", ylab="", main="")
par(new=T)

mtext("Days",side=1,line=1.5,cex=1)
mtext("sc-Entropy",side=2,line=2.0,cex=1)
legend(10,5.2,legend=c("Subgroup1","Subgroup2"),cex=1.0,x.intersp=0.5,y.intersp=1.5,col=c('red','black'),bty="n",lty=c(16,16))




#u function
x1=4.8
x2=5.2
x3=5.7
B=000.1
A=1
x=seq(0,10,0.01)
upie=(x-x1)*(x-x2)*(x-x3)
ufunction=0.25*x^4-1/3*x3*x^3-1/3*(x1+x2)*x^3+0.5*(x1+x2)*x3*x^2+x1*x2*0.5*x^2-x1*x2*x3*x
ufunction=-ufunction/B+log(B)
pro=A*exp(-ufunction/B)
pro=pro/(sum(pro))*100
plot(x,pro,cex=0.2)
plot(x,ufunction,cex=0.2)

#u simulation
x1=5.6
B=0.05

ynumeric=vector()

if((runif(1)>0.75)){
  ynumeric[1]=5.6
  color='red'
}else 
{ynumeric[1]=4.7
color='black'
}

steps=10
sddeltat=(1/steps)^0.5

for(step in 1:(steps*16))
{ 
  x=ynumeric[step]
  if(step<100)
  deltax=-(x-x1)*(2-1/(1+100/step))
  deltax=deltax*1/steps
  deltaguassian=rnorm(1,mean=0,sd=sddeltat*(B))
  ynumeric[step+1]=ynumeric[step]+deltax*1/steps+deltaguassian
}
xnumeric=seq(0,length(ynumeric)-1,1)
par(mar=c(3,3,1,1))
plot(xnumeric/steps,ynumeric,type='l',col=color,xlim=c(0,16),ylim=c(4.5,6),tck=0.03,mgp=c(0,0.5,0),las=1,xlab="", ylab="", main="",xaxt='n',yaxt='n')
par(new=T)

axis(side=1,las=1,at=c(0,5,10,15),mgp=c(0,0.5,0),tck=0.03,las=1,labels=c("0","5","10","15"),cex.axis=1)
axis(side=2,las=1,at=c(4.5,5,5.5,6),mgp=c(0,0.5,0),tck=0.03,las=1,labels=c("4.5","5","5.5","6"),cex.axis=1)

mtext("Days",side=1,line=1.5,cex=1)
mtext("scEntropy",side=2,line=2.0,cex=1)

#xun zhao寻找 gui xian圭贤

scentropysimu=ynumeric[seq(1,10,1)]
for(i in 1:10)
{
  scentropysimu[i]=which(abs(distancestat-scentropysimu[i])==min(abs(distancestat-scentropysimu[i])))
}
scentropysimutotal=data.frame(scentropysimu)
for(i in 1:15)
{
  for(j in 1:10)
  {scentropysimu[j]=which(abs(distancestat-ynumeric[i*10+j])==min(abs(distancestat-ynumeric[i*10+j])))}
  scentropysimutotal=data.frame(scentropysimutotal,scentropysimu)
}