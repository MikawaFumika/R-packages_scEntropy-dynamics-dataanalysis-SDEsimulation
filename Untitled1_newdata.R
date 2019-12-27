setwd("/Users/yeyusong/Desktop/在投/entropydynamics/shuju")
mouse <- read.table("GSE103221_normalized_counts.csv",fill=TRUE,header=T,sep=',')
templatemouse <- read.table("mart_export-1.txt",fill=TRUE,header=T,sep=',')

#templatemouse预处理


#exon length template1 intersect
geneid <- mouse[,1]
geneid <- (as.character(geneid))
exontemplate2 <- matrix(geneid)
exonlength2 <- vector()
for(i in 2956:20273)
{
  generaw=grep(geneid[i],templatemouse[,7])
  if(length(generaw)>0)
  {
  start=as.numeric(as.character((templatemouse[generaw,1])))
  end=as.numeric(as.character(templatemouse[generaw,2]))
  genelengthtotal=seq(start[1],end[1],1)
  if(length(generaw)>1)
  {
    for(j in 2:(length(generaw)))
      genelengthtotal=union(genelengthtotal,seq(start[j],end[j],1))
  }
  exonlength2[i]=length(genelengthtotal)
  }
  else
  {
  exonlength2[i]=0
  }
  print(i)
}
exontemplate2=data.frame(exontemplate2,exonlength2)
write.table(exontemplate2,"exontemplate2",quote = FALSE,row.names = FALSE,col.names = FALSE,sep="\t")

setwd("/Users/yeyusong/Desktop/在投/entropydynamics/shuju")
exontemplate2 <- read.table("exontemplate2",fill=TRUE,header=F,sep='\t')
mouse <- read.table("GSE103221_normalized_counts.csv",fill=TRUE,header=T,sep=',')

#process tpm
exonlength2=exontemplate2[,2]
mouse=mouse[which(exonlength2!=0),]
exonlength2=exonlength2[which(exonlength2!=0)]

rpkmdatamouse=as.numeric(as.character(mouse[,2]))/exonlength2
rpkmdatamouse=rpkmdatamouse/sum(rpkmdatamouse)*1000000
rpkmdatamouse=log((rpkmdatamouse+1),2)
for(j in 3:913)
{  
  rpkm=as.numeric(as.character(mouse[,j]))/exonlength2
  rpkm=rpkm/sum(rpkm)*1000000
  rpkm=log((rpkm+1),2)
  rpkmdatamouse=data.frame(rpkmdatamouse,rpkm)
  print(j)
}



#define days
cells=vector()
cellaaa=0
for(cell in 1:912)
{
  if(substring(colnames(mouse)[cell+1],1,3)=='mef')
  {cellaaa=cellaaa+1
  cells[cellaaa]=cell}
}  
Dayscell=list(cells)


cells=vector()
cellaaa=0
for(cell in 1:912)
{
  if(substring(colnames(mouse)[cell+1],1,6)=='osk_d0')
  {cellaaa=cellaaa+1
  cells[cellaaa]=cell}
}  
cells=list(cells)
Dayscell=c(Dayscell,cells)

cells=vector()
cellaaa=0
for(cell in 1:912)
{
  if(substring(colnames(mouse)[cell+1],1,6)=='osk_d1')
  {cellaaa=cellaaa+1
  cells[cellaaa]=cell}
}  
cells=list(cells)
Dayscell=c(Dayscell,cells)

cells=vector()
cellaaa=0
for(cell in 1:912)
{
  if(substring(colnames(mouse)[cell+1],1,6)=='osk_d2')
  {cellaaa=cellaaa+1
  cells[cellaaa]=cell}
}  
cells=list(cells)
Dayscell=c(Dayscell,cells)

cells=vector()
cellaaa=0
for(cell in 1:912)
{
  if(substring(colnames(mouse)[cell+1],1,6)=='osk_d3')
  {cellaaa=cellaaa+1
  cells[cellaaa]=cell}
}  
cells=list(cells)
Dayscell=c(Dayscell,cells)


cells=vector()
cellaaa=0
for(cell in 1:912)
{
  if(substring(colnames(mouse)[cell+1],1,6)=='osk_d5')
  {cellaaa=cellaaa+1
  cells[cellaaa]=cell}
}  
cells=list(cells)
Dayscell=c(Dayscell,cells)

cells=vector()
cellaaa=0
for(cell in 1:912)
{
  if(substring(colnames(mouse)[cell+1],1,6)=='osk_d7')
  {cellaaa=cellaaa+1
  cells[cellaaa]=cell}
}  
cells=list(cells)
Dayscell=c(Dayscell,cells)

cells=vector()
cellaaa=0
for(cell in 1:912)
{
  if(substring(colnames(mouse)[cell+1],1,6)=='osk_d8')
  {cellaaa=cellaaa+1
  cells[cellaaa]=cell}
}  
cells=list(cells)
Dayscell=c(Dayscell,cells)

cells=vector()
cellaaa=0
for(cell in 1:912)
{
  if(substring(colnames(mouse)[cell+1],1,3)=='ips')
  {cellaaa=cellaaa+1
  cells[cellaaa]=cell}
}  
cells=list(cells)
Dayscell=c(Dayscell,cells)

cells=vector()
cellaaa=0
for(cell in 1:912)
{
  if(substring(colnames(mouse)[cell+1],1,3)=='esc')
  {cellaaa=cellaaa+1
  cells[cellaaa]=cell}
}  
cells=list(cells)
Dayscell=c(Dayscell,cells)

#difine refer
cells=unlist(Dayscell[10])
w=vector()
for(i in 1:912)
  w[i]=0
for(i in cells)
  w[i]=1
w=w/length(cells)
#define reference set
referenceset=as.numeric(as.character(rpkmdatamouse[,1]))*w[1]
for(i in 2:912)
{
  referenceset=referenceset+as.numeric(as.character(rpkmdatamouse[,i]))*w[i]
}
hist(referenceset)


#distance calculation
distancestat=vector()
for(i in 1:912)
{
  y=as.numeric(as.character(rpkmdatamouse[,i]))-referenceset
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

#plot distance all
par(mar=c(3,3,1,1))
colors <- colorRampPalette(c(rgb(238/255,44/255,44/255), rgb(0/255,144/255,255/255)))(10)
for(i in 0:9)
{
  cells=unlist(Dayscell[i+1])
  plot((runif(length(distancestat[cells][which(distancestat[cells]>4.7)]))/2+0.75+i),distancestat[cells][which(distancestat[cells]>4.7)],xlim=c(0,11),ylim=c(4.5,6),pch=19,cex=0.5,col=colors[i+1],xlab="", ylab="", main="",xaxt="n",yaxt="n",mgp=c(0,0.5,0),tck=0.03)
  par(new=TRUE)
}
axis(side=1,las=1,at=c(1,2,3,4,5,6,7,8,9,10),mgp=c(0,0.5,0),tck=0.03,las=1,labels=c("MEF","0","1","2","3","5","7","8","iPSC","ESC"),cex.axis=0.8)
axis(side=2,las=1,at=c(4.5,5,5.5,6),mgp=c(0,0.5,0),tck=0.03,las=1,labels=c("4.5","5","5.5","6"),cex.axis=1)
mtext("Cell group (Days)",side=1,line=1.5,cex=1)
mtext("scEntropy",side=2,line=1.7,cex=1)


#jun zhi均值fang cha
meandistancea=vector()
vardistancea=vector()
for(i in 1:10)
{  
  meandistancea[i]=mean(distancestat[unlist(Dayscell[i])][which(distancestat[unlist(Dayscell[i])]>4.7)])
  vardistancea[i]=var(distancestat[unlist(Dayscell[i])][which(distancestat[unlist(Dayscell[i])]>4.7)])
}

par(mar=c(3,3,1,1))
colors <- colorRampPalette(c(rgb(238/255,44/255,44/255), rgb(0/255,144/255,255/255)))(10)
plot(seq(0,9),meandistancea,type="p",pch=19,tck=0.03,cex=0.8,las=1,xlab="",col=colors[rank(1:10)], ylab="", main="",xlim=c(0,9),ylim=c(4.5,6),xaxt="n",yaxt="n",bty="l")
mtext("Cell group (Days)",side=1,line=0.5,cex=1)
mtext("Average scEntropy",side=2,line=1.5,cex=1)
axis(side=1,las=1,at=seq(0,9,1),mgp=c(0,0.5,0),tck=0.03,las=1,labels=c("","","","","","","","","",""),cex.axis=0.4)
axis(side=2,las=1,at=c(4.5,5,5.5,6),mgp=c(0,0.5,0),tck=0.03,las=1,labels=c("4.5","5","5.5","6"),cex.axis=1)

par(mar=c(3,3.5,1,1))
colors <- colorRampPalette(c(rgb(238/255,44/255,44/255), rgb(0/255,144/255,255/255)))(10)
plot(seq(0,9),vardistancea,type="p",pch=19,tck=0.03,cex=0.8,las=1,xlab="",col=colors[rank(1:10)], ylab="", main="",xlim=c(0,9),ylim=c(0,0.12),xaxt="n",yaxt="n",bty="l")
mtext("Cell group (Days)",side=1,line=0.5,cex=1)
mtext("Variance of scEntropy",side=2,line=2,cex=1)
axis(side=1,las=1,at=seq(0,9,1),mgp=c(0,0.5,0),tck=0.03,las=1,labels=c("","","","","","","","","",""),cex.axis=0.4)
axis(side=2,las=1,at=c(0,0.05,0.1),mgp=c(0,0.5,0),tck=0.03,las=1,labels=c(0,0.05,0.1),cex.axis=1)



#forceatlas


axismatrix=matrix(nrow=912,ncol=2)

for(i in 1:912)
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

iteration=100000

for(iter in 1:iteration)
{ 
  samplex=sample(seq(1,912,1),1)
  sampley=sample(seq(1,912,1),1)
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



colors <- colorRampPalette(c(rgb(0/255,144/255,255/255), rgb(238/255,44/255,44/255)))(10)

for(i in seq(0,9,1))
{
  i=i+1
  cells=unlist(Dayscell[i+1])
  plot(axismatrix[cells,1],axismatrix[cells,2],cex=0.2,xlim=c(min(axismatrix[,1]),max(axismatrix[,1])),ylim=c(min(axismatrix[,2]),max(axismatrix[,2])),col=colors[i+1],xlab="",ylab="")
  par(new=TRUE)
  #Sys.sleep(1)
}
#pingjun guixian
for(i in seq(0,9,1))
{
  i=i+1
  cells=unlist(Dayscell[i+1])
  cells1=unlist(Dayscell[i+2])
  plot(c(mean(axismatrix[cells,1]),mean(axismatrix[cells1,1])),c(mean(axismatrix[cells,2]),mean(axismatrix[cells1,2])),type='l',lwd=2.0,xlim=c(min(axismatrix[,1]),max(axismatrix[,1])),ylim=c(min(axismatrix[,2]),max(axismatrix[,2])),col='black',xlab="",ylab="")
  par(new=TRUE)
  Sys.sleep(1)
}
