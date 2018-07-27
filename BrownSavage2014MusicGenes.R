#The following code was used for the analyses reported in:
#Brown, S., Savage, P. E., Ko, A. M.-S., Stoneking, M., Ko, Y.-C., Loo, J.-H., & Trejaut, J. A. (2014). Correlations in the population structure of music, genes and language. Proceedings of the Royal Society B: Biological Sciences, 281(1774), 20132072.

#The following packages must be installed/loaded
library(ade4)
library(vegan)

#set working drive
setwd("/Users/pesavage/Documents/Research/Papers/Published/Brown et al (2014) Proc R Soc B/Code for Github")


#Correlations with CantoCore +Cantometrics averaged distance matrices (PEARSON r)
taiwan.ling9.dist<- as.dist(read.csv("a400patristic9.csv",header=TRUE,row.names=1))
taiwangeo.dist<-as.dist(read.table("9_geo_dist_converted_km.txt",row.names=1,header=TRUE))
taiwan259CCCMFst.dist<-as.dist(read.csv("Taiwan9DistCCCMCombined.csv",header=TRUE,row.names=1))
taiwanFst.dist<-as.dist(read.csv("Taiwan9DistNewJeanAlbertCombined.csv",header=TRUE,row.names=1))

mantel(taiwanFst.dist,taiwan259CCCMFst.dist,permutations=10000)
#Mantel statistic r: 0.4169 
 #     Significance: 0.014899
mantel(taiwan.ling9.dist,taiwan259CCCMFst.dist,permutations=10000)
#Mantel statistic r: 0.3411 
#      Significance: 0.12909
 mantel(taiwangeo.dist,taiwan259CCCMFst.dist,permutations=10000)
#Mantel statistic r: 0.1742 
#      Significance: 0.23398 

mantel.partial(taiwan259CCCMFst.dist,taiwanFst.dist, taiwan.ling9.dist,permutations=10000)
#Mantel statistic r: 0.3013 
 #     Significance: 0.046095

mantel.partial(taiwan259CCCMFst.dist,taiwanFst.dist, taiwangeo.dist,permutations=10000)
#Mantel statistic r: 0.3854 
#      Significance: 0.031597 
mantel.partial(taiwan259CCCMFst.dist,taiwan.ling9.dist,taiwangeo.dist,permutations=10000)
#Mantel statistic r: 0.2989 
 #     Significance: 0.14879 

mantel.partial(taiwanFst.dist,taiwan.ling9.dist, taiwan259CCCMFst.dist,permutations=10000)
#Mantel statistic r: 0.4257 
 #     Significance: 0.016998 

mantel.partial(taiwanFst.dist,taiwan.ling9.dist, taiwangeo.dist,permutations=10000)

# Mantel statistic r: 0.3289 
#      Significance: 0.057194

#Music-gene-language regressions:
plot(taiwan259CCCMFst.dist,taiwanFst.dist,pch=19,xlim=c(0,0.12),ylim=c(0,0.2),cex.axis=2,cex=2)
music.gene.lm<-lm(taiwanFst.dist~taiwan259CCCMFst.dist)
abline(music.gene.lm)
plot(taiwan259CCCMFst.dist,taiwan.ling9.dist,pch=19,xlim=c(0,0.12),ylim=c(0,0.1), cex.axis=2,cex=2)
music.lang.lm<-lm(taiwan.ling9.dist~taiwan259CCCMFst.dist)
abline(music.lang.lm)
plot(taiwanFst.dist,taiwan.ling9.dist,pch=19,xlim=c(0,0.2),ylim=c(0,0.1), cex.axis=2,cex=2)
gene.lang.lm<-lm(taiwan.ling9.dist~taiwanFst.dist)
abline(gene.lang.lm)


#Creating musical distance matrices

#Distances among songs
#First, check spreadsheet for errors, save as .csv spreadsheet with first row=CantoCoreID and rows 2-44= CC1-26+CMPerformance+2 Instrumental
#import spreadsheet
taiwan<-as.matrix(read.csv("Taiwan259SongFinalSampleCCCMInstHemitonicityNominal.csv",row.names=1,header=TRUE))

#add distance matrix algorithm
ordinal.fn<-function (x,y) {
if(is.na(x)|is.na(y))"NA" else
abs(x-y)}


weightedv6.dist<-function(d,ord,nom) {

d[d==""]<-NA

nominal<-cbind(d[,nom])

nominal.fields<-vector("list",length=length(nominal[1,]))

for (i in 1:length(nominal[1,])) {nominal.fields[[i]]<-matrix(nrow=length(nominal[,1]),ncol=13)}

for (i in 1:length(nominal[1,])) {rownames(nominal.fields[[i]])<-rownames(nominal)}

for (i in 1:length(nominal[1,])) {colnames(nominal.fields[[i]])<-c("a","b","c","d","e","f","g","h","i","j","k","l","m")}

for (j in 1:13){
for (k in 1:length(nominal[,1])){
for (i in 1:length(nominal[1,])){nominal.fields[[i]][k,j]<-(if(is.na(nominal[k,i]))"NA" else if(substr(nominal[k,i],1,1)==colnames(nominal.fields[[i]])[j] | substr(nominal[k,i],2,2)==colnames(nominal.fields[[i]])[j] |substr(nominal[k,i],3,3)==colnames(nominal.fields[[i]])[j] |substr(nominal[k,i],4,4)==colnames(nominal.fields[[i]])[j] |substr(nominal[k,i],5,5)==colnames(nominal.fields[[i]])[j] |substr(nominal[k,i],6,6)==colnames(nominal.fields[[i]])[j] |substr(nominal[k,i],7,7)==colnames(nominal.fields[[i]])[j] |substr(nominal[k,i],8,8)==colnames(nominal.fields[[i]])[j] |substr(nominal[k,i],9,9)==colnames(nominal.fields[[i]])[j] |substr(nominal[k,i],10,10)==colnames(nominal.fields[[i]])[j] |substr(nominal[k,i],11,11)==colnames(nominal.fields[[i]])[j] |substr(nominal[k,i],12,12)==colnames(nominal.fields[[i]])[j] |substr(nominal[k,i],13,13)==colnames(nominal.fields[[i]])[j])1 else 0) }}}

suppressWarnings(for (i in 1:length(nominal[1,])){storage.mode(nominal.fields[[i]])<-"numeric"})

nominal.dist<-vector("list",length=length(nominal[1,]))

for (i in 1:length(nominal[1,])) {nominal.dist[[i]]<-matrix(nrow=length(nominal[,1]),ncol=length(nominal[,1]))}

nominal.result<-matrix(nrow=length(d[,1]),ncol=length(d[,1]),c(rep(0,(length(d[,1])*length(d[,1])))))

ordinal<-cbind(d[,ord])
suppressWarnings(storage.mode(ordinal)<-"numeric")

ordinal.dist<-vector("list",length=length(ordinal[1,]))

for (i in 1:length(ordinal[1,])) {ordinal.dist[[i]]<-matrix(nrow=length(ordinal[,1]),ncol=length(ordinal[,1]))}

ordinal.result<-matrix(nrow=length(d[,1]),ncol=length(d[,1]),c(rep(0,(length(d[,1])*length(d[,1])))))

result<-matrix(nrow=length(d[,1]),ncol=length(d[,1]))



for (k in 1:length(nominal[,1])){
for (j in 1:length(nominal[,1])){
for (i in 1:length(nominal[1,])){nominal.dist[[i]][k,j]<-if(is.na(nominal.fields[[i]][k,1])|is.na(nominal.fields[[i]][j,1]))"NA" else ((if(nominal.fields[[i]][k,1]==nominal.fields[[i]][j,1])0 else 1)+(if(nominal.fields[[i]][k,2]==nominal.fields[[i]][j,2])0 else 1)+(if(nominal.fields[[i]][k,3]==nominal.fields[[i]][j,3])0 else 1)+(if(nominal.fields[[i]][k,4]==nominal.fields[[i]][j,4])0 else 1)+(if(nominal.fields[[i]][k,5]==nominal.fields[[i]][j,5])0 else 1)+(if(nominal.fields[[i]][k,6]==nominal.fields[[i]][j,6])0 else 1)+(if(nominal.fields[[i]][k,7]==nominal.fields[[i]][j,7])0 else 1)+(if(nominal.fields[[i]][k,8]==nominal.fields[[i]][j,8])0 else 1)+(if(nominal.fields[[i]][k,9]==nominal.fields[[i]][j,9])0 else 1)+(if(nominal.fields[[i]][k,10]==nominal.fields[[i]][j,10])0 else 1)+(if(nominal.fields[[i]][k,11]==nominal.fields[[i]][j,11])0 else 1)+(if(nominal.fields[[i]][k,12]==nominal.fields[[i]][j,12])0 else 1))/((if(nominal.fields[[i]][k,1]==1 | nominal.fields[[i]][[j,1]]==1)1 else 0)+(if(nominal.fields[[i]][k,2]==1 | nominal.fields[[i]][[j,2]]==1)1 else 0)+(if(nominal.fields[[i]][k,3]==1 | nominal.fields[[i]][[j,3]]==1)1 else 0)+(if(nominal.fields[[i]][k,4]==1 | nominal.fields[[i]][[j,4]]==1)1 else 0)+(if(nominal.fields[[i]][k,5]==1 | nominal.fields[[i]][[j,5]]==1)1 else 0)+(if(nominal.fields[[i]][k,6]==1 | nominal.fields[[i]][[j,6]]==1)1 else 0)+(if(nominal.fields[[i]][k,7]==1 | nominal.fields[[i]][[j,7]]==1)1 else 0)+(if(nominal.fields[[i]][k,8]==1 | nominal.fields[[i]][[j,8]]==1)1 else 0)+(if(nominal.fields[[i]][k,9]==1 | nominal.fields[[i]][[j,9]]==1)1 else 0)+(if(nominal.fields[[i]][k,10]==1 | nominal.fields[[i]][[j,10]]==1)1 else 0)+(if(nominal.fields[[i]][k,11]==1 | nominal.fields[[i]][[j,11]]==1)1 else 0)+(if(nominal.fields[[i]][k,12]==1 | nominal.fields[[i]][[j,12]]==1)1 else 0))}}}
suppressWarnings(for (i in 1:length(nominal[1,])){storage.mode(nominal.dist[[i]])<-"numeric"})

vnom<-vector(mode="numeric",length=length(nominal[1,]))
for (k in 1:length(nominal[,1])){
for (j in 1:length(nominal[,1])){
for (i in 1:length(nominal[1,])){
  vnom[i]<-nominal.dist[[i]][k,j]
  }
  nominal.result[k,j]<-mean(suppressWarnings(as.numeric(vnom)),na.rm=TRUE)
  }}

for (k in 1:length(ordinal[,1])){
for (j in 1:length(ordinal[,1])){
for (i in 1:length(ordinal[1,])){ordinal.dist[[i]][k,j]<-ordinal.fn(x=ordinal[k,i],y=ordinal[j,i]) }}}
suppressWarnings(for (i in 1:length(ordinal[1,])){storage.mode(ordinal.dist[[i]])<-"numeric"})

vord<-vector(mode="numeric",length=length(ordinal[1,]))
for (k in 1:length(ordinal[,1])){
for (j in 1:length(ordinal[,1])){
for (i in 1:length(ordinal[1,])){
  vord[i]<-ordinal.dist[[i]][k,j]
  }
  ordinal.result[k,j]<-mean(suppressWarnings(as.numeric(vord)),na.rm=TRUE)
  }}


for (k in 1:length(nominal[,1])){
for (j in 1:length(nominal[,1])){
result[k,j]<-if(is.na(nominal.result[k,j]))ordinal.result[k,j] else if(is.na(ordinal.result[k,j])) nominal.result[k,j] else (ordinal.result[k,j]*length(ord) + nominal.result[k,j]*length(nom))/(length(ord)+length(nom))
}}

row.names(result)<-row.names(d)
colnames(result)<-row.names(d)
as.dist(result)
}
#w/out instruments
taiwan.dist<-weightedv6.dist(taiwan,c(5:7,10,12:13,15:17,19,21:23,26:41),c(1:4,8:9,11,14,18,20,24:25))
#NB: This treats hemitonicity as nominal, not ordinal, because that’s how it used to be coded
taiwan.euc.dist<-lingoes(taiwan.dist)
taiwan.sq.euc.dist<-(taiwan.euc.dist)^2
taiwan.sq.euc.frame<-as.data.frame(as.matrix(taiwan.sq.euc.dist))
taiwan.sq.euc.frame[upper.tri(taiwan.sq.euc.frame)] <- NA
write.table(taiwan.sq.euc.frame,"TaiwanCCEucDist.txt",na="", row.names=FALSE)

#Distances among populations (AMOVA)

#COPY THIS DISTANCE MATRIX INTO THE DISTNACE MATRIX SECTION OF THE ARLEQUIN .ARP FILE AND RUN AMOVA USING PROVIDED DISTANCE MATRIX
#After saving this .arp (or other files, e.g., .txt) MAKE SURE TO OPEN IT ONCE WITH EXCEL IN WINDOWS AND RE-SAVE IT AS .TXT, OTHERWISE THE FORMATTING GETS ALL MESSED UP!!!
# THEN COPY RESULTING FST TABLES AND SAVE AS FOLLOWING .CSV FILES

#Inter-rater reliability (New w Sakai Yamashita Shimozaki) [NB: This is not the exact code/file used in Brown et al., 2014, but a more streamlined version used in other later analyses]
library(irr)
cc.ep.mat<-as.matrix(read.csv("ConsensusTapeCantoCoreCantometricsEmilyPat41.csv",header=TRUE,row.names=1))
#(manually edit out R code pertaining to empty columns for the Lomax comparisons)
#Kappa (ignoring NAs)
cc.sakai.savage.kappa<-cbind(kappa2(cc.ep.mat[,c(1,49)]), kappa2(cc.ep.mat[,c(2,50)],weight="squared"),kappa2(cc.ep.mat[,c(3,51)]),
kappa2(cc.ep.mat[,c(4,52)],weight="squared"),
kappa2(cc.ep.mat[,c(5,53)],weight="squared"),
kappa2(cc.ep.mat[,c(6,54)],weight="squared"),
kappa2(cc.ep.mat[,c(7,55)],weight="squared"),
kappa2(cc.ep.mat[,c(8,56)]),
kappa2(cc.ep.mat[,c(9,57)]),
kappa2(cc.ep.mat[,c(10,58)],weight="squared"),
kappa2(cc.ep.mat[,c(11,59)]),
kappa2(cc.ep.mat[,c(12,60)],weight="squared"),
kappa2(cc.ep.mat[,c(13,61)],weight="squared"),
kappa2(cc.ep.mat[,c(14,62)]),
kappa2(cc.ep.mat[,c(15,63)],weight="squared"),
kappa2(cc.ep.mat[,c(16,64)],weight="squared"),
kappa2(cc.ep.mat[,c(17,65)],weight="squared"),
kappa2(cc.ep.mat[,c(18,66)]),
kappa2(cc.ep.mat[,c(19,67)],weight="squared"),
kappa2(cc.ep.mat[,c(20,68)]),
kappa2(cc.ep.mat[,c(21,69)],weight="squared"),
kappa2(cc.ep.mat[,c(22,70)],weight="squared"),
kappa2(cc.ep.mat[,c(23,71)],weight="squared"),
kappa2(cc.ep.mat[,c(24,72)]),
kappa2(cc.ep.mat[,c(25,73)]),
kappa2(cc.ep.mat[,c(26,74)],weight="squared"),
kappa2(cc.ep.mat[,c(27,75)],weight="squared"),
kappa2(cc.ep.mat[,c(28,76)],weight="squared"),
kappa2(cc.ep.mat[,c(29,77)],weight="squared"),
kappa2(cc.ep.mat[,c(30,78)],weight="squared"),
kappa2(cc.ep.mat[,c(31,79)],weight="squared"),
kappa2(cc.ep.mat[,c(32,80)],weight="squared"),
kappa2(cc.ep.mat[,c(33,81)],weight="squared"),
kappa2(cc.ep.mat[,c(34,82)],weight="squared"),
kappa2(cc.ep.mat[,c(35,83)],weight="squared"),
kappa2(cc.ep.mat[,c(36,84)],weight="squared"),
kappa2(cc.ep.mat[,c(37,85)],weight="squared"),
kappa2(cc.ep.mat[,c(38,86)],weight="squared"),
kappa2(cc.ep.mat[,c(39,87)],weight="squared"),
kappa2(cc.ep.mat[,c(40,88)],weight="squared"),
kappa2(cc.ep.mat[,c(41,89)],weight="squared"))

write.csv(cc.sakai.savage.kappa,"cc.sakai.savage.kappa.csv")

#%Agreement (including NAs)

#First, replace NAs with 0
cc.ep.mat[is.na(cc.ep.mat)] <- 0

#Then check % agreement

cc.sakai.savage.kappa<-cbind(agree(cc.ep.mat[,c(1,49)]), agree(cc.ep.mat[,c(2,50)]),agree(cc.ep.mat[,c(3,51)]),
agree(cc.ep.mat[,c(4,52)]),
agree(cc.ep.mat[,c(5,53)]),
agree(cc.ep.mat[,c(6,54)]),
agree(cc.ep.mat[,c(7,55)]),
agree(cc.ep.mat[,c(8,56)]),
agree(cc.ep.mat[,c(9,57)]),
agree(cc.ep.mat[,c(10,58)]),
agree(cc.ep.mat[,c(11,59)]),
agree(cc.ep.mat[,c(12,60)]),
agree(cc.ep.mat[,c(13,61)]),
agree(cc.ep.mat[,c(14,62)]),
agree(cc.ep.mat[,c(15,63)]),
agree(cc.ep.mat[,c(16,64)]),
agree(cc.ep.mat[,c(17,65)]),
agree(cc.ep.mat[,c(18,66)]),
agree(cc.ep.mat[,c(19,67)]),
agree(cc.ep.mat[,c(20,68)]),
agree(cc.ep.mat[,c(21,69)]),
agree(cc.ep.mat[,c(22,70)]),
agree(cc.ep.mat[,c(23,71)]),
agree(cc.ep.mat[,c(24,72)]),
agree(cc.ep.mat[,c(25,73)]),
agree(cc.ep.mat[,c(26,74)]),
agree(cc.ep.mat[,c(27,75)]),
agree(cc.ep.mat[,c(28,76)]),
agree(cc.ep.mat[,c(29,77)]),
agree(cc.ep.mat[,c(30,78)]),
agree(cc.ep.mat[,c(31,79)]),
agree(cc.ep.mat[,c(32,80)]),
agree(cc.ep.mat[,c(33,81)]),
agree(cc.ep.mat[,c(34,82)]),
agree(cc.ep.mat[,c(35,83)]),
agree(cc.ep.mat[,c(36,84)]),
agree(cc.ep.mat[,c(37,85)]),
agree(cc.ep.mat[,c(38,86)]),
agree(cc.ep.mat[,c(39,87)]),
agree(cc.ep.mat[,c(40,88)]),
agree(cc.ep.mat[,c(41,89)]))

write.csv(cc.sakai.savage.kappa,"cc.sakai.savage.kappa.csv")
