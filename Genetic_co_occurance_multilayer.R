#set up data for multilayer network

####Data loading ########


#load data SK

setwd("C:/Users/Documents/EC PDF/CMRnet")
gn<-read.csv("SK_15alleles_relatedness.csv", header=T)
sn_sri<-read.csv("SK_allyears_edgelist_SRI.csv", header=T)
all_nodes<-read.csv("SK_nodes.csv", header=T)


############################################



##### set up multilayer file #############

setwd("C:/Users/Documents/EC PDF/muxViz-master/data/SK")

sn_sri$pairid<-rep(NA,nrow(sn_sri))
for (i in (1:nrow(sn_sri))){
  sn_sri$pairid[i]<-paste(min(c(sn_sri$id1[i], sn_sri$id2[i])), "-", max(c(sn_sri$id1[i], sn_sri$id2[i])), sep="")  
}

sn_sri$weight2<-sn_sri$weight/5

#threshold genetic network 

gn2<-gn[which (gn$Weight>=0.25),]
gn3<-gn[which (gn$Weight>=0.4),]


sn2<-cbind.data.frame(sn_sri$id1, sn_sri$id2, sn_sri$weight2)
names(sn2)<-c("id1", "id2","weight")

all_ids<-unique(all_nodes$id)


write.table(gn2, "gn25.edges", row.names=F, col.names=F)
write.table(gn3, "gn4.edges", row.names=F, col.names=F)
write.table(gn, "gnall.edges", row.names=F, col.names=F, )


sn<-sn_sri[,1:3]
sn_sri2<-sn_sri[,c(1:2,6)]

write.table(sn, "sn_num.edges", row.names=F, col.names=F)
write.table(sn_sri2, "sn_sri.edges", row.names=F, col.names=F)
write.table(sn2, "sn_stand.edges", row.names=F, col.names=F)

num_id<-1:length(all_ids)


#SK layout
SKlayout<-cbind.data.frame(num_id, all_ids)
names(SKlayout)<-c("nodeID", "nodeLabel")
write.table(SKlayout, "SK_layout.txt", row.names=F)


#check correlation between social association and genetic association 

sn_sri$pairid<-rep(NA,nrow(sn_sri))
for (i in (1:nrow(sn_sri))){
  sn_sri$pairid[i]<-paste(min(c(sn_sri$id1[i], sn_sri$id2[i])), "-", max(c(sn_sri$id1[i], sn_sri$id2[i])), sep="")  
}


gn2$pairid<-rep(NA,nrow(gn2))
for (i in (1:nrow(gn2))){
  gn2$pairid[i]<-paste(min(c(gn2$Target[i], gn2$Source[i])), "-", max(c(gn2$Target[i], gn2$Source[i])), sep="")  
}

all_dat<-merge(gn2, sn_sri, by="pairid", all.y=T, all.x=T)

all_dat[which (is.na(all_dat$SRI)),]$SRI<-0
all_dat[which (is.na(all_dat$Weight)),]$Weight<-0


#remove INF values
check<-all_dat[which (all_dat$weight>1),]
all_dat2<-all_dat[!all_dat$pairid %in% check$pairid,]
all_dat2<-all_dat
plot(all_dat2$SRI~all_dat2$Weight)

cor.test(all_dat2$SRI, all_dat2$Weight)

library(lme4)
mod1<-lm(log(weight+1)~log(Weight+1), data=all_dat2)




#standardize sn to 0-1
#sn$weight2<-sn$weight/(max(sn$weight))



#add missing ind to end of edge list
allnew<-data.frame(matrix(0,nrow=0,ncol=3))
for (i in unique(check)){
  newline<-cbind(i, i, 0)
  allnew<-rbind(allnew, newline)
}
colnames(allnew)<-c("Source", "Target", "Weight")
gn_new<-rbind(gn3, allnew)



#add missing ind to end of edge list
allnew2<-data.frame(matrix(0,nrow=0,ncol=3))
for (i in unique(check2)){
  newline<-cbind(i, i, 0)
  allnew2<-rbind(allnew2, newline)
}
colnames(allnew2)<-c("id1", "id2", "SRI")
sn_new<-rbind(sn_sri2, allnew2)


#order edgelist
sn_new2<-sn_new[order(sn_new[1:2]),]
gn_new2<-gn_new[order(gn_new[1:2]),]

sn_new3<-na.omit(sn_new2)


gn_new3<-na.omit(gn_new2)

#make networks

library(igraph)


net_sn <- graph.data.frame(sn_new3, directed=F)
net_sn2<- get.adjacency(net_sn, type="both", attr='SRI')


getorder<-colnames(net_sn2) 
getorder2<-cbind.data.frame(getorder, c(getorder[2:length(getorder)],"1"))
colnames(getorder2)<-c("Source", "Target")
getorder2$Weight<-rep(0, length(getorder2$Source))
alledges3<-rbind(getorder2, gn_new3)

net_gn <- graph.data.frame(alledges3, directed=F)
net_gn2<-get.adjacency(net_gn, type="both", attr='Weight')

#write matrices
sn_net<-as.matrix(net_sn2)  
gn_net<-as.matrix(net_gn2)

write.csv(sn_net, "sn_net.csv", row.names=T)
write.csv(gn_net, "gn_net.csv", row.names=T)




#run muxviz 

setwd("C:/Users/Documents/EC PDF/muxViz-master")
source('muxVizGUI.R')


