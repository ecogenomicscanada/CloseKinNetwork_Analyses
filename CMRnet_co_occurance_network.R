
library(lubridate)
library(CMRnet)
library(glmmTMB)
library(igraph)

#set directory
setwd("C:/Users/Documents/EC PDF/CMRnet")
CMRdata<-read.csv("CaptureHistory.csv", header=T)



#make unique ID -- for 0 PCID value is unique for all others cluster id is unique
CMRdata$id<-CMRdata$Cluster
CMRdata[which (CMRdata$Cluster==0),]$id<-CMRdata[which (CMRdata$Cluster==0),]$PCID
head(CMRdata)
length(unique(CMRdata$id))

summary(CMRdata)

#add number for year for permutations
CMRdata$tp<-rep(NA, nrow(CMRdata))
CMRdata[which (CMRdata$Year=="2017"),]$tp<-1
CMRdata[which (CMRdata$Year=="2018"),]$tp<-2
CMRdata[which (CMRdata$Year=="2019"),]$tp<-3



#get x and y data
CMRdata$UTM_Easting<-as.integer(CMRdata$utm_easting)
CMRdata$UTM_Northing<-as.integer(CMRdata$utm_northing)

#make unique site name (combine x and y to get a unique site label for each sampling coordinates)
CMRdata$loc2<-paste(CMRdata$UTM_Easting,CMRdata$UTM_Northing, sep="_")

#make single date column 
CMRdata$date<-paste(CMRdata$Year, "-", CMRdata$Month, "-", CMRdata$Day, sep="")
CMRdata$date<-as.Date(CMRdata$date)

CMRdata3<-data.frame(CMRdata$id,CMRdata$loc2,CMRdata$UTM_Easting,CMRdata$UTM_Northing,CMRdata$date)
names(CMRdata3)<-c("id","loc","x","y","date")

#check and remove duplicate values
check<-duplicated(CMRdata3)

CMRdata4<-unique(CMRdata3)

CMRdata5<-CMRdata4[,c(1,3,4)]
write.csv(CMRdata5, "SK_nodes.csv", row.names=F)

#####make social networks######

#set up CMRnet for data - need to change time period for each sample data 
#make one network list for all data


#years of data collection 
Yearlist<-sort(unique(CMRdata$Year))
minyear<-min(Yearlist)
maxyear<-max(Yearlist)
netdats<-list()
attlist<-list()

#make each network -- 
for (i in Yearlist){
  mindate<-paste0(i,"-01-01")
  maxdate<-paste0(i+1,"-01-01")
  intwindow<-1
  netwindow<-12
  overlap<-0
  spacewindow<-0
  
  
  netdat1<-DynamicNetCreate(data=CMRdata4,
                            intwindow=intwindow,
                            mindate=mindate,
                            maxdate=maxdate, 
                            netwindow=netwindow, 
                            overlap=overlap, 
                            spacewindow=spacewindow)
  
  
  ##Convert social networks into a list of igraph networks
  
  net<-graph.adjacency(netdat1[[2]][,,1],weighted=NULL,mode="undirected")
  vertex.attributes(net)$name
  
  atdata<-data.frame(vertex.attributes(net)$name)
  names(atdata)<-c("id")
  atdata$sex<-rep(NA, nrow(atdata))
  #atdata$pop<-rep(NA, nrow(atdata))
  atdata$fit<-rep(NA, nrow(atdata))
  atdata$preg<-rep(NA, nrow(atdata))
  for(j in 1:nrow(atdata)){
    atdata$sex[j] <- unique(as.character(CMRdata$Sex[CMRdata$id == as.character(atdata$id[j])]))
    #atdata$pop[j] <- unique(as.character(CMRdata$Population_id[CMRdata$id == as.character(atdata$id[j])]))
    atdata$fit[j] <- unique(as.numeric(CMRdata$fit[CMRdata$id == as.character(atdata$id[j])]))
    CMRyear<-CMRdata[which (CMRdata$Year==i),]
    atdata$preg[j] <- unique(as.character(CMRyear$preg[CMRyear$id == as.character(atdata$id[j])]))
    
  }
  
  
  atdata$col<-rep("grey", nrow(atdata))
  atdata[which (atdata$sex=="M"),]$col<-"blue"
  atdata[which (atdata$sex=="F"),]$col<-"purple"

  
  atdata$year<-rep(i,length(vertex_attr(net)$name))
  atdata$tp<-rep((i-(min(Yearlist)-1)),length(vertex_attr(net)$name))
  
  atdata$deg<-igraph::degree(net)
  atdata$wdeg<-igraph::strength(net)
  atdata$eig<-igraph::eigen_centrality(net)$vector
  atdata$bet<-igraph::betweenness(net)
  atdata$pr<-igraph::page.rank(net)$vector
  
  
  
  netdats[[i-(Yearlist[1]-1)]]<-net 
  attlist[[i-(Yearlist[1]-1)]]<-atdata
  
}


#make full study network for full layout
#make one social network for all years
nyear<-12*length(Yearlist)

mindate<-paste0(minyear,"-01-01")
maxdate<-paste0(maxyear+1,"-01-01")
intwindow<-0
netwindow<-36
overlap<-0
spacewindow<-0


netdatcomb<-DynamicNetCreate(data=CMRdata4,
                             intwindow=intwindow,
                             mindate=mindate,
                             maxdate=maxdate, 
                             netwindow=netwindow, 
                             overlap=overlap, 
                             spacewindow=spacewindow)



netall <-graph.adjacency(netdatcomb[[2]][,,1],weighted=NULL,mode="undirected")
vertex.attributes(netall)$name

dataall<-data.frame(vertex.attributes(netall)$name)
names(dataall)<-c("id")

l<-layout.fruchterman.reingold(netall)
l2<-data.frame(l)

dataall$x<-l2$X1
dataall$y<-l2$X2

#get layout for each year
attlist2<-list()
for (i in 1:length(attlist)){
  data<-data.frame(attlist[i])
  if(length(data)>0){
    data$x<-rep(NA, nrow(data))
    data$y<-rep(NA, nrow(data))
    for(j in 1:nrow(data)){
      data$x[j] <- unique(as.character(dataall$x[dataall$id == as.character(data$id[j])]))
      data$y[j] <- unique(as.character(dataall$y[dataall$id == as.character(data$id[j])]))
    }
  }
  attlist2[[i]]<-data
}






#make and save an edge list 

netall <-graph.adjacency(netdatcomb[[2]][,,1],weighted=T,mode="undirected")
e <- get.edgelist(netall)
df <- as.data.frame(cbind(e,E(netall)$weight))

names(df)<-c("id1","id2","weight")

#get count of each ind number of detections 
freq<-as.data.frame(table(CMRdata4$id))
CMRdata4<-merge(CMRdata4, freq, by.x="id", by.y="Var1")
freqdf<-CMRdata4[,c(1,6)]

df$weight<-as.numeric(df$weight)
df$freq1<-rep(NA, nrow(df))
df$freq2<-rep(NA, nrow(df))
for(j in 1:nrow(df)){
  df$freq1[j] <- unique(as.numeric(freqdf$Freq[freqdf$id == as.character(df$id1[j])]))
  df$freq2[j] <- unique(as.numeric(freqdf$Freq[freqdf$id == as.character(df$id2[j])]))
}

df$SRI<-df$weight/(df$freq1+df$freq2-df$weight)

write.csv(df, "SK_allyears_edgelist_SRI.csv", row.names=F)














