#setwd("~/Dropbox/Documents/dixon_2017_RT-GBM/reciprocal_transplant_methylationV2/clean_july12/scripts")

data1=read.table("MWU_MF_heat_lpv.csv",header=T)
data2=read.table("MWU_MF_surv_lpv.csv",header=T)

goods=intersect(data1$term,data2$term)
#goods=unique(as.character(c(data1$term[data1$p.adj<=0.1],data2$term[data2$p.adj<=0.1])))
length(goods)

data1=data1[data1$term %in% goods,]
data2=data2[data2$term %in% goods,]

# all overlapping GO terms
ress=merge(data1,data2,by="term")
plot(delta.rank.x~delta.rank.y,ress,xlab="GBM", ylab="GE",mgp=c(2.3,1,0))
abline(v=0,lty=3)
abline(h=0,lty=3)

# GO terms highly signifcant in any of the two datasets
sigs=(ress$p.adj.x<=0.01 | ress$p.adj.y<=0.01)
sum(sigs) # 71
plot(delta.rank.x~delta.rank.y,ress[sigs,],xlab="GBM", ylab="GE",mgp=c(2.3,1,0))
abline(v=0,lty=3)
abline(h=0,lty=3)

# GO terms signifcant in both datasets
sigs=(ress$p.adj.x<=0.1 & ress$p.adj.y<=0.1)
sum(sigs) # 20
plot(delta.rank.x~delta.rank.y,ress[sigs,],xlab="GBM", ylab="GE",mgp=c(2.3,1,0))
abline(v=0,lty=3)
abline(h=0,lty=3)

