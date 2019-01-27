#setwd("~/Dropbox/Documents/dixon_2017_RT-GBM/reciprocal_transplant_methylationV2/clean_july12/scripts")

geres=read.table("MWU_CC_heat_lpv.csv",header=T)
gbmres=read.table("MWU_CC_surv_lpv.csv",header=T)

goods=intersect(geres$term,gbmres$term)
#goods=unique(as.character(c(geres$term[geres$p.adj<=0.1],gbmres$term[gbmres$p.adj<=0.1])))
length(goods)

geres=geres[geres$term %in% goods,]
gbmres=gbmres[gbmres$term %in% goods,]

# all overlapping GO terms
ress=merge(geres,gbmres,by="term")
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

