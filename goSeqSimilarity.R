go=read.table("MF_heats.csv",header=T)
nseq=length(levels(go$seq))
sim=matrix(nrow=nseq,ncol=nseq,0)
colnames(sim)=row.names(sim)=levels(go$seq)
pb=txtProgressBar(0,length(levels(go$term)))
for(ii in 1:length(levels(go$term))){
	setTxtProgressBar(pb,ii)
	i=levels(go$term)[ii]
	g=subset(go,term==i)
	sim[as.character(g$seq),as.character(g$seq)]=sim[as.character(g$seq),as.character(g$seq)]+1
}

signedPs=read.csv("heats.csv")
signedPs$pvals=10^(-abs(signedPs$logP))
similarity=sim
signedPs=signedPs[signedPs$gene %in% colnames(similarity),]
signedPs=signedPs[signedPs$gene==colnames(similarity),]
save(similarity,signedPs,file="Pvals_Similarities.RData")
nrow(pvalues)