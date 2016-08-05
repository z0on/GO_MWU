go=read.table("MF_heats.csv",header=T)
seqs=c();gos=c()
i="isogroup100"
pb=txtProgressBar(0,length(levels(go$seq)))
for(ii in 1:length(levels(go$seq))){
	setTxtProgressBar(pb,ii)
	i=levels(go$seq)[ii]
	g=subset(go,seq==i)
	seqs=append(seqs,as.character(i))
	gos=append(gos,paste(g$term,collapse=";"))
}
gotab=data.frame(cbind("gene"=seqs,"terms"=gos))
write.table(gotab,file="MF_heats_iso2go.tab",quote=F,row.names=F)
