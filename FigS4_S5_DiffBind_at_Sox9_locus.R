library(DiffBind)
dbObj <- dba(sampleSheet="enhancer_MC_Fore_20220826_Sox9_Candidate_enhancers.csv")
dbObj <- dba.count(dbObj, bUseSummarizeOverlaps=TRUE,summits=FALSE,filter=0,minCount=1)
dba.plotPCA(dbObj,  attributes=DBA_FACTOR, label=DBA_ID)
plot(dbObj)
dbObj
# Establishing a contrast 
dbObj <- dba.contrast(dbObj, categories=DBA_FACTOR,minMembers = 2)
dbObj <- dba.analyze(dbObj, method=DBA_EDGER)
#  summary of results
dba.show(dbObj, bContrasts=T)
#  overlapping peaks identified by the two different tools (DESeq2 and edgeR)
dba.plotVenn(dbObj,contrast=1,method=DBA_ALL_METHODS)
comp1.edgeR <- dba.report(dbObj, method=DBA_EDGER, contrast = 1, th=1)
comp1.deseq <- dba.report(dbObj, method=DBA_DESEQ2, contrast = 1, th=1)
# EdgeR
out <- as.data.frame(comp1.edgeR)
write.table(out, file="results_20220827_Whole_genome_SEACR/Forelimb_vs_MC_edgeR.txt", sep="\t", quote=F, col.names = NA)
# DESeq2
out <- as.data.frame(comp1.deseq)
write.table(out, file="results_20220827_Whole_genome_SEACR/Forelimb_vs_MC_deseq2.txt", sep="\t", quote=F, col.names = NA)
# Create bed files for each keeping only significant peaks (p < 0.05)
# EdgeR
out <- as.data.frame(comp1.edgeR)
edge.bed <- out[ which(out$FDR < 0.05), 
                 c("seqnames", "start", "end", "strand", "Fold")]
write.table(edge.bed, file="results_20220827_Whole_genome_SEACR/Forelimb_vs_MC_edgeR_sig.bed", sep="\t", quote=F, row.names=F, col.names=F)

# DESeq2
out <- as.data.frame(comp1.deseq)
deseq.bed <- out[ which(out$FDR < 0.05), 
                  c("seqnames", "start", "end", "strand", "Fold")]
write.table(deseq.bed, file="results_20220827_Whole_genome_SEACR/Forelimb_vs_MC_deseq2_sig.bed", sep="\t", quote=F, row.names=F, col.names=F)

