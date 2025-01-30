library(ggplot2)
output='./pca_snps.pdf'

v=read.table('faw_bialSNP_pca.eigenvec')
p=read.table('faw_bialSNP_pca.eigenval')
sample=read.table("pop_faw-Copy2.txt",sep="\t",header=F)
colnames(v)[1]='sample_ID'
colnames(sample)[1]='sample_ID'
colnames(sample)[2]='pop'
colnames(sample)[3]='status'


# Identify the repeating pattern (first quarter of the parts)
#sample$sample_ID <- sapply(v$sample_ID, function(x) {
 # parts <- unlist(strsplit(as.character(x), "_"))
 # pattern_length <- length(parts) / 4
 # pattern <- paste(parts[1:pattern_length], collapse = "_")
 # return(pattern)
#})

vs=merge(v,sample,by='sample_ID')

p=ggplot(vs,aes(x=V2,y=V3,col=pop,shape=status))+geom_point()+xlab(paste("PC1, ",format(p$V1[1]*100/sum(p$V1),digits=4),"%",sep=''))+ylab(paste("PC2, ",format(p$V1[2]*100/sum(p$V1),digits=4),"%",sep=''))+ theme(axis.text.x = element_text(size=8),axis.text.y = element_text(size=8),legend.text = element_text(size = 15))+theme_bw()+theme(legend.title = element_blank())

pdf(output,height=5,width=7.2)
p
dev.off()
