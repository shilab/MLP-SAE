setwd($dir);
library(ggplot2);
library(reshape);
library(grid);
library(futile.logger);
library(VennDiagram);

true <- read.table("output_true.txt",header = F);
predic <- read.table("output dropout.txt",header = F);
predic1 <- read.table("output.txt",header=F);
true <- as.matrix(true);
predic <- as.matrix(predic);
predic1 <- as.matrix(predic1);

posi_withdropout = 0;
nega_withdropout = 0;
for(i in 1:nrow(true)){
  for(j in 1:ncol(true)){
    if (true[i,j] * predic[i,j] >= 0) {posi_withdropout = posi_withdropout +1;}
    if (true[i,j] * predic[i,j] < 0)  {nega_withdropout = nega_withdropout +1;}
  }
}
###posi_withdropout: 54305; nega_withdropout: 25027

posi_withoutdropout = 0;
nega_withoutdropout = 0;
for(i in 1:nrow(true)){
  for(j in 1:ncol(true)){
    if (true[i,j] * predic1[i,j] >= 0) {posi_withoutdropout = posi_withoutdropout +1;}
    if (true[i,j] * predic1[i,j] < 0)  {nega_withoutdropout = nega_withoutdropout +1;}
  }
}
###posi_withoutdropout:54284;nega_withoutdropout:25048


corr <- matrix(NA, nrow(true),1);
for(i in 1:nrow(true)){
  corr[i] <- cor.test(true[i,],predic[i,],alternative = "two.sided",method = "pearson")$estimate;
}
outlier <- which(is.na(corr));

corr1 <- matrix(NA, nrow(true),1);
for(i in 1:nrow(true)){
  corr1[i] <- cor.test(true[i,],predic1[i,],alternative = "two.sided",method = "pearson")$estimate;
}
outlier1 <- which(is.na(corr1));
index_corr<-NULL;
index_corr1<-NULL;
index_corr <- cbind(matrix(seq(1,6611),nrow(as.matrix(corr)),1),as.matrix(corr));
index_corr1 <- cbind(matrix(seq(1,6611),nrow(as.matrix(corr1)),1),as.matrix(corr1));
sub <- index_corr[which(index_corr[,2] <=0.05 & index_corr[,2] > 0),1];

corr_squared <- corr*corr; #### with drop out;
corr_squared1 <- corr1*corr1;#### without drop out;
index_corr<-NULL;
index_corr1<-NULL;
index_corr <- cbind(matrix(seq(1,6611),nrow(as.matrix(corr_squared)),1),as.matrix(corr_squared));
index_corr1 <- cbind(matrix(seq(1,6611),nrow(as.matrix(corr_squared1)),1),as.matrix(corr_squared1));
sub = sort(corr_squared,decreasing = TRUE,index.return=TRUE);
sub1 = sort(corr_squared1,decreasing = TRUE,index.return=TRUE);
subindex_corr = sub$ix[1:200];
subindex_corr1 = sub1$ix[1:200];
common <- intersect(sub,sub1);

df.cut_Dropout <- as.data.frame(table(cut(corr_squared[,1], breaks=seq(0,1, by=0.10))));
df.cut <- as.data.frame(table(cut(corr_squared1[,1], breaks=seq(0,1, by=0.10))));

me1 <- read.table("corr1_log",header=T,row.names=1);
qnames = c("without_Dropout","with_Dropout");
bin = rownames(me1);
me1 <- as.matrix(me1);
#me = structure(c(me$pre,me$true),.Dim=c(6611,2))
me1.df = melt(me1);
pdf("correlation8_log.pdf");
ps <- ggplot(me1.df,aes(x=X1, y=value)) + 
  geom_line(aes(group=X2,color=X2),size=1) + 
  geom_point(aes(shape=X2,color=X2),size=1.5) +
  scale_color_discrete(name="",labels=c("MLP-SAE with dropout","MLP-SAE")) +
  scale_shape(name="",labels=c("MLP-SAE with dropout","MLP-SAE")) +
  scale_y_continuous(breaks = seq(0,4, by = 1),limits = c(0, 4),expand = c(0.01,0.01)) + 
  scale_x_discrete(breaks = unique(me1.df$X1)[seq(1,8,1)]) +
  xlab(bquote("R"^"2")) + ylab('log(Number of genes)\n') + 
  theme(legend.text = element_text(size = 8),legend.background = element_rect(fill = 'white'), legend.position = "bottom", legend.title = element_blank()) + 
  theme(panel.grid.major = element_line(colour = "grey", linetype = "dotted",size = 0.3), panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank()) +
  theme(panel.background = element_rect(fill = "white",colour = "grey50")) +
  theme(panel.border = element_rect(linetype = "solid", fill = NA)) +
  theme(plot.background = element_rect(fill = "white", colour = "grey50")) + 
  theme(plot.caption=element_text(size=8, hjust=1.0, margin=margin(t=10)));
ps
dev.off()

pdf("R-squared.pdf");
histinfo <- hist(corr_squared, breaks=seq(0,1,by=0.1),main='R-squared distribution',xlab="R-squared",ylim=c(0,4800),xaxt="n",yaxt="n");
axis(side=1, at=seq(0,1,by=0.1));
axis(side=2, at=seq(0,4800,by=400));
plot(histinfo)
dev.off();

gene_full <- read.table("~/Dropbox (UNC Charlotte)/Xie2017BMCgenomics/inputDataBren2005/expression_matrix",header = T,row.names = 1);
annotation <- read.table("~/Dropbox (UNC Charlotte)/Xie2017BMCgenomics/inputDataBren2005/gene_anno2",header = F,fill=TRUE,sep="\t",col.names=c(1:3));
gene_full <- as.matrix(gene_full);
set <- cbind(gene_full,annotation);

sub_gene <- NULL;
index <- NULL;
for(i in 1:nrow(set)){
  #if (sum(is.na(t(gene_full[i,]))) != ncol(gene_full) ){sub_gene <- cbind(sub_gene,gene_full[i,])}
  if (sum(is.na(t(set[i,1:112]))) != ncol(gene_full) && sum(is.na(t(set[i,1:112]))) != 111 && sum(is.na(t(set[i,1:112]))) != 110 && sum(is.na(t(set[i,1:112]))) != 109 && sum(is.na(t(set[i,1:112]))) != 108 && sum(is.na(t(set[i,1:112]))) != 107 && sum(is.na(t(set[i,1:112]))) != 106)
      {sub_gene <- rbind(sub_gene,set[i,]);}
}

sub_gene2<-NULL;
for(i in 1:nrow(sub_gene)){
  if((is.na(sub_gene[i,ncol(sub_gene)]) && nchar(as.character(sub_gene[i,ncol(sub_gene)-2])) == 2) || (is.na(sub_gene[i,ncol(sub_gene)]) && as.character(sub_gene[i,ncol(sub_gene)-2]) == "empty")){sub_gene2 <- rbind(sub_gene2,sub_gene[i,]);}
}
write.table(sub_gene,"sub_gene",quote = F,sep="\t",col.names = F,row.names = F);

aver_true= matrix(NA,nrow = nrow(true),1);
aver_pre = matrix(NA,nrow = nrow(predic),1);
for(k in 1:nrow(true)){
  aver_true[k] = mean(true[k,]);
  aver_pre[k] = mean(predic[k,]);
}

me <- read.table("fig7",header=T,row.names = 1);
qnames = c("pre","tru");
index = rownames(me);
me <- as.matrix(me);
me.df = melt(me)


pdf("fig7-8.pdf",width = 14, height = 5);
p <- ggplot() + 
  geom_bar(data = me.df, aes(x=X1, y=value,fill=X2),size=1.55, stat = "identity",position='identity') + 
  scale_fill_discrete(name="",breaks=c("pre","tru"),labels=c("Estimated Gene Expression","True Gene Expression")) +
  scale_y_continuous(breaks = seq(-5,6, by = 1),expand = c(0.005,0.005)) + scale_x_continuous(breaks = seq(0,6620, by = 200),expand = c(0.005,0.005)) + 
  scale_color_manual(values=c("#0000FF", "#FF00FF")) +
  xlab('\nGenes') + ylab('Average of Gene Expression Value across Samples\n') + 
  ggtitle("Comparison between True Gene Expression and Estimated Expression\n") + 
  theme(plot.title = element_text(face="bold",size = 10,hjust = 0.5)) + 
  theme(axis.text.x=element_text(size=7, color = "black")) + 
  theme(axis.text.y=element_text(size=9, color = "black")) + 
  theme(axis.title.x = element_text(vjust=-0.70,size=9), axis.title.y = element_text(vjust=0.55,size=9)) + 
  theme(legend.position = "bottom",legend.margin = margin(2, 2, 2, 2)) + 
  theme(legend.text = element_text(size = 8, colour = "black")) +
  theme(plot.caption=element_text(size=6, hjust=0, margin=margin(t=10))) + 
  theme(panel.background = element_rect(fill = "white",colour = "grey50")) +
  theme(panel.border = element_rect(linetype = "solid", fill = NA)) +
  theme(plot.background = element_rect(fill = "white", colour = "grey50")) + 
  theme(plot.caption=element_text(size=8, hjust=1.0, margin=margin(t=10))) 
p
dev.off()

me1 <- read.table("subfig7",header=T,row.names = 1);
qnames = c("pre","tru");
index = rownames(me1);
me1 <- as.matrix(me1);
#me = structure(c(me$pre,me$true),.Dim=c(6611,2))
me1.df = melt(me1);
pdf("fig7-8_sub.pdf",width = 14, height = 5);
ps <- ggplot(me1.df,aes(x=X1, y=value)) + 
  geom_line(aes(group=X2,color=X2)) + 
  geom_point(aes(shape=X2,color=X2)) +
  geom_hline(yintercept = 0,size=0.4) +
  scale_color_discrete(name="",labels=c("Estimated Gene Expression","True Gene Expression")) +
  scale_shape(name="",labels=c("Estimated Gene Expression","True Gene Expression")) +
  scale_y_continuous(breaks = seq(-1.0,1.5, by = 0.5),limits = c(-1.0, 1.5),expand = c(0.005,0.005)) + 
  scale_x_continuous(breaks = seq(800,880, by = 10),expand = c(0.005,0.005)) + 
  xlab('\nSelected Genes') + ylab('Average of Gene Expression Value across Samples\n') + 
  theme(legend.text = element_text(size = 8),legend.background = element_rect(fill = 'white'), legend.position = "bottom", legend.title = element_blank()) + 
  theme(panel.background = element_rect(fill = "white",colour = "grey50")) +
  theme(panel.border = element_rect(linetype = "solid", fill = NA)) +
  theme(plot.background = element_rect(fill = "white", colour = "grey50")) + 
  theme(plot.caption=element_text(size=8, hjust=1.0, margin=margin(t=10)));
ps
dev.off()

