graphe_aire<-function(path, sorting_col, graphe_name, nameX, nameY, grapheTitle)
#graphe_aire<-function(path, sorting_col, graphe_name)
{
	library(ggplot2)
	library(reshape)
  nameGraphePDF<-paste(graphe_name, "pdf", sep = "")
  nameGrapheJPG<-paste(graphe_name, "jpg", sep = "")
	file = read.table(path, sep="\t", header=T)
  file[1] <- NULL
  #file <- file[order(file[,as.numeric(sorting_col)]),]
  #file <- file[order(file[,as.numeric(sorting_col)], -file[,as.numeric(2)], -file[,as.numeric(3)],-file[,as.numeric(4)],-file[,as.numeric(5)],decreasing=T),]
  #file <- file[order(file[,as.numeric(5)], -file[,as.numeric(4)], -file[,as.numeric(3)],-file[,as.numeric(2)],-file[,as.numeric(1)],decreasing=T),]
  file <- file[order(file[,as.numeric(5)], file[,as.numeric(4)], file[,as.numeric(3)],file[,as.numeric(2)],file[,as.numeric(1)],decreasing=T),]
  #file <- file[order(file[,as.numeric(5)], decreasing=F, -file[,as.numeric(4)],decreasing=F),]
  #file <- file[order(file[,as.numeric(5)], file[,as.numeric(4)]),]
	sample.id <- row.names(file)
	row.names(file) <- seq(1:length(sample.id))
	p.dat <- data.frame(step=row.names(file),file,stringsAsFactors=F)
	p.dat <- melt(p.dat,id='step')
	p.dat$step <- as.numeric(p.dat$step)
	p <- ggplot(p.dat, aes(x=step,y=value)) + theme(legend.justification=c(0,1), legend.position=c(0,1), legend.text=element_text(size = 14)) + scale_y_continuous(sec.axis = dup_axis()) + theme(axis.text.x= element_text(face="bold", size=16)) + theme(axis.text.y= element_text(face="bold", size=16), axis.title=element_text(size=18,face="bold")) + guides(fill=guide_legend(title=NULL)) + labs(x=nameX, size=16) + labs(y=nameY, size=16) + labs(title=grapheTitle)
	p + geom_area(aes(fill=variable))
	#ggsave(graphe_name, width = 16, height = 9)
	ggsave(nameGraphePDF, width = 16, height = 9)
	ggsave(nameGrapheJPG, width = 16, height = 9)
return(p + geom_area(aes(fill=variable)))
}
