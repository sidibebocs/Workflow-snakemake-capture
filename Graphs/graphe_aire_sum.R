graphe_aire<-function(path)
{
	library(ggplot2)
	library(reshape)
	file = read.table(path, sep="\t", header=T)
  file[1] <- NULL
  somme <- rowSums(file)
  file <- file[order(rowSums(file), decreasing=F),]
	sample.id <- row.names(file)
	row.names(file) <- seq(1:length(sample.id))
	p.dat <- data.frame(step=row.names(file),file,stringsAsFactors=F)
	p.dat <- melt(p.dat,id='step')
	p.dat$step <- as.numeric(p.dat$step)
	p <- ggplot(p.dat, aes(x=step,y=value)) + theme(legend.justification=c(0,1), legend.position=c(0,1)) +  scale_fill_brewer(palette="Set2")
	p + geom_area(aes(fill=variable))
	ggsave("graphe.pdf", width = 15, height = 11)
	#ggsave("graphe.pdf", width = 15, height = 11)
return(p + geom_area(aes(fill=variable)))
}
