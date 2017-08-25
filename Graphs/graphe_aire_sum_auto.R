#This function will take as an input a csv file with values in different columns. The output graphe will consist of the sum of each line, so we first need to modify the data compare to what they look like first. 

graphe_aire<-function(path)
{
	library(ggplot2)
	library(reshape)
  library(methods)
	file = read.table(path, sep="\t", header=T)
  file[1] <- NULL
  for (i in 1:nrow(file)){
    for ( j in 1:(ncol(file)-1) ){
      file[i,j] <- file[i,j]-file[i,j+1]
  }}
  somme <- rowSums(file)
  file <- file[order(rowSums(file), decreasing=F),]
	sample.id <- row.names(file)
	row.names(file) <- seq(1:length(sample.id))
	p.dat <- data.frame(step=row.names(file),file,stringsAsFactors=F)
	p.dat <- melt(p.dat,id='step')
	p.dat$step <- as.numeric(p.dat$step)
	#p <- ggplot(p.dat, aes(x="Sample",y="Number_of_Reads")) + theme(legend.justification=c(0,1), legend.position=c(0,1)) +  scale_fill_brewer(palette="Set2")
	p <- ggplot(p.dat, aes(x=step,y=value)) + theme(legend.justification=c(0,1),  legend.position=c(0,1), legend.text=element_text(size = 14)) +  scale_fill_brewer(palette="Set2") + theme(axis.text.x= element_text(face="bold", size=16)) + theme(axis.text.y= element_text(face="bold", size=16), axis.title=element_text(size=18,face="bold")) + labs(x = "Sample", size=20) + labs(y = "Number of reads", size=20) + labs(title = "Number of reads depending of the different steps of the mapping") +  theme(axis.ticks = element_blank()) + guides(fill=guide_legend(title=NULL))# + scale_y_continuous(position = "right")
	#p <- ggplot(p.dat, aes(x=step,y=value)) + theme(legend.justification=c(0,1), legend.position=c(0,1)) +  scale_fill_brewer(palette="Set2") + theme(axis.text.x= element_text(face="bold", size=16)) + theme(axis.text.y= element_text(face="bold", size=16), axis.title=element_text(size=18,face="bold")) + labs(x = "Sample", size=20) + labs(y = "Number of reads", size=20) + labs(title = "Number of reads depending of the different steps of the mapping") +  theme(axis.ticks = element_blank())
  #p + theme(axis.text.x= element_text(face="bold", size=15))
  #p + theme(axis.text.y= element_text(face="bold", size=15))
	p + geom_area(aes(fill=variable))
	#ggsave("graphe_mapping.pdf", width = 15, height = 11)
	ggsave("graphe_mapping.pdf", width = 16, height = 9)
	ggsave("graphe_mapping.jpg", width = 16, height = 9)
return(p + geom_area(aes(fill=variable)))
}
