graphe_sbe<-function(path, outputName)
{
  library(colorspace)
  library(RColorBrewer)
  library(plotrix)
  library(methods)

  fichier = read.table(path, sep="\t", header=T)

  pdf(outputName)

  matplot(fichier, col=brewer.pal(6,"Dark2") ,type="l",lty=1, lwd=1.5, xlab="Base relative position", ylab="Depth", xaxt="n", main="Depth of sequencing at different targeted regions")

  legend('topright', names(fichier), lty=1,  lwd=2, col=brewer.pal(6,"Dark2"), bty='n', cex=.75)

  mtext("Beginning", side=1, adj=0.15, padj=1)
  mtext("Middle", side=1, adj=0.5, padj=1)
  mtext("End", side=1, adj=0.85, padj=1)

  axis.break(axis=1,breakpos=56,bgcol="white",breakcol="black",style="slash",brw=0.02)
  axis.break(axis=1,breakpos=115,bgcol="white",breakcol="black",style="slash",brw=0.02)
  #dev.off ();
}
