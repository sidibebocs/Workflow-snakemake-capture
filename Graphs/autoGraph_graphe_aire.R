#!/usr/bin/env Rscript

library(methods)

args = commandArgs(trailingOnly=TRUE)

source(args[7])
graphe_aire(path=args[1], sorting_col=args[2], graphe_name=args[3], nameX=args[4], nameY=args[5], grapheTitle=args[6])
#source(args[4])
#graphe_aire(path=args[1], sorting_col=args[2], graphe_name=args[3])
