#!/usr/bin/env Rscript

library(methods)

args = commandArgs(trailingOnly=TRUE)

#if (args[2] == "graphe_aire_sum_auto.R"){
  #source(args[2])
  #graphe_aire(args[1])
#} else{
  #source(args[2])
  #graphe_aire(args[1])
#}
source(args[2])
graphe_aire(args[1])
