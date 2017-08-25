#!/usr/bin/env Rscript

library(methods)

args = commandArgs(trailingOnly=TRUE)

source(args[3])
graphe_sbe(args[1], args[2])
