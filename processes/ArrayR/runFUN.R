x <- readRDS("script/x.RDS")
commonData <- readRDS("script/commonData.RDS")
FUN <- readRDS("script/FUN.RDS")

i <- as.integer(commandArgs(trailingOnly = T)[1])

FUN(x[[i]])
