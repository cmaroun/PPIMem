source ("matched.R")

data=read.table("validSummary.txt" , sep="\t", header=TRUE, stringsAsFactors=FALSE)

getMatched(data)