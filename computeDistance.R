#Get all the fixed residues and identify all those that are within 30% similar and group them

motifFile=read.table("./fixedResMotifs.txt" , stringsAsFactors=FALSE)

source("./motifDistances.R")

computeDistance(motifFile$V1)

