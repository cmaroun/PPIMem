#***********************************************************************************
#This script calls the function computeMotifs with different error rates at each run
#in order to get putative motifs from all proteins in the database
#Last updated: February 19, 2017
#***********************************************************************************

source("motif.R")

data=read.table("validSummary.txt" , sep="\t", header=TRUE, stringsAsFactors=FALSE)

sequences=read.table("pfams_ids.txt" , sep="\t" , header=TRUE,stringsAsFactors=FALSE)

computeMotifs(0, data , sequences)
print ("0 done")

computeMotifs(0.05, data , sequences)

print ("0.05 done")

computeMotifs(0.1, data , sequences)

print ("0.1 done")

computeMotifs(0.15, data , sequences)

print ("0.15 done")

computeMotifs(0.2, data , sequences)

print ("0.2 done")

computeMotifs(0.25, data , sequences)

print ("0.25 done")

#computeMotifs(0.3, data , sequences)

#print ("0.3 done")

#computeMotifs(0.35, data , sequences)

#print ("0.35 done")

#computeMotifs(0.4, data , sequences)

#print ("0.4 done")

#computeMotifs(0.45, data , sequences)

#print ("0.45 done")

#computeMotifs(0.5, data , sequences)

#print ("0.5 done")
