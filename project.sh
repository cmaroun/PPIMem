#!/usr/bin/bash

#***********************************************************************************
#This script calls all other scripts required to update the TRANSINT database
#Last updated: Jan 3, 2020
#***********************************************************************************

#bash backupDatabase.sh

#echo "Previous version backed up";

perl proteins.pl #Checked

echo "Proteins added to db";
date

perl getPDBSum.pl #Checked

echo "Motifs found\n";
date

perl get_residues.pl # Checked this script get the fixed residues total number or residues of each motif and the sequences and gerenates residues1.txt

echo "Residue of motifs found\n";
date

perl pfams.pl  # Checked this step is needed to get all sequences of the proteins in the database

echo "Proteins and sequences found\n";
date

perl getPdbUrl.pl # Checked This code gets the URL of the pdb website to inlcude as a link in the interaction page

echo "Pdb links\n";
date

perl createClusters.pl #Checked Remove the X in the motif to get only the Fixed residues

echo "Cluster created\n";
date

Rscript computeDistance.R #Checked this compares the fixedresiduesmotifs and groups the ones that have a 30% or more similarity into one line (ClustersAdjusted.txt)

echo "Distances Computed\n";
date

perl createFastaFile.pl #Checked This codes generates a consensus 

echo "Fasta created\n";
date

perl getRegexCons.pl #Checked converts consensus to regex format and add the consensus to the database

echo "Regex Consensus\n";
date

perl GetValidConsensus.pl # Checked replace motifs residues1.txt into the consensus motif !!!!!!	 but fixed residues count are not updated!!!!!!!

echo "Replace consensus\n";
echo "Consenus sequences found\n"
date


Rscript valid_steps.R # finds the valid interactions

echo "Valid interactions found\n";
date 

perl add_interact_db.pl #insert the interaction into the database

echo "Valid interactions added to db\n";
date

Rscript predicted_steps.R #Look for sequences with similar motifs

echo "Putative motifs found\n";
date 

perl predictedProteins.pl

echo "Putative motifs added to db\n";
date

perl aa_score.pl #add the score based on the matrix

echo "Putative motifs score calculated\n";
date 

perl get_files.pl #write the interactions into a file

echo "Interaction files found\n";
date

perl getPdb.pl #check which predicted interactions have valid pdbs

echo "Interactions with pdb found\n";
date

perl get_unique.pl #Supposed to get unique interactions a-b b-a removed. TO CHECK!!!!

echo "Unique interactions found\n";
date