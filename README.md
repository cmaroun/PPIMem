# PPIMem
# Welcome to PPIMem

PPIMem is an algorithm for predicting the interactions between transmembrane (TM) proteins of the plasma membrane of eukaryotes.
PPIMem is a knowledge-based approach based on the detection of α-helical TM non-bonded contact residues between different chains in experimental 3D structures of membrane protein multimers. The pipeline of programs collects, filters, and processes the input data and, based on the assumption that homologs of the benchmark template interface contact motifs are expected to interact analogously, proposes TM proteins forming a putative complex. The database with the annotated predicted interactions is implemented as a web application that supports sorting and filtering, as well as exporting of the results as a csv file.


# Understanding PPIMem results

Each line of the PPIMem database webpage corresponds to a predicted binary complex between Protein A and Protein B. The menu on the left allows filtering the columns on the right. The first column corresponds to the Organism name. The next eight columns are for Protein A, and represent its UniProtKB AC, its Name, its Gene, its contact Consensus Motif in regular expression format, the Number of Contact Residues in the Consensus Motif, the applied percent Mutation Rate, the Cost parameter (number of allowed amino acid mutations for a given number of contact residues and mutation rate), and the Valid parameter (1 if the experimental structure of the TM protein is available in the PDB; 0 otherwise). The values of several of these variables may be modified through a sliding window. The following eight columns correspond to the same variables for Protein B. The value of the last column is 0 if Valid A = Valid B = 0; 1 if Valid A or Valid B are 1; and 2 if Valid A = Valid B = 1, that is, when each of the two TM proteins possesses an experimental structure. In this latter case, the user can proceed to a docking simulation of Protein A and Protein B. Two options are available in the top of the webpage: “Show X entries,” that is, the number of entries that can be shown in the page (up to 100), and a Search engine for searching any character (motif, digit, rate, protein name, etc.) in the entire database. The bottom of the webpage lists the number of entries resulting from a given choice of parameters (left), and the corresponding number of pages (left). This last option allows to pass from one page of the results to another.


# Why is PPIMem original ?

Most PPI databases are mostly concerned with soluble globular proteins. The originality of PPIMem is that of proposing binary molecular complexes between TM proteins based on membrane-exposed packing motifs at the interface. PPIMem reveals 98 α-helix interacting motifs and 1504 unique membrane proteins across 39 species involved in the predicted complexes. The number of PPIMem-predicted membrane heteromer interactions is of 21,544 (9,797 human).


# Running the pipeline that predicts the membrane protein - membrane protein complexes.  

The project.sh file calls all the scripts required to update the database. The scripts have been uploaded in this project. The system requirements are a Linux computer.
The different databases, such as UniProt, PDB, PDBsum, OPM, and Pfam should be downloaded locally.


# The PPIMem results are implemented as a web application

The resulting database with the annotated predicted interactions is implemented as a web application that supports sorting and filtering. The output data can be downloaded as a csv file and the predictions can be accessed for the time being at https://transint.shinyapps.io/transint/.
