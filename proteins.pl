#!/usr/bin/perl

#***********************************************************************************
#This script retrieves information from the Uniprot website concerning transmembrane
#proteins with both IDs GO:0016021 and GO:0005886 or with ID GO:0005887 and stores 
#them in a database.
#Last updated: April 29, 2017
#***********************************************************************************

use strict;

use warnings;

use DBI;

#Each update require a new list of proteins from uniprot 

#system "curl -o uniprot.txt 'http://www.uniprot.org/uniprot/?query=taxonomy%3A%22Eukaryota+%5B2759%5D%22+AND+%28%28go%3A%22integral+component+of+membrane+%5B0016021%5D%22+AND+go%3A%22plasma+membrane+%5B0005886%5D%22+%29OR+go%3A%22integral+component+of+plasma+membrane+%5B0005887%5D%22%29+AND+reviewed%3Ayes&sort=score&format=txt'";

#system "curl -o list_ids.txt 'http://www.uniprot.org/uniprot/?query=taxonomy%3A%22Eukaryota+%5B2759%5D%22+AND+%28%28go%3A%22integral+component+of+membrane+%5B0016021%5D%22+AND+go%3A%22plasma+membrane+%5B0005886%5D%22+%29OR+go%3A%22integral+component+of+plasma+membrane+%5B0005887%5D%22%29+AND+reviewed%3Ayes&sort=score&format=tab&columns=id'";

#open ( WRITE , ">>uniprot.txt");

#print WRITE "ID";

#close WRITE;

#The script starts by emptying the database and recreating the tables in order to
#refill it every time the script runs and update all the information

#my $dbh=DBI ->connect('dbi:mysql:db_protein;host=localhost;mysql_socket=/Applications/XAMPP/xamppfiles/var/mysql/mysql.sock' , 'root' , '');
our ($dbh);
require 'connection.pl' || die "error $!" || die " error $@ ";


#$dbh->do("DROP DATABASE db_protein");
$dbh->do("DROP DATABASE transint");

#$dbh->do("CREATE DATABASE db_protein"); 
$dbh->do("CREATE DATABASE transint"); 

#$dbh->do("USE db_protein"); 
$dbh->do("USE transint"); 

my $create_proteins='CREATE TABLE proteins_ids (protein_id VARCHAR(10) , name VARCHAR(250), Organism 
VARCHAR(250) , sequence TEXT , endo_mito INT , PRIMARY KEY (protein_id), 
INDEX (protein_id) , INDEX (endo_mito))';

$dbh->do("$create_proteins");

my $create_pdbs='CREATE TABLE pdbs ( id INT AUTO_INCREMENT, p_id VARCHAR(10), pdb_id 
VARCHAR(15) ,start INT, end INT, valid INT(1), resolution FLOAT, type VARCHAR (45), 
chain VARCHAR (200) ,PRIMARY KEY(id), FOREIGN KEY(p_id) REFERENCES  proteins_ids 
(protein_id) , INDEX(p_id) , INDEX (pdb_id) , INDEX (valid))';

$dbh->do("$create_pdbs");

my $create_regions='CREATE TABLE regions ( region_id INT AUTO_INCREMENT, p_id VARCHAR(10),
 start INT, end INT, type VARCHAR(15), PRIMARY KEY(region_id), FOREIGN KEY(p_id) 
 REFERENCES proteins_ids (protein_id), INDEX(p_id))';

$dbh->do("$create_regions");

my $create_subunit='CREATE TABLE subunit_type (id INT AUTO_INCREMENT, p_id VARCHAR(10) , 
type VARCHAR(20) , PRIMARY KEY (id) ,  FOREIGN KEY(p_id) REFERENCES proteins_ids 
(protein_id) , INDEX (p_id))';

$dbh->do("$create_subunit");

my $create_interactions='CREATE TABLE interactions ( intact_id INT AUTO_INCREMENT , 
protein1_id VARCHAR(10) , protein2_id VARCHAR(10) , pdb1_valid INT DEFAULT 0 , pdb2_valid 
INT DEFAULT 0, PRIMARY KEY (intact_id) , FOREIGN  KEY (protein1_id) REFERENCES
proteins_ids (protein_id) , INDEX ( protein1_id))';

$dbh->do("$create_interactions");

my $create_pfams_proteins='CREATE TABLE pfams ( table_id INT AUTO_INCREMENT , protein_id 
VARCHAR(10) , pfam_id VARCHAR(15) , PRIMARY KEY ( table_id), FOREIGN KEY (protein_id)
 REFERENCES proteins_ids ( protein_id) , INDEX ( pfam_id) , INDEX ( protein_id))';

$dbh->do("$create_pfams_proteins");

my $createGO='CREATE TABLE go (id INT AUTO_INCREMENT, protein_id VARCHAR(10), go_id VARCHAR(10), 
PRIMARY KEY (id) , FOREIGN KEY (protein_id) REFERENCES proteins_ids (protein_id), INDEX (protein_id), INDEX (go_id))';

$dbh->do("$createGO");

open(FILE , "uniprot.txt");

open (IDS , "list_ids.txt");

open ( PFAMS , ">allPfams.txt");

my @temp=<IDS>; #list of all transmembrane ids

my @pdb_regions=(); #to store all pdb regions of a protein

my @trans_region=(); #to store all transmemebrane regions of a protein

my $id=""; #to store the id of the protein 

my $check="false"; #flag to check if we got an id for the protein or no 

my $organism=""; #to store the organism of the protein

my $sequence=""; #to store the sequence of the protein

my $id_counter=0; #to store the number of proteins

my $pdb_counter=0; #to store the number of pdbs

my $valid_pdb=0; #to store the number of valid pbds 

my %species=(); #to store the number of proteins in each species

my $pdb_id=""; # to store the pdb id we are working with 

my $pdb_check=0; #flag to check if the pdb is valid or no 

my $homo_counter=0; # to store the number of homomeric proteins

my $pdb_homo=0; #to store the number of homomeric proteins with pdbs

my $hetero_counter=0; #to store the number of heteromeric proteins

my $pdb_hetero=0; #to store the number of heteromeric proteins with pdbs

my $homomer=0; #flag to check if protein is homomeric

my $heteromer=0; #flag to check if protein is heteromeric

my %pdbs=(); #to store a pdb in case it is valid

my $added=0; #flag to check whether we already counted the pdb id or no

my $endo=0; #flag to check whether protein is endoplasmi, golgi or mitochondrion

my $pfamsOfProtein=""; #string to store all pfams of protein for later use

print PFAMS "ProteinID\tPfams\n";

while (my $line=<FILE>){

    chomp($line); 
    
    if ( $line =~ m/^ID/) { 

	#lines starting with ID represents the start of a new protein so we add to the database the information of the previous protein

	#we reset the pdb variables to start comparing
	
	$pdb_id=""; 

	$pdb_check=0;

	#save all pfams of protein to file to later filter 

	print PFAMS $id ."\t". $pfamsOfProtein . "\n";

	$pfamsOfProtein="";

	#if the previous protein has a pdb region, we check if it includes the transmembrane region or no 

	if(@pdb_regions){ #we check that the protein has pdb regions

	    #we split each pdb region by the : sign to get the id, region limits and resolution
	    
	    for ( my $x=0 ; $x <@pdb_regions; $x++){

		my @temp=split(/:/,$pdb_regions[$x]); #we split by : because it is the format chosen later in the script, we stored pdb id then start and end and in the end the resolution

		#my $resolution=$temp[2]; 
		my ($resolution, $unit ) = split / /, $temp[2]; 
		
		my $exp=$temp[3];
		
		my $chain = $temp[4];
		
		
		#to make sure we are not counting the same pdb twice for statistical purposes.

		if ($pdb_id ne $temp[0]) { #if the pdbn doesn't already exist, we increment the pdb counter 
		    
		    $pdb_id=$temp[0];

		    $added=0;

		    $pdb_counter++; 

		}

		    $pdb_check=0; #0 means not valid, once we found a valid pdb we change it to one

		my @pdb_ends=split(/-/ , $temp[1]); #we split by - because this is how the start and end were separated later in the script

		#for each pdb region, we split the trans region and compare, if the trans is included in the pdb region it is added in the database with valid=1, otherwise valid=0
		
		for ( my $y=0 ; $y<@trans_region ; $y++ ) {

		    my @split=split(/-/ , $trans_region[$y]); #the format used for each trans region is start-end

		    my $start=$split[0];

		    my $end=$split[1];

		    if ( ($pdb_ends[0]<=$start ) &&( $pdb_ends[1]>=$end)){

			#if the trans region falls in the pdb region then we change the pdb_check to 1
			#we add the pdb region to a hash so that we can keep track of the valid pdbs to later check which interacting proteins have valid pdbs

			$pdb_check=1;	
			
			if ( defined $pdbs{$id} ) {

			    $pdbs{$id}=$pdbs{$id}+1;

			}

			else {
			    
			    $pdbs{$id}=1;

			}
			
		    }

		}

		if ( $exp =~ /Model/ || ($exp =~ /EM/ && $resolution > 3.5 ) ){
		
			$pdb_check=0;
		
		
		}

		#we add the pdb to the database
		
		#$dbh->do("INSERT INTO pdbs ( p_id , pdb_id , start ,end , valid , resolution , type , chain ) VALUES ('$id ' , '$pdb_id' , '$pdb_ends[0]','$pdb_ends[1]','$pdb_check' , '$resolution' , '$exp' , '$chain')");
		#Georges Added this
		if($resolution ne "-"){
			$dbh->do("INSERT INTO pdbs ( p_id , pdb_id , start ,end , valid , resolution , type , chain ) VALUES ('$id' , '$pdb_id' , '$pdb_ends[0]','$pdb_ends[1]','$pdb_check' , '$resolution' , '$exp' , '$chain')");

			#print ("INSERT INTO pdbs ( p_id , pdb_id , start ,end , valid , resolution , type , chain ) VALUES ('". $id ."' , '".$pdb_id."' , '".$pdb_ends[0]."','".$pdb_ends[1]."','".$pdb_check."' , '". $resolution."' , '".$exp."' , '".$chain."')");
		}
		else{

			$dbh->do("INSERT INTO pdbs ( p_id , pdb_id , start ,end , valid ,  type , chain ) VALUES ('$id' , '$pdb_id' , '$pdb_ends[0]','$pdb_ends[1]','$pdb_check' ,  '$exp' , '$chain')");
		}

		if ( $pdb_check==1 && $added==0) {

		    $added=1;
		    $valid_pdb++;
		    

		}

	    
	    }
        
	}
	
	$dbh->do("UPDATE proteins_ids SET endo_mito='$endo' WHERE protein_id='$id'");

	#we reset all variables to null

	@trans_region=();
	
	@pdb_regions=();
		
	$id="";

	$check="false";

	$organism="";   

	$sequence="";

	$homomer=0;

	$heteromer=0;
	
	$endo=0;
	

    }

    if ( index ( $line , "AC ") == 0 ) {

	#each protein can have many accession id and we only want the first one (one protein can have many AC lines, hence we check if we already saved the first one

	if ($check eq "false") {

	    my @ids=split(/;/ , $line); #to split the ids on the same line

	    my @first=split(/   / , $ids[0]); #to take the first id only

	    $id=$first[1];

	    $id=~s/^\s+|\s+$//g;

	    $check="true"; #to stop looping through the rest of the ids

 	    $dbh->do("INSERT INTO proteins_ids ( protein_id ) VALUES ('$id')"); 

	    $id_counter++; #we increment the counter at each protein to eventually get the total number of proteins
	}

  }

  	if ( index ($line, "DE   RecName:")==0) {

  		#One protein can have many name, we only take the first one and not the alterantive ones

  		my @get_name=split(/=/, $line);

  		my $name=substr($get_name[1], 0, -1);

  		$dbh->do("UPDATE proteins_ids SET name=\"$name\" WHERE protein_id='$id'");

  	}
    
    if ( index($line,"DR   PDB;")==0){

	#this line contains all the pdb regions

	#A pdb line contains separated by ; the name , resolution , the way it was discovered(NMR,X-RAY,Model) and the region in case it is not a model

	$pdb_counter++; 

	my @get_pdb=split(/; /,$line); #splits the line to get the pdb id

	my $pdb_id=$get_pdb[1]; #array at 1 contains the id
	
	my $exp=$get_pdb[2];
	
	if (index($get_pdb[4] , ",") != -1) { #when we have many regions they are separated by "," 

	    my @get_regions=split(/,/, $get_pdb[4]); #to get all regions

	    for ( my $x=0 ; $x<@get_regions ; $x++) {

		if ( $x == scalar(@get_regions)-1 ) { #indicating the last region of this pdb

		    my @get_ends=split(/=/, $get_regions[$x]);

		    #each region has the format start-end 

		    my $ends=substr($get_ends[1] , 0, -1);  #we remove the dot at the end of the region
		    		     
		    my $chains=$get_ends[0];
	    
		    $chains =~ s/^\s+|\s+$//g;

		    push @pdb_regions , $pdb_id .":" . $ends.":" . $get_pdb[3] . ":" . $exp . ":" . $chains; #for each region we only want the id, start-end and the resolution
		
		}

		else {

		    my @get_ends=split(/=/, $get_regions[$x]);
		    
		    my $chains=$get_ends[0];
	    
	   		$chains =~ s/^\s+|\s+$//g;

		    push @pdb_regions , $pdb_id .":" . $get_ends[1] .":".$get_pdb[3] . ":" . $exp . ":" . $chains; 

		}

	    }

	}

	else {

	    if ( index ($get_pdb[4] , "-.") != -1) { #indicating no region for this pdb so we skip it

		next;

	    }

	    else { #the pdb has only one region

	    my @get_ends=split(/=/, $get_pdb[4]); #only one region exists

	    my $ends=substr($get_ends[1] , 0, -1); #we remove the dot at the end
	    
	    my $chains=$get_ends[0];
	    
	    $chains =~ s/^\s+|\s+$//g;

	    push @pdb_regions , $pdb_id .":" . $ends .":" . $get_pdb[3] . ":" . $exp . ":" . $chains;  #for each region we only want the id, start-end and the resolution

		 }
	}

    }
    
    if ( $line =~ /DR   GO/ ) {
    
    	if ( $line =~ /[Ee]ndoplasmic reticulum/){
    
    		$endo=1;
    	}
    	
    	if ( $line =~ /[Mm]itochondrion/ ) {
    	
    		$endo=1;
    	
    	}
    	
    	if ( $line =~ /[Gg]olgi/ ) {
    	
    	
    		$endo=1;
    	
    	}

    	if ( $line =~ /C:/ ) { 

    		#we get all GO ids of cellular component that doesn't correspond to our GOs of interest

    		my @line_split=split(/;/ , $line);

    		my @get_id=split(/:/ , $line_split[1]);

    		$line_split[1]=~s/^\s+|\s+$//g; #in case there is any leading or trailing spaces

    		if ( $get_id[1] ne "0016021" && $get_id[1] ne "0005886" && $get_id[1] ne "0005887") {

    			$dbh->do("INSERT INTO go (protein_id , go_id) VALUES ('$id' , '$get_id[1]')");

    		}

    	}
    
    
    }
    
    if ( index ( $line , "DR   Pfam") == 0 ) {
    
    
    	my @line_split=split(/; /, $line );

    	if ( $pfamsOfProtein eq "") {

    		$pfamsOfProtein=$line_split[1];

    	}

    	else {

    	$pfamsOfProtein= $pfamsOfProtein . ";" . $line_split[1];

    	}
		
		$dbh->do("INSERT INTO pfams (protein_id , pfam_id) VALUES ('$id' , '$line_split[1]')");
    
    }
    
    if (index ( $line , "FT   TRANSMEM") == 0) {

	#this line contains the transmembrane region

	my $start=substr($line , 13, 7);

	$start=~s/^\s+|\s+$//g; #we remove the leading and trailing spaces from the start 

	my $end=substr($line ,20 , 7);

	$end=~s/^\s+|\s+$//g; #we remove the leading and trailing spaces from the end

	my $type="Transmembrane"; #it can be either transmembrane or intramembrane

	#sometimes, the start and end regins can have "<" or "<" so we remove them in order to be able to compare them later on

	if ( index ( $start, "<") != -1){ 

	    $start= substr($start,1);
 
	}

	if ( index ($end , ">") != -1 ) {

	    $end=substr($end,1);
     
	}

	#we add the region to the database and to the array of regions to use them to validate proteins at the end 

	$dbh->do("INSERT INTO regions(p_id , start , end , type ) VALUES ('$id' , '$start' , '$end' , '$type')");

	push @trans_region , $start."-".$end;

    }

    if ( index ( $line , "FT   INTRAMEM") == 0){

	#this line contains the intramembrane, we get the start and the end.
	
	my $type="Intramembrane";

	my $start=substr($line , 13, 7); 

	$start=~s/^\s+|\s+$//g; #we remove the leading and trailing spaces from the start

	my $end=substr($line ,20 , 7);

	$end=~s/^\s+|\s+$//g; #we remove the leading and trailing spaces from the end

	#we add the region to the database

	$dbh->do("INSERT INTO regions(p_id , start , end , type ) VALUES ('$id' , '$start' , '$end' , '$type')");




    }

    if ( index ($line , "OS " ) == 0) {

	#this line is to get the organism of the protein in question

	$organism=substr($line , 5);

	$organism=~ s/^\s+|\s+$//g; #we remove the leading and trailing spaces

	my $next_line=<FILE>;

	chomp($next_line);

	#in case the organism name is too long and spans more than one line

	while ( index ( $next_line , "OS ") == 0) {

	    my $temp=substr($next_line,5);
	    
	    $temp=~s/^\s+|\s+$//g;

	    $organism=$organism." ".$temp;

	    $next_line=<FILE>;

	    chomp($next_line);

	}

	$organism=substr($organism , 0, -1);

	#we update the proteins_ids table since we previously didn't add the organism name

	$dbh->do("UPDATE proteins_ids SET organism=\"$organism\" WHERE protein_id='$id'");

	if ( defined $species{$organism}) { #we add the organism to the hash to keep track of the number of occurrences of each one

	    #if the organism exists, we increment its value

	    $species{$organism}=$species{$organism}+1;

	}

	else {

	    #else we set it to one

	    $species{$organism}=1;

	}

    }

    if ( index ( $line , "SQ   SEQUENCE ") == 0 ) {

	# this line contains the sequence of the protein

	my $new_line=<FILE>;

	chomp ($new_line);
	
	while ( index ( $new_line , "     ") == 0) {

	    #in case the sequence spans more than one line

	    chomp ($new_line);

	    $new_line=~ s/\s+//g;

	    $sequence=$sequence.$new_line;

	    $new_line=<FILE>;	    

	}

	$dbh->do("UPDATE proteins_ids SET sequence='$sequence' WHERE protein_id='$id'");
	

    }

    if ( index ( $line , "CC   -!- INTERACTION")==0){

	#this line signals the start of the interacting proteins

	my $new_line=<FILE>;

	while ( index ($new_line , "CC      ")==0) { #we keep on looping to get all the interacting proteins since each line contains only one 

	    #we loop through the next lines to get the id of all the interacting proteins ( each one contains one id)

	    	my @int_id=split(/:/, $new_line); 

		my @ids=split(/       / , $int_id[0]); #the line starts by spaces

		#ids[1] contains the id of the protein in question

	        if ( index ($ids[1] , "Self" ) == 0){ #it interacts with itself

		    my @self_int=split(/;/,$ids[1]);

 		    $dbh->do("INSERT INTO interactions ( protein1_id , protein2_id) VALUES ( '$id' , '$id')");

		}

		else {

		    #we only add the interactions between transmembrane proteins

		    if ( grep {/$ids[1]/} @temp ) {  

 			$dbh->do("INSERT INTO interactions( protein1_id ,protein2_id) VALUES ('$id' , '$ids[1]')");

		    }

		}

		$new_line=<FILE>; 


	}


    }

    if ( index ($line , "CC   -!- SUBUNIT") == 0 ) {

	#this line gives us information about the protein, if it is homomeric or heterometic

	my $homo="homomer";

	my $hetero="heteromer";

 	if ( index ( $line , "homo")!= -1) { 

	    $dbh->do("INSERT INTO subunit_type (p_id, type) VALUES ( '$id' , '$homo')");
 
	    $homomer=1;

	    $homo_counter++;

	}

	if ( index ( $line , "hetero") != -1 ) {

	    $dbh->do("INSERT INTO subunit_type (p_id, type) VALUES ( '$id' , '$hetero')");

	    $heteromer=1;

	    $hetero_counter++;

	}

    }
	

}
 
for ( keys %pdbs ) {

    #this is used to check which interacting proteins have valid pdbs

    $dbh->do("UPDATE interactions SET pdb1_valid='1' WHERE protein1_id='$_'"); #this updates all interacting proteins 1

    $dbh->do("UPDATE interactions SET pdb2_valid='1' WHERE protein2_id='$_'"); #this updates all interacting proteins 2

}

my $size = keys %pdbs;


close FILE;

close IDS;

