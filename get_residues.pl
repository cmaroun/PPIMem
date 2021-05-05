#!/usr/bin/perl

use strict;

use warnings;

use DBI;

#*****************************************************************************
#This script gets the number of fixed residues and total number of residues of 
#each motif of the non bonded contacts as well as their sequence.
#*****************************************************************************

open (FILE , "interaction_list1.txt");

open(WRITE , ">residues1.txt");

open ( WRITE1 , ">residues1_with_missing.txt");

#my $dbh=DBI ->connect('dbi:mysql:db_protein;host=localhost;mysql_socket=/Applications/XAMPP/xamppfiles/var/mysql/mysql.sock' , 'root' , '');
our ($dbh);
require 'connection.pl' || die "error $!" || die " error $@ ";

print WRITE "PDB_ID\tPROTEIN1\tPROTEIN2\tCHAIN1\tCHAIN2\tPATTERN1\tPATTERN2\tFIXED_RESIDUES_P1\tFIXED_RESIDUES_P2\tTOTAL_RESIDUES_P1\tTOTAL_RESIDUES_P2\tSequence1\tSequence2";

print WRITE1 "PDB_ID\tPROTEIN1\tPROTEIN2\tCHAIN1\tCHAIN2\tPATTERN1\tPATTERN2\tFIXED_RESIDUES_P1\tFIXED_RESIDUES_P2\tTOTAL_RESIDUES_P1\tTOTAL_RESIDUES_P2\tSequence1\tSequence2";

sub get_fixed_length {

#this sub gets the fixed length of every motif
    
    my $length=0;

    my ($id, $pattern)=@_;

    my @residues=split(//, $pattern); 

    #we split each pattern to get every residue alone and loop through each one. if the residue is an amino acid we add it to the length else we skip it

    $length=0;
    
   for ( my $x =0 ; $x<@residues ; $x++ ) {

	if ( $residues[$x] =~ /[ARNDCQEGHILKMFPSTWYV]/ ) { 

	    if ( $x<@residues-1 && $residues[$x+1] eq "{" ) { #this if checks if the residue is found two or more times consecutively 

		$length=$length+$residues[$x+2];
		
	    }

	    else {

		$length=$length+1;	        


	    }
	    
	}
	   
   }

    #we add the length to the line

	if ( $id ne "-" ) { 

  	 print WRITE  $length."\t";
  	 
  	}
   
	$dbh->do("UPDATE motifs SET fixed_residues='$length' WHERE motif='$pattern'");
  
   print WRITE1  $length."\t";

}

sub get_total_residues {

#this sub gets the total length of every motif

    my $length=0;
    
    my ($id,$pattern)=@_;

    my @residues=split(//, $pattern);

    $length=0;
    
   for ( my $x =0 ; $x<@residues ; $x++ ) {

       if ( $x<@residues-1 && $residues[$x+1] eq "{" ) { #if the residue or . is found two or more times consecutively

	   my $sum="";

	   my $y=$x+2;

	   while ( $residues[$y] ne "}" ) {

	       $sum=$sum.$residues[$y];

	       $y++;

	   }

	   $x=$y;
	   
	 $length=$length+$sum;  
		
       }

	    else {

		$length=$length+1;	        


	    }
	    
   }

    # we add the length to the line
    
   if ( $id ne "-" ) { 
   
   	print WRITE  $length."\t";
   
   }
   
   print WRITE1  $length."\t";
	   
}

sub get_sequence {

    #we use curl to go to the UniProt database and only pick the proteins that are reviewed and are transmembrane
    
    my ($protein_id1, $protein_id2) = @_;
    
    my $sequence1= "NA";
    
    my  $sequence2= "NA" ;
    
    if($protein_id1 ne "-" ){

   		my $sth=$dbh->prepare("SELECT sequence FROM proteins_ids WHERE protein_id=?");
   	
  	 	$sth->bind_param(1, $protein_id1);
   	
  		$sth->execute;
   	
  	 	$sequence1=$sth->fetchrow;
  	 	
  	 	$sequence1 =~s/\s+//g;
   	
   	}
   	
   	if($protein_id2 ne "-" ){
   	
   		my $sth=$dbh->prepare("SELECT sequence FROM proteins_ids WHERE protein_id=?");
   	
  	 	$sth->bind_param(1, $protein_id2);
   	
  	 	$sth->execute;
   	
   		$sequence2=$sth->fetchrow;
   		
   		$sequence2 =~ s/\s+//g;
   		
   	} 	
   	
   	print WRITE1 $sequence1 . "\t" . $sequence2;
   	
   	if ( $sequence1 ne "NA" && $sequence2 ne "NA" ) {
   	  	
   		print WRITE $sequence1 . "\t" . $sequence2;	
   	
   	}
}

while ( my $line=<FILE>) {

    chomp($line);  

    my @temp=split(/\t/ , $line);

    #we get each interaction, split it and only consider the non bonded interactions

    
		if ( $temp[1] ne "-" && $temp[0] ne "PDB_ID" && $temp[2] ne "-" ) {

			$temp[1] =~ s/\s+//g;
			
			$temp[2] =~ s/\s+//g;

			print WRITE "\n" . $temp[0]."\t".$temp[1]."\t".$temp[2]."\t".$temp[3]."\t".$temp[4]."\t".$temp[5]."\t". $temp[6]."\t";
			
		}
	
  if($temp[0] ne "PDB_ID"){
	   
     print WRITE1 "\n" . $temp[0]."\t".$temp[1]."\t".$temp[2]."\t".$temp[3]."\t".$temp[4]."\t".$temp[5]."\t". $temp[6]."\t";
	
  }
	
  if ( $temp[1] ne "-" && $temp[0] ne "PDB_ID" && $temp[2] ne "-" ) {

    print $temp[1] . "\n";

	 get_fixed_length($temp[1],$temp[5]);  #fixed length of pattern 1

	 get_fixed_length($temp[2] , $temp[6]); #total length of pattern 1

	 get_total_residues($temp[1], $temp[5]); #fixed length of pattern 2

	 get_total_residues($temp[2],$temp[6]); #total length of pattern 2
	 
	 get_sequence($temp[1] , $temp[2]); #sequence of the protein
	 
	 }
	 
	 #if ( $temp[1] ne "-" && $temp[2] ne "-" ) {

		#print WRITE "\n";
		
	#}
	
	#print WRITE1 "\n";
	   
	
}


close WRITE;

close WRITE1;