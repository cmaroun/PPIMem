#!/usr/bin/perl

#***********************************************************************************
#This script adds the putative interactions from all error rates 
#to the database and checks which one are valid : motif in trans region and motif
#in pdb region 
#Last updated: April 30, 2017
#***********************************************************************************

use strict;

use warnings;

use DBI;

#my $dbh=DBI ->connect('dbi:mysql:db_protein;host=localhost;mysql_socket=/Applications/XAMPP/xamppfiles/var/mysql/mysql.sock' , 'root' , '');
our ($dbh);
require 'connection.pl' || die "error $!" || die " error $@ ";

sub getInteractions {

my $cost=$_[0];

my $pattern1="./" .$cost . "_hits_pattern1.txt";

my $pattern2="./". $cost . "_hits_pattern2.txt";

my $output ="./". $cost . "_hits_stats.txt";

print $output ."\n";

open (PAT1 , $pattern1);

open (PAT2 , $pattern2);

my $header1=<PAT1>;

my $header2=<PAT2>;

while ( my $line=<PAT1> ) {

	my @line_elements=split(/\t/, $line);
	chomp $line_elements[9];

	if ( $line_elements[0] ne $line_elements[4] ) {
    
    my $motif1=$line_elements[1];
    
    my $chain1=$line_elements[2];
    
    my $sth=$dbh->prepare("SELECT id FROM motifs WHERE consensus=? AND chain=?");	
	
	$sth->bind_param(1, $motif1);
	
	$sth->bind_param(2, $chain1);
	
	$sth->execute;
		
	my $motif1_id=$sth->fetchrow;
	
	my $seq=$dbh->prepare("SELECT sequence FROM proteins_ids WHERE protein_id=?");
	
	$seq->bind_param(1, $line_elements[4]);
	
	$seq->execute;
	
	my $seq1=$seq->fetchrow;
	
	$seq1 =~ s/\s//g;
	
	my $seq1_start=index ( $seq1 , $line_elements[6]);
	
	my $reg=$dbh->prepare("SELECT start , end FROM regions WHERE p_id=? AND type='transmembrane'");
	
	$reg->bind_param(1, $line_elements[4]);
	
	$reg->execute;

	my $seq1_end=$seq1_start+length($line_elements[6])-1;
	
	my $residues_trans=0;
	
	my $half=$line_elements[5]/2;
	
	my @trans_region=();
	
	SEARCH: while ( my @row= $reg->fetchrow_array) {
	
			  my $test= join(", ", @row);
			  
			  push @trans_region , $test;		    
	}
	
	for ( my $x =$seq1_start; $x<$seq1_end; $x++ ) {
	
		INNERLOOP: for ( my $y=0 ; $y<scalar(@trans_region) ; $y++) {
		
			my @ends=split(/,/ , $trans_region[$y]);
			
			if ( ($ends[0] <= $x) && ($ends[1] >=$x)) {
				
				$residues_trans++;
			
				last INNERLOOP;
						
			}
		
		}	
	
	}
	
	my $test_available=$dbh->prepare( "SELECT count(*) FROM binding_motifs where motif_id=? AND protein_id=? AND matchedseq=? AND pattern=? ");
	
	$test_available->bind_param(1, $motif1_id);
	
	$test_available->bind_param(2, $line_elements[4]);
	
	$test_available->bind_param(3, $line_elements[6]);
	
	$test_available->bind_param(4, $line_elements[7]);
	
	$test_available->execute;
	
	my $available=$test_available->fetchrow;
	
	if ( $available == 0 ) {
	
	if ( $residues_trans >= 6 ) {
		
		
		
		$dbh->do("INSERT INTO binding_motifs( motif_id , protein_id , cost , matchedseq , error_rate , pattern, trans) VALUES ( '$motif1_id' , '$line_elements[4]' , '$line_elements[8]' , '$line_elements[6]', '$line_elements[9]' , '$line_elements[7]' , '1')");
	
	}
	
	else {
		
		$dbh->do("INSERT INTO binding_motifs ( motif_id , protein_id , cost , matchedseq , error_rate , pattern, trans) VALUES ( '$motif1_id' , '$line_elements[4]' , '$line_elements[8]' , '$line_elements[6]', '$line_elements[9]' , '$line_elements[7]' , '0')");
	
	
	}
	
	}
	
	my $pdb_valid=$dbh->prepare("SELECT pdb_id , start , end FROM pdbs WHERE valid=1 AND p_id=?");
	
	$pdb_valid->bind_param(1,$line_elements[4]);
	
	$pdb_valid->execute;
	
	my @pdb_regions=();
	
	my $residues_pdb=0;
	
	SEARCH: while ( my @row= $pdb_valid->fetchrow_array) {
	
			  my $test= join(", ", @row);
			  
			  push @pdb_regions , $test;		    
	}
	
	for ( my $y=0 ; $y<scalar(@pdb_regions) ; $y++) {
	
		my @ends=split(/,/ , $pdb_regions[$y]);
		
		$residues_pdb=0;
	
		INNERLOOP: for ( my $x =$seq1_start; $x<=$seq1_end; $x++ ) {
		
				my @ends=split(/,/ , $pdb_regions[$y]);
			
				if ( ($ends[1] <= $x) && ($ends[2] >= $x)) {
				
					$residues_pdb++;
						
				}
			
			}
			
			if ( $residues_pdb >= 6) {
			
				$dbh->do("UPDATE binding_motifs SET pdb_valid=1 WHERE motif_id='$motif1_id' AND protein_id='$line_elements[4]' AND error_rate='$line_elements[9]'  AND pattern='$line_elements[7]'"); 
		
    #    		$dbh->do("INSERT INTO pdb_valid  (pdb_id, motif_id , protein_id) VALUES ('$ends[0]' , '$motif1_id' , '$line_elements[4]')");
	
			}
		
		}
	
	}
	   
}

while ( my $line=<PAT2> ) {

	my @line_elements=split(/\t/, $line);
	chomp $line_elements[9];

	if ( $line_elements[0] ne $line_elements[4] ) {
    
    my $motif2=$line_elements[1];
    
    my $chain2=$line_elements[2];
    
    my $sth=$dbh->prepare("SELECT id FROM motifs WHERE consensus=? AND chain=?");	
	
	$sth->bind_param(1, $motif2);
	
	$sth->bind_param(2, $chain2);
	
	$sth->execute;
		
	my $motif2_id=$sth->fetchrow;
	
	my $seq=$dbh->prepare("SELECT sequence FROM proteins_ids WHERE protein_id=?");
	
	$seq->bind_param(1, $line_elements[4]);
	
	$seq->execute;
	
	my $seq2=$seq->fetchrow;
	
	$seq2 =~ s/\s//g;
	
	my $seq2_start=index ( $seq2 , $line_elements[6]);
	
	my $reg=$dbh->prepare("SELECT start , end FROM regions WHERE p_id=? AND type='transmembrane'");
	
	$reg->bind_param(1, $line_elements[4]);
	
	$reg->execute;

	my $seq2_end=$seq2_start+length($line_elements[6])-1;
	
	my $residues_trans=0;
	
	my $half=$line_elements[5]/2;
	
	my @trans_region=();
	
	SEARCH: while ( my @row= $reg->fetchrow_array) {
	
			  my $test= join(", ", @row);
			  
			  push @trans_region , $test;		    
	}
	
	for ( my $x =$seq2_start; $x<=$seq2_end; $x++ ) {
	
		INNERLOOP: for ( my $y=0 ; $y<scalar(@trans_region) ; $y++) {
		
			my @ends=split(/,/ , $trans_region[$y]);
			
			if ( ($ends[0] <= $x) && ($ends[1] >= $x)) {
				
				$residues_trans++;
			
				last INNERLOOP;
						
			}
		
		}	
	
	}
	
	my $test_available=$dbh->prepare( "SELECT count(*) FROM binding_motifs where motif_id=? AND protein_id=? AND matchedseq=? AND pattern=? ");
	
	$test_available->bind_param(1, $motif2_id);
	
	$test_available->bind_param(2, $line_elements[4]);
	
	$test_available->bind_param(3, $line_elements[6]);
	
	$test_available->bind_param(4, $line_elements[7]);
	
	$test_available->execute;
	
	my $available=$test_available->fetchrow;

	if ( $available == 0 ) {
	
	if ( $residues_trans >= 6 ) {
	
		$dbh->do("INSERT INTO binding_motifs ( motif_id , protein_id , cost , matchedseq , error_rate , pattern, trans) VALUES ( '$motif2_id' , '$line_elements[4]' , '$line_elements[8]' , '$line_elements[6]', '$line_elements[9]' , '$line_elements[7]' , '1')");
	
	}
	
	else {
	
		$dbh->do("INSERT INTO binding_motifs ( motif_id , protein_id , cost , matchedseq , error_rate , pattern, trans) VALUES ( '$motif2_id' , '$line_elements[4]' , '$line_elements[8]' , '$line_elements[6]', '$line_elements[9]' , '$line_elements[7]' , '0')");
	
	}
	
	}
	my $pdb_valid=$dbh->prepare("SELECT pdb_id , start , end FROM pdbs WHERE valid=1 AND p_id=?");
	
	$pdb_valid->bind_param(1,$line_elements[4]);
	
	$pdb_valid->execute;
	
	my @pdb_regions=();
	
	my $residues_pdb=0;
	
	SEARCH: while ( my @row= $pdb_valid->fetchrow_array) {
	
			  my $test= join(", ", @row);
			  
			  push @pdb_regions , $test;		    
	}
	
	for ( my $y=0 ; $y<scalar(@pdb_regions) ; $y++) {
	
		my @ends=split(/,/ , $pdb_regions[$y]);
	
		INNERLOOP: 	for ( my $x =$seq2_start; $x<=$seq2_end; $x++ ) {
					
			if ( ($ends[1] <= $x) && ($ends[2] >= $x)) {
				
				$residues_pdb++;
			
						
			}
			
		}
				
		if ( $residues_pdb >= 6) {
			
				$dbh->do("UPDATE binding_motifs SET pdb_valid=1 WHERE motif_id='$motif2_id' AND protein_id='$line_elements[4]' AND error_rate='$line_elements[9]' AND pattern='$line_elements[7]'"); 
			
     #   		$dbh->do("INSERT INTO pdb_valid (pdb_id, motif_id , protein_id) VALUES ('$ends[0]' , '$motif2_id' , '$line_elements[4]')");
			}

		}
	
	}
	
}

print "done\n";

}
#for ( my $percent=0 ; $percent<=0 ; $percent=$percent+0.05) {
for ( my $percent=0 ; $percent<=0.25 ; $percent=$percent+0.05) {

	getInteractions($percent);

}
   