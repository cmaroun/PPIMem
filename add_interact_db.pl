#!/usr/bin/perl

#***********************************************************************************
#This script adds all motifs and valid interactions to the database and calculates
#a score for each interaction
#Last updated: December 18, 2016
#***********************************************************************************

use strict;

use warnings;

use DBI;

#my $dbh=DBI ->connect('dbi:mysql:db_protein;host=localhost;mysql_socket=/Applications/XAMPP/xamppfiles/var/mysql/mysql.sock' , 'root' , '');
our ($dbh);
require 'connection.pl' || die "error $!" || die " error $@ ";

my $drop_valid='DROP TABLE pdb_valid';

$dbh->do("$drop_valid");

my $drop_table='DROP TABLE binding_motifs';

$dbh->do("$drop_table");

my $create_table='CREATE TABLE binding_motifs ( interact_id INT AUTO_INCREMENT, motif_id INT , protein_id VARCHAR(10) , cost INT , matchedseq TEXT, error_rate DOUBLE , pattern INT , prop_cost INT , trans BOOLEAN , pdb_valid INT default 0, PRIMARY KEY ( interact_id) , FOREIGN KEY (motif_id ) REFERENCES motifs (id) , FOREIGN KEY (protein_id) REFERENCES proteins_ids (protein_id) , index ( motif_id) , index( protein_id))';

$dbh->do("$create_table");

my $drop_interact='DROP TABLE valid_interact';

$dbh->do("$drop_interact");

my $create_interact='CREATE TABLE valid_interact ( binding_id INT AUTO_INCREMENT , motif1_id INT , motif2_id INT , PRIMARY KEY ( binding_id))';

$dbh->do("$create_interact");

my $create_valid='CREATE TABLE pdb_valid (table_id INT AUTO_INCREMENT, PRIMARY KEY (table_id) , pdb_id VARCHAR (5) , motif_id INT , protein_id VARCHAR(10) , foreign key (motif_id) references binding_motifs ( motif_id) , foreign key (protein_id) references binding_motifs ( protein_id) , foreign key (pdb_id) references pdbs ( pdb_id))';

$dbh->do("$create_valid");

open ( INTER , "validated_interactions.txt");

my $header = <INTER>;

while (<INTER>) {

	chomp ($_);
	
	my @elements=split(/\t/ , $_);
	
	my $motif1=$elements[3];
	
	my $motif2= $elements[12];
	
	my $sth=$dbh->prepare("SELECT id FROM motifs WHERE consensus=? AND chain=?");
	
	$sth->bind_param(1, $motif1);
	
	$sth->bind_param(2, $elements[2]);
	
	$sth->execute;
		
	my $motif1_id=$sth->fetchrow;
	
	$sth=$dbh->prepare("SELECT FIXED_RESIDUES FROM motifs where id=?");
	
	$sth->bind_param(1, $motif1_id);
	
	$sth->execute;
	
	my $count=$sth->fetchrow;

	print $motif1_id . "\n";
	
	my $prop_cost=$count*2;
		
	if ( $elements[5] eq "NA" ) {
	
		$dbh->do("INSERT INTO binding_motifs ( motif_id , protein_id , cost , matchedseq , error_rate , pattern, trans , prop_cost) VALUES ( '$motif1_id' , '$elements[0]' , '$elements[6]' , '$elements[5]', '0' , '1' , '0' , '$prop_cost')");
	
	}
	
	else {
	
		my $seq=$dbh->prepare("SELECT sequence FROM proteins_ids WHERE protein_id=?");
	
		$seq->bind_param(1, $elements[0]);
	
		$seq->execute;
	
		my $seq1=$seq->fetchrow;
	
		$seq1 =~ s/\s//g;
	
		my $seq1_start=index ( $seq1 , $elements[5]);
	
		my $reg=$dbh->prepare("SELECT start , end FROM regions WHERE p_id=? AND type='transmembrane'");
	
		$reg->bind_param(1, $elements[0]);
	
		$reg->execute;

		my $seq1_end=$seq1_start+length($elements[5])-1;
	
		my $half=$elements[4]/2;
	
		my @trans_region=();
		
		my $residues_trans=0;
	
		SEARCH: while ( my @row= $reg->fetchrow_array) {
	
				  	my $test= join(", ", @row);
			  
			  		push @trans_region , $test;		    

		}
	
		for ( my $x = $seq1_start; $x<=$seq1_end; $x++ ) {
	
			INNERLOOP: for ( my $y=0 ; $y<scalar(@trans_region) ; $y++) {
		
				my @ends=split(/,/ , $trans_region[$y]);
			
				if ( ($ends[0] <= $x) && ($ends[1] >= $x)) {
				
					$residues_trans++;
			
					last INNERLOOP;
						
				}
		
			}	
	
		}
	
		
		if ( $residues_trans >= 15) {
	
			$dbh->do("INSERT INTO binding_motifs ( motif_id , protein_id , cost , matchedseq , error_rate , pattern, trans , prop_cost) VALUES ( '$motif1_id' , '$elements[0]' , '$elements[6]' , '$elements[5]', '0' , '1' , '1' , '$prop_cost')");
	
		}
	
		else {

			$dbh->do("INSERT INTO binding_motifs ( motif_id , protein_id , cost , matchedseq , error_rate , pattern, trans, prop_cost) VALUES ( '$motif1_id' , '$elements[0]' , '$elements[6]' , '$elements[5]', '0' , '1' , '0' , '$prop_cost')");
	
		}
		
	my $pdb_valid=$dbh->prepare("SELECT pdb_id , start , end FROM pdbs WHERE valid=1 AND p_id=?");
	
	$pdb_valid->bind_param(1,$elements[0]);
	
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
			
			if ( ($ends[1] <= $x) && ($ends[2] >= $x)) {
			
				$residues_pdb++;
						
			}
			
		}
			
			#if ( $residues_pdb >= 15) {
			if ( $residues_pdb >= 10) {
	#here uncommented			
				$dbh->do("UPDATE binding_motifs SET pdb_valid=1 WHERE motif_id='$motif1_id' AND protein_id='$elements[0]' AND error_rate='0' AND pattern='1'"); 
	#here uncommented
  				$dbh->do("INSERT INTO pdb_valid  (pdb_id, motif_id , protein_id) VALUES ('$ends[0]' , '$motif1_id' , '$elements[0]')");

			}
			
		}
		
	}
	
	$sth=$dbh->prepare("SELECT id FROM motifs WHERE consensus=? AND chain=?");
	
	$sth->bind_param(1, $motif2);
	
	$sth->bind_param(2, $elements[11]);
	
	$sth->execute;
	
	my $motif2_id=$sth->fetchrow;
	
	$sth=$dbh->prepare("SELECT FIXED_RESIDUES FROM motifs where id=?");
	
	$sth->bind_param(1, $motif2_id);
	
	$sth->execute;
	
	$count=$sth->fetchrow;
	
	$prop_cost=$count*2;
	
	if ( $elements[12] eq "NA" ) {
	
		$dbh->do("INSERT INTO binding_motifs ( motif_id , protein_id , cost , matchedseq , error_rate , pattern, trans, prop_cost) VALUES ( '$motif1_id' , '$elements[0]' , '$elements[6]' , '$elements[5]', '0' , '1' , '0', '$prop_cost')");
	
	}
	
	else {
	
		my $seq=$dbh->prepare("SELECT sequence FROM proteins_ids WHERE protein_id=?");
	
		$seq->bind_param(1, $elements[9]);
	
		$seq->execute;
	
		my $seq2=$seq->fetchrow;
	
		$seq2 =~ s/\s//g;
	
		my $seq2_start=index ( $seq2 , $elements[12]);
	
		my $reg1=$dbh->prepare("SELECT start , end FROM regions WHERE p_id=? AND type='transmembrane'");
	
		$reg1->bind_param(1, $elements[9]);
	
		$reg1->execute;
	
		my $seq2_end=$seq2_start+length($elements[15])-1;
	
		my $residues_trans=0;
	
		my $half=$elements[13]/2;
	
		my @trans_region=();
	
		SEARCH: while ( my @row= $reg1->fetchrow_array) {
	
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
		
		if (  $residues_trans >= 15 ) {
		
			$dbh->do("INSERT INTO binding_motifs ( motif_id , protein_id , cost , matchedseq , error_rate , pattern, trans, prop_cost ) VALUES ( '$motif2_id' , '$elements[9]' , '$elements[15]' , '$elements[14]', '0' , '2' , '1' , '$prop_cost')");	
	
		}
	
		else {
	
			$dbh->do("INSERT INTO binding_motifs ( motif_id , protein_id , cost , matchedseq , error_rate , pattern , trans , prop_cost) VALUES ( '$motif2_id' , '$elements[9]' , '$elements[15]' , '$elements[14]', '0' , '2' , '0', '$prop_cost')");
		
		}
		
	my $pdb_valid=$dbh->prepare("SELECT pdb_id , start , end FROM pdbs WHERE valid=1 AND p_id=?");
	
	$pdb_valid->bind_param(1,$elements[9]);
	
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
		
		INNERLOOP: for ( my $x =$seq2_start; $x<=$seq2_end; $x++ ) {
			
			if ( ($ends[1] <= $x) && ($ends[2] >= $x)) {
				
				$residues_pdb++;		
						
			}
			
		}
			
			#if ( $residues_pdb >= 15) {
			if ( $residues_pdb >= 10) {
			
				$dbh->do("UPDATE binding_motifs SET pdb_valid=1 WHERE motif_id='$motif2_id' AND protein_id='$elements[9]' AND error_rate='0'  AND pattern='2'"); 
     
     			$dbh->do("INSERT INTO pdb_valid  (pdb_id, motif_id , protein_id) VALUES ('$ends[0]' , '$motif2_id' , '$elements[9]')");
	
			}
			
		}
	
		
	}

		$dbh->do("INSERT INTO valid_interact ( motif1_id , motif2_id ) VALUES ( '$motif1_id' , '$motif2_id' )");	

}	
