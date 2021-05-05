#!/usr/bin/perl

#***********************************************************************************
#This script checks which if the putative interactions have valid pdb structures
#covering the motifs up to the error rate specified in the input file
#Last updated: December 31, 2018
#***********************************************************************************

use strict;

use warnings;

use DBI;

#my $dbh=DBI ->connect('dbi:mysql:db_protein;host=localhost;mysql_socket=/Applications/XAMPP/xamppfiles/var/mysql/mysql.sock' , 'root' , '');
our ($dbh);
require 'connection.pl' || die "error $!" || die " error $@ ";

my @files = <./hits1/*_hits_stats.txt>;

foreach my $x (@files) {

    open (FILE , $x);
    
    my $write="New_Interactions_pdbs_".substr($x, 8,4);
    
    my $write2= "./temp2/Valid_interactions_pdbs_".substr($x,8,4);
    
    print $write . "\n";

    if ( substr( $write , -1 ) eq "_") {

	$write=substr($write, 0, -1);

    }
    
    $write=$write.".txt";
    
    if ( substr( $write2 , -1 ) eq "_") {

	$write2=substr($write2, 0, -1);

    }
    
    $write2=$write2.".txt";
     
    open ( WRITE , ">$write");
    
    open ( WRITE2 , ">$write2");
    
    my $line=<FILE>;

    chomp($line);

    my @header=split(/\t/ , $line);

    print WRITE "protein1\torganism1\tpfam1motif1\tmatched_motif1\tfixed_residues1\tmatched_motif1\tcost1\terror_rate1\tprotein2\torganism2\tpfam2\tmotif2pfam2\tfixed_residues2\tmatched_motif2\tcost2\terror_rate\tpdb1\tpdb2\tvalid\n";

	print WRITE2 "protein1\torganism1\tmotif1\tpfam1\tfixed_residues1\tmatched_motif1\tcost1\terror_rate1\tprotein2\torganism2\tmotif2\tpfam2\tfixed_residues2\tmatched_motif2\tcost2\t\error_rate2\tpdb1\tpdb2\tvalid\n";

    while ( $line = <FILE> ) {

	chomp ($line);

	my @ids= split ( /\t/ , $line);

    my $sth=$dbh->prepare("SELECT distinct p.pdb_id FROM pdb_valid p, motifs m WHERE m.id=p.motif_id AND p.protein_id=? AND motif=?");

	$sth->bind_param(1, $ids[0]);
	
	$sth->bind_param(2, $ids[2]);

	$sth->execute;

	my $pdbs1="";

	while ( my $row= $sth->fetchrow_array) {
	
	   $pdbs1 = $pdbs1 . $row . ";" ;  

	}

	$pdbs1 =substr($pdbs1 , 0 , -1);

	if (!($pdbs1)) {

	    $pdbs1="-";

	}

    $sth=$dbh->prepare("SELECT distinct p.pdb_id FROM pdb_valid p, motifs m WHERE m.id=p.motif_id AND p.protein_id=? AND motif=?");

	$sth->bind_param(1, $ids[8]);
	
	$sth->bind_param(2, $ids[10]);

	$sth->execute;

	my $pdbs2="";

	while ( my $row= $sth->fetchrow_array) {

	    $pdbs2 = $pdbs2 . $row . ";" ;  

	}

	$pdbs2 =substr($pdbs2 , 0 , -1);

	if (!($pdbs2)) {


	    $pdbs2="-";

	}

   $sth=$dbh->prepare("SELECT distinct pf.pfam_id FROM pfams pf WHERE pf.protein_id=?");

	$sth->bind_param(1, $ids[0]);

	$sth->execute;

	my $pfams1="";

	while ( my $row= $sth->fetchrow_array) {

	    $pfams1 = $pfams1 . $row . ";" ;  

	}

	$pfams1 =substr($pfams1 , 0 , -1);
	
	$sth=$dbh->prepare("SELECT distinct pf.pfam_id FROM pfams pf WHERE pf.protein_id=?");

	$sth->bind_param(1, $ids[8]);

	$sth->execute;

	my $pfams2="";

	while ( my $row= $sth->fetchrow_array) {

	    $pfams2 = $pfams2 . $row . ";" ;  

	}

	$pfams2 =substr($pfams2 , 0 , -1);

	$sth=$dbh->prepare("SELECT distinct p.Organism FROM proteins_ids p WHERE p.protein_id=?");

	$sth->bind_param(1, $ids[0]);

	$sth->execute;

	my $organism1="";


	while ( my $row= $sth->fetchrow_array) {

	    $organism1 = $organism1 . $row . ";" ;  

	}
	$organism1 =substr($organism1 , 0 , -1);

	
       	if (!($organism1)) {


	    $organism1="NA";

	}

	$sth=$dbh->prepare("SELECT distinct p.Organism FROM proteins_ids p WHERE p.protein_id=?");

	$sth->bind_param(1, $ids[8]);

	$sth->execute;

	my $organism2="";

	while ( my $row= $sth->fetchrow_array) {

	    $organism2 = $organism2 . $row . ";" ;  

	}
	
	$organism2 =substr($organism2 , 0 , -1);

	if (!($organism2)) {

	    $organism2="NA";

	}
	
	$sth=$dbh->prepare("SELECT distinct p.matchedseq FROM binding_motifs p , motifs m WHERE p.protein_id=? and p.motif_id=id and motif=?");

	$sth->bind_param(1, $ids[0]);
	
	$sth->bind_param(2, $ids[2]);

	$sth->execute;

	my $matched1="";
	
	while ( my $row= $sth->fetchrow_array) {

	    $matched1 = $matched1 . $row . ";" ;  

	}
	
	$matched1 =substr($matched1 , 0 , -1);
		
	$sth=$dbh->prepare("SELECT distinct p.matchedseq FROM binding_motifs p , motifs m WHERE p.protein_id=? and p.motif_id=id and motif=?");

	$sth->bind_param(1, $ids[8]);
	
	$sth->bind_param(2, $ids[10]);

	$sth->execute;

	my $matched2="";
	
	while ( my $row= $sth->fetchrow_array) {

	    $matched2 = $matched2 . $row . ";" ;  

	}
	
	$matched2 =substr($matched2 , 0 , -1);

	if ($pdbs1 eq "-" || $pdbs2 eq "-") {

	    print WRITE $ids[0]."\t" . $organism1. "\t". $pfams1 . "\t". $ids[2]."\t".$ids[3]."\t". $matched1 . "\t". $ids[6] . "\t" .$ids[5]."\t". $ids[8] . "\t" . $organism2 ."\t" . $pfams2. "\t". $ids[10].  "\t". $ids[11]."\t". $matched2. "\t". $ids[14]. "\t". $ids[13] . "\t". $pdbs1 . "\t" . $pdbs2 ."\t". 0 ."\n";
	    

	}

	else {
	
		print WRITE2 $ids[0]."\t" . $organism1. "\t". $pfams1 . "\t". $ids[2]."\t".$ids[3]."\t". $matched1 . "\t". $ids[6] . "\t" .$ids[5]."\t". $ids[8] . "\t" . $organism2 ."\t" . $pfams2. "\t". $ids[10].  "\t". $ids[11]."\t". $matched2. "\t". $ids[14]. "\t". $ids[13] . "\t". $pdbs1 . "\t" . $pdbs2 ."\t". 1 ."\n";
		
          print WRITE $ids[0]."\t" . $organism1. "\t". $pfams1 . "\t". $ids[2]."\t".$ids[3]."\t". $matched1 . "\t". $ids[6] . "\t" .$ids[5]."\t". $ids[8] . "\t" . $organism2 ."\t" . $pfams2. "\t". $ids[10].  "\t". $ids[11]."\t". $matched2. "\t". $ids[14]. "\t". $ids[13] . "\t". $pdbs1 . "\t" . $pdbs2 ."\t". 1 ."\n";
	  
	}

   }

}
