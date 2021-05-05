#!/usr/bin/perl

 #***********************************************************************************
#This script is needed to get the URL for an image for each transmembrane PDB to be
#displayed on the interface of TRANSINT
#Last updated: January 6, 2017
#***********************************************************************************

use strict;

use warnings;

use DBI;

use LWP::Simple;

#my $dbh=DBI ->connect('dbi:mysql:db_protein;host=localhost;mysql_socket=/Applications/XAMPP/xamppfiles/var/mysql/mysql.sock' , 'root' , '');
our ($dbh);
require 'connection.pl' || die "error $!" || die " error $@ ";

$dbh->do("ALTER TABLE pdbs ADD pdb_url VARCHAR(250)");

my $query= "SELECT DISTINCT pdb_id FROM pdbs WHERE valid=1";

my $sth=$dbh->prepare($query);
 
$sth->execute;
 
my @pdbs=();
 
 while ( my $row= $sth->fetchrow_array) {
 
 	push @pdbs , $row;
 
 }
 
for ( my $x=0 ; $x<scalar(@pdbs) ; $x++ ) {

	my $id=lc $pdbs[$x];
	
	my $id_sub=substr($id , 1, 2);

	my $url = "http://cdn.rcsb.org/images/hd/".$id_sub."/".$id."/".$id.".0__chimera_tm_surface_ligand_75_75.png";
	
	my $url1="http://www.rcsb.org/pdb/images/".$pdbs[$x]."_asym_r_500.jpg";

	if (head($url)) {

		$dbh->do("UPDATE pdbs SET pdb_url=\"$url\" where pdb_id='$pdbs[$x]'");
	
	}

	elsif ( head($url1)) {

		$dbh->do("UPDATE pdbs SET pdb_url=\"$url1\" where pdb_id='$pdbs[$x]'");

	}

	else {


	}

}