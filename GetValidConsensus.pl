#!/usr/bin/perl

use strict;

use warnings;

use DBI;

#my $dbh=DBI ->connect('dbi:mysql:db_protein;host=localhost;mysql_socket=/Applications/XAMPP/xamppfiles/var/mysql/mysql.sock' , 'root' , '');
our ($dbh);
require 'connection.pl' || die "error $!" || die " error $@ ";

open(FILE , "residues1.txt");

open(WRITE, ">validSummary.txt");

print WRITE "PDB_ID\tPROTEIN1\tPROTEIN2\tCHAIN1\tCHAIN2\tCONSENSUS1\tCONSENSUS2\tFIXED_RESIDUES_P1\tFIXED_RESIDUES_P2\tTOTAL_RESIDUES_P1\tTOTAL_RESIDUES_P2\tSequence1\tSequence2\n";

my $header=<FILE>;

while (<FILE>) {

	chomp $_;

	my @line=split(/\t/,$_);

	my $consensus1="";

	my $sth=$dbh->prepare("SELECT consensus FROM motifs WHERE motif=?");
   	
  	$sth->bind_param(1, $line[5]);
   	
  	$sth->execute;
   	
  	$consensus1=$sth->fetchrow;
  	
  	my $consensus2="";

	$sth=$dbh->prepare("SELECT consensus FROM motifs WHERE motif=?");
   	
  	$sth->bind_param(1, $line[6]);
   	
  	$sth->execute;
   	
  	$consensus2=$sth->fetchrow;

  	print WRITE $line[0] . "\t" . $line[1] . "\t" . $line[2] . "\t" . $line[3] . "\t" . $line[4] . "\t" . $consensus1 . "\t" . $consensus2 . "\t" . $line[7] . "\t" .$line[8] . "\t" . $line[9] . "\t" . $line[10] . "\t" . $line[11] . "\t" . $line[12]. "\n";
}