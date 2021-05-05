#!/usr/bin/perl

use strict;

use warnings;

use DBI;

#my $dbh=DBI ->connect('dbi:mysql:db_protein;host=localhost;mysql_socket=/Applications/XAMPP/xamppfiles/var/mysql/mysql.sock' , 'root' , '');
our ($dbh);
require 'connection.pl' || die "error $!" || die " error $@ ";

open ( WRITE , ">pfams.txt");

my $sth=$dbh->prepare("SELECT distinct pfam_id FROM pfams");
	
$sth->execute;
	
while ( my $row= $sth->fetchrow_array) {

	print WRITE $row . "\n";
	
}

close WRITE;

open ( FILE , "pfams.txt");

open ( WRITE , ">pfams_ids.txt");

print WRITE "Pfam\tUniprotId\tSequence\n";

while ( <FILE> ) {

	chomp ($_);
	
	my $sth=$dbh->prepare("SELECT pf.pfam_id , pf.protein_id , p.sequence FROM pfams pf , proteins_ids p  WHERE  p.protein_id=pf.protein_id AND pf.pfam_id=?");
	
	$sth->bind_param(1, $_);
		
	$sth->execute;
	
	while ( my @row= $sth->fetchrow_array) {
	
				  my $test= join("\t", @row);
			  
					print WRITE $test . "\n";
	
	
	}
}


