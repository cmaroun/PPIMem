#!/usr/bin/perl

#***********************************************************************************
#This script calls the database and join all tables needed to construct the
#interactions up until the error rate specified
#Last updated: February 22, 2017
#***********************************************************************************

use strict;

use warnings;

use DBI;

#my $dbh=DBI ->connect('dbi:mysql:transint;host=localhost;mysql_socket=/Applications/XAMPP/xamppfiles/var/mysql/mysql.sock' , 'root' , '');
our ($dbh);
require 'connection.pl' || die "error $!" || die " error $@ ";

sub getInteractions {

my $cost=$_[0];

my $output ="./hits1/". $cost . "_hits_stats.txt";

print $output ."\n";

open (WRITE , ">$output");

print WRITE "protein1\torganism1\tmotif1\tfixed_residues1\tscore\tError_rate1\tcost1\tpdb_valid1\tprotein2\torganism2\tmotif2\tfixed_residues2\tscore\terror_rate2\tcost2\tpdb_valid2\n";

my $query="SELECT
distinct  
`transint`.`p1`.`protein_id` AS `protein_id_A`, `transint`.`p1`.`Organism` AS `organism_A`,
`transint`.`m1`.`motif` AS `motif_A`, `transint`.`m1`.`FIXED_RESIDUES` AS `fixed_residues_A`,
`transint`.`bm1`.`prop_cost` AS `prop_cost_A`, `bm1`.`error_rate` AS `error_rate_A`,
`transint`.`bm1`.`cost` AS `cost_A`, `transint`.`bm1`.`pdb_valid` AS `Valid_A`,
`transint`.`p2`.`protein_id` AS `protein_id_B`, `transint`.`p2`.`Organism` AS `organism_B`,`transint`.`m2`.`motif` AS `motif_B`,
`transint`.`m2`.`FIXED_RESIDUES` AS `fixed_residues_B`,`transint`.`bm2`.`prop_cost` AS `prop_cost_B`,
`transint`.`bm2`.`error_rate` AS `error_rate_B`,`transint`.`bm2`.`cost` AS `cost_B` ,
`transint`.`bm2`.`pdb_valid` AS `Valid_B`
FROM `transint`.`valid_interact` `v`
INNER JOIN `transint`.`motifs` `m1` on `v`.`motif1_id` =`m1`.`id`
INNER JOIN `transint`.`motifs` `m2` on `v`.`motif2_id` =`m2`.`id`
INNER JOIN `transint`.`binding_motifs` `bm1` on `m1`.`id`= `bm1`.`motif_id`
INNER JOIN `transint`.`binding_motifs` `bm2` on `m2`.`id`= `bm2`.`motif_id`
INNER JOIN `transint`.`proteins_ids` `p1` on `bm1`.`protein_id` = `p1`.`protein_id`
INNER JOIN `transint`.`proteins_ids` `p2` on `bm2`.`protein_id` = `p2`.`protein_id`
where
`transint`.`p1`.`Organism` = `transint`.`p2`.`Organism`
and `transint`.`bm2`.`pattern` = 2
and `transint`.`bm1`.`pattern` = 1
and `transint`.`bm1`.`error_rate` <= '$cost'
and `transint`.`bm2`.`error_rate` <= '$cost'
and `transint`.`bm1`.`trans`=1
and `transint`.`bm2`.`trans`=1";

my $sth= $dbh->prepare("$query");

$sth->execute;

while ( my @row=$sth->fetchrow_array) {

    my $line=join ( "\t" , @row);

  	 print WRITE $line ."\n"

}

}

for ( my $percent=0; $percent<=0.25 ; $percent=$percent+0.05) {

	getInteractions($percent);

}
