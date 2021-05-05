#!/usr/bin/perl

use strict;

use warnings;

open(MOTIFS , "originalMotifs.txt");

open ( FIXED, ">fixedResMotifs.txt");

open ( ALL, ">allResMotifs.txt");

while ( <MOTIFS>) {

	chomp $_;

	my @getMotif=split(/\t/, $_);

	$getMotif[1] =~ s/\.|[0-9]//g;

	$getMotif[1] =~ s/{}//g;

	$getMotif[0] =~ s/\./X/g;

	print FIXED $getMotif[1] . "\n";

	print ALL $getMotif[0] . "\t" . $getMotif[1] . "\n";

}

