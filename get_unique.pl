#!/usr/bin/perl

#***********************************************************************************
#This script removes all duplicate interactions from the 50% error rate which 
#includes all the interactions found in TRANSINT
#Last updated: November 6, 2016
#***********************************************************************************


use strict;

use warnings;

#open ( FILE , "./0.25_hits_stats.txt");
open ( FILE , "./New_Interactions_pdbs_0.25.txt");

open ( WRITE , ">./unique_valid_0.25.txt");

open ( WRITE1 , ">./duplicates_valid_0.25.txt");

my %duplicates=();

while (<FILE>) {

	chomp $_;
	
	my $found=0;
	
	my @line_elements=split(/\t/ , $_);
	
	my $motif1=$line_elements[0].":".$line_elements[3]."-".$line_elements[8].":".$line_elements[11];
	print $motif1 ."\n";


			
	my $motif2=$line_elements[8].":".$line_elements[11]."-".$line_elements[0].":".$line_elements[3];
	
	print $motif2 ."\n";		
	if ( defined ($duplicates{$motif1})) {
			
		$found=1;				
			
	}
			
	if (defined ($duplicates{$motif2})) {
		
		$found=1;			
			
	}
			
	if ( $found == 0) {
			
		print WRITE $_ . "\n";
				
		$duplicates{$motif1}=$motif2;
			
	}
			
	else {
			
		print WRITE1 $_."\n";
			
	}


}