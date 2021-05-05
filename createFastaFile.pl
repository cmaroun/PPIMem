	#!/usr/bin/perl

	use warnings;

	use strict;

	open(FILE, "./clustersAdjusted.txt");

	open (MOTIF , "./allResMotifs.txt");

	open(WRITE , ">finalCons.txt");

	my %motifs=();

	while (<MOTIF>) {

		chomp($_);

		my @line=split(/\t/, $_);

		push(@{$motifs{$line[1]}}, $line[0]);

	}

	#print "Printing Motifs\n";
	#print map { "$_ => @{$motifs{$_}}\n" } keys %motifs;

	close MOTIF;

	while (<FILE>){

		chomp $_;

		open(TEMP, ">./tempfasta.fasta");

		my @line=split(/;/, $_);

		if ( scalar(@line) == 1) {

			foreach my $motif (@{$motifs{$line[0]}}) {

			print WRITE $_ . "\t" . $motif ."\n";

			}

		}

		else {

			for ( my $x=0; $x<scalar(@line) ; $x++ ){
				#print "In the else: @line\n";
				my $tempKey= "";
				foreach my $key (@{$motifs{$line[$x]}}) {
					
					if($key ne $tempKey){
					
					print TEMP ">Protein" . $x . "\n" . $key . "\n";
					}
					$tempKey=$key;
				}
			}

			system ("Rscript getConsensus.R");

			open ( CONS , "tempCons.txt");

			my $line1=<CONS>;

			my $line2=<CONS>;

			chomp $line2;

			print WRITE $_ . "\t" . substr($line2,2) ."\n";	

		}
		
	}