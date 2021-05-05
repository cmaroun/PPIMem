#!/usr/bin/perl

#***********************************************************************************
#This script calculates a score based on chemical properties of all putative motifs
#by comparing their fixed residues to those of the original motifs 
#Last updated: April 27, 2017
#***********************************************************************************

use strict;

use warnings;

use DBI;

#my $dbh=DBI ->connect('dbi:mysql:db_protein;host=localhost;mysql_socket=/Applications/XAMPP/xamppfiles/var/mysql/mysql.sock' , 'root' , '');
our ($dbh);
require 'connection.pl' || die "error $!" || die " error $@ ";

open ( MAT , "./aa_matrix.txt");

my $header=<MAT>;

#chomp ($header);

my @elements=split(/\t/ , $header);

#chomp $elements[20];

print  "elements size $#elements @elements\n";	
print "$elements[20]\n";
$elements[20]="V";

my %aa_matrix=();

while (<MAT>) {

	chomp ($_);
 
	my @line=split(/\t/ , $_);
	
	my $a_a=$line[0];
	
	for ( my $x=1; $x<=$#line ; $x++) {
	
	#print "element $x: $elements[$x] & line: $line[0] & value = $line[$x]\n";
		$aa_matrix{$line[0]}{$elements[$x]}=$line[$x];	
	
	}

}


close MAT;

sub getScore{

	my $cost=$_[0];
	
	my $pattern=$_[1];

	my $file ="./" . $cost . "_hits_pattern".$pattern.".txt";

	print $file . "\n";

	open ( FILE , "$file");

	while ( <FILE> ) {

		chomp ($_);
	
		if ( index ($_ , "Uniprot" ) == 0 ) {
	
			next;

		}
		
		else {
		
			my @line_split=split(/\t/ , $_); 
					
			my @motif_split=split(//,$line_split[1]);
			
			my @matched_pattern=split(// , $line_split[6]);
			
			my $temp_position=-1;
			
			my $motif_score=0;
						
			for ( my $x=0; $x<scalar (@motif_split) ; $x++)  {
			
				if ( $motif_split[$x] =~ /[A-Z]/ ) {
					
					$temp_position++;
					
					#print ("in the if: $motif_split[$x] \t $matched_pattern[$temp_position] \t $aa_matrix{$motif_split[$x]}{$matched_pattern[$temp_position]}\n");
					$motif_score = $motif_score + $aa_matrix{$motif_split[$x]}{$matched_pattern[$temp_position]};
					
					}
					
				if ( $motif_split[$x] eq ".") {

					if ( $x == (scalar (@motif_split)-1)) {

						$temp_position++;
					}
				
					elsif  ($motif_split[$x+1] =~ /[A-Z]/) {
					
						$temp_position++;
								
					}
				
					else {
					
						my $number="";
						
						INNER:for ( my $y=$x ; $y<scalar(@motif_split) ; $y++) {
						
							if ( $motif_split[$y] eq "}" ) {				
							
								$temp_position= $temp_position+$number;
								
								$number="";
					
								last INNER;
							
							}
						
							if ( $motif_split[$y] =~ /[0-9]/ ) {
											
								$number= $number . $motif_split[$y];	
							
							}
						
						
						}
					
					
					}
			
				} 
							
			}	
			
		
		my $sth=$dbh->prepare("SELECT id FROM motifs where consensus=? AND chain=?");
		
		$sth->bind_param(1,$line_split[1]);
		
		$sth->bind_param(2,$line_split[2]);
		
		$sth->execute;
		
		my $id=$sth->fetchrow;
		
		$dbh->do("UPDATE binding_motifs SET prop_cost='$motif_score' WHERE protein_id='$line_split[4]' AND motif_id='$id' AND error_rate='$cost'");

	}

	
	}

}

for ( my $percent=0.05; $percent<=0.25 ; $percent=$percent+0.05) {

	getScore($percent, "1");

	getScore($percent, "2");
}
