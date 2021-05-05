	use strict;

	use warnings;

	use DBI;

	#my $dbh=DBI ->connect('dbi:mysql:db_protein;host=localhost;mysql_socket=/Applications/XAMPP/xamppfiles/var/mysql/mysql.sock' , 'root' , '');
    our ($dbh);
    require 'connection.pl' || die "error $!" || die " error $@ ";

	open ( FILE , "./finalCons.txt");

	open (WRITE , ">finalCountMotifs.txt");

	open (MOTIF , "./allResMotifs.txt");

	my %motifs=();

	#Read the motif file and fill in the motif hash
	while (<MOTIF>) {

		chomp($_);

		my @line=split(/\t/, $_);

		# if (defined($motifs{$line[1]})) {

			push(@{$motifs{$line[1]}}, $line[0]);
		#}

		#else {

		#$motifs{$line[1]}=$line[0];

		#}

	}

	close MOTIF;

	while ( <FILE>) {

		chomp $_;

		my @line=split(/\t/ , $_);

		my @numOfMotifs=split(/;/, $line[0]);

		my @motifChars=split(// , $line[1]);

		my $motif="";

		if ( scalar(@numOfMotifs) > 0) {

			my $counter=0;

			for (my $char=0; $char<scalar(@motifChars); $char++) { 

				#if ($motifChars[$char] =~ /[A-Z]/) { CHANGED IT TO THE BELOW TO ACCOUNT ONLY FOR AAs
				if ($motifChars[$char] =~ /[GALMFWKQESPVICYHRNDT]/) {

					$motif .= $motifChars[$char];

					$counter++;
				}

				else {

					$motif .= ".";

				}

			}

			$motif=~s/((\.)\2+)/$2 . "{".length($1)."}"/ge;
			$motif=~s/^\.\{\d+\}/''/ge;
			$motif=~s/\.\{\d+\}$/''/ge;
			$motif=~s/^\./''/ge;

			print WRITE $line[0] . "\t" . scalar(@numOfMotifs) . "\t" . $motif . "\n";
			#print "BEFORE THE IF @numOfMotifs\n";
			for ( my $originalMotif=0 ; $originalMotif<scalar(@numOfMotifs) ; $originalMotif++) {
				
				if (scalar(@{$motifs{$numOfMotifs[$originalMotif]}})==1) {

					my $motifInQuestion=@{$motifs{$numOfMotifs[$originalMotif]}}[0];
					
					$motifInQuestion=~s/X/"."/ge;
					$motifInQuestion=~s/((\.)\2+)/$2 . "{".length($1)."}"/ge;
					$motifInQuestion=~s/^\.\{\d+\}//ge;
					$motifInQuestion=~s/\.\{\d+\}$//ge;
					$motifInQuestion=~s/^\.//ge;
					

					
					#print "in the if !!!!!!!!!!!!!!!!!!!!!!!!!!!!$motifInQuestion \n";
					#print $counter . "\n";

					$dbh->do("UPDATE motifs SET consensus='$motif' where motif='$motifInQuestion'");

					$dbh->do("UPDATE motifs SET fixed_residues='$counter' where motif='$motifInQuestion'");

				}

				else {

					#print "in the else\n";
					foreach (@{$motifs{$numOfMotifs[$originalMotif]}}) {

						#print "$motif\n";

						my $motifInQuestion = $_;
						$motifInQuestion=~s/X/"."/ge;
						$motifInQuestion=~s/((\.)\2+)/$2 . "{".length($1)."}"/ge;
						$motifInQuestion=~s/^\.\{\d+\}//ge;
						$motifInQuestion=~s/\.\{\d+\}$//ge;
						$motifInQuestion=~s/^\.//ge;
						
						#print "Original : $motifInQuestion\n";
						

						$dbh->do("UPDATE motifs SET consensus='$motif' where motif='$motifInQuestion'");

						$dbh->do("UPDATE motifs SET fixed_residues='$counter' where motif='$motifInQuestion'");

					}
				}

			}
		}

		else {

			$motif=$line[1];

			print WRITE $line[0] . "\t" . scalar(@numOfMotifs) . "\t" . $motif . "\n";

			for ( my $originalMotif=0 ; $originalMotif<scalar(@numOfMotifs) ; $originalMotif++) {

				my $motifInQuestion=@{motifs{$numOfMotifs[$originalMotif]}};

				$dbh->do("UPDATE motifs SET consensus='$motif' where motif='$motifInQuestion'");

			}

		}


	}