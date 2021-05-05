	#!/usr/bin/perl

	use strict;

	use warnings;

	use DBI;

	use LWP::Simple;

	#my $dbh=DBI ->connect('dbi:mysql:db_protein;host=localhost;mysql_socket=/Applications/XAMPP/xamppfiles/var/mysql/mysql.sock' , 'root' , '');
	#my $dbh=DBI ->connect('dbi:mysql:db_protein:host=linuxdev.accbyblos.lau.edu.lb', 'gkhazen' , 'GKh@Z1');
    our ($dbh);
    require 'connection.pl' || die "error $!" || die " error $@ ";

	$dbh->do("DROP TABLE IF EXISTS binding_motifs");

	$dbh->do("DROP TABLE IF EXISTS valid_interact");

	$dbh->do("DROP TABLE IF EXISTS motifs");

	$dbh->do("CREATE TABLE motifs (id INT AUTO_INCREMENT , motif TEXT, CHAIN VARCHAR(5), FIXED_RESIDUES INT, consensus TEXT, PRIMARY KEY (id)  , INDEX (id) )");

	open (LIST, ">interaction_list1.txt");

	print LIST "PDB_ID\tUNIPROT_ID1\tUNIPROT_ID2\tCHAIN1\tCHAIN2\tPATTERN1\tPATTERN2\n"; 

	open ( ORIG , ">originalMotifs.txt");

	open(FILEOPM , "pdbsOPM.txt" );

	my %TMRegions=();

	my %mapping=();

	sub getPDBsumFile {

		my ($pdb , @chains)=@_;

		for ( my $chain1=0 ; $chain1<(scalar(@chains)-1) ; $chain1++) {

			for ( my $chain2=$chain1+1; $chain2<scalar(@chains); $chain2++) {

				my $url="http://www.ebi.ac.uk/thornton-srv/databases/cgi-bin/pdbsum/GetIface.pl?pdb=" . $pdb . "&chain1=". $chains[$chain1] ."&chain2=" .$chains[$chain2];

				my $outputfile= "./PDBsumFiles/" .$pdb."-".$chains[$chain1].$chains[$chain2].".txt";

	#			system "curl -o $outputfile '$url'";

			}
		}
	}

	sub checkCorrectPDBsum {

		my @filenames=<./PDBsumFiles/*.txt>;

		foreach my $file ( @filenames) {

			my $size=-s $file;

			if ( $size == 238 ) {

				unlink $file;
			}
		}
	}

	sub getPDBFile {

		my ( $pdb_id , $chain) = @_;

		my $url="http://www.ebi.ac.uk/pdbe/entry-files/download/pdb". $pdb_id .".ent";

		my $outputfile=$pdb_id. ".ent";

		system "curl -o $outputfile '$url'";

		open ( PDB , $outputfile);

		my @PDBFile=<PDB>;

		my @atoms=grep /^ATOM/i, @PDBFile;

		my @chain_atoms=grep /\s$chain\s/i , @atoms;

		my $chain_start=substr($chain_atoms[1], 22, 6);

		$chain_start=~ s/^\s+|\s+$//g;

		my $chain_end=substr($chain_atoms[scalar(@chain_atoms)-1], 22,6);

		$chain_end=~ s/^\s+|\s+$//g;

		my $returnString=$chain_start . "-" .$chain_end;

		return $returnString;

	}

sub getMotifs {

	my %a_a=("ALA" =>"A","ARG" =>"R", "ASN" =>"N" , "ASP" =>"D", "CYS" => "C", "GLU" =>"E", "GLN" => "Q", "GLY" =>"G" ,"HIS"=>"H" ,"ILE"=>"I", "LEU"=>"L","LYS"=>"K","MET"=>"M","PHE"=>"F", "PRO"=>"P","SER"=>"S", "THR"=>"T", "TRP"=>"W","TYR"=>"Y","VAL"=>"V", "UNK"=>".");

	my @filenames=<./PDBsumFiles/*.txt>;

	foreach my $file ( @filenames) {

		my $chain1=substr($file, 19,1);

		my $chain2=substr($file,20,1);

		my $pdb_id=substr($file , 14,4);

		my $toSearch1=lc($pdb_id) . " " . uc($chain1);

		my $toSearch2=lc($pdb_id) . " " . uc($chain2);

		my $mol1=""; my $mol2="";

		my %residues_1=(); my $mol1_start=0; my $check=0;

		my %residues_2=(); my $mol2_start=0; my $bond="";

		my @uniprot_ids1=(); my @uniprot_ids2=();

		my @regions1=(); my @tm_regions1=();

		my @regions2=(); my @tm_regions2=();

		open ( FILETMSegments , $file ) || die ( "Error opening file!");

		OUTER:while ( my $line=<FILETMSegments> ) {

			chomp ($line);

			$line =~ s/$1/\n/g if ($line =~ m/(\r\n?|\n\r?)/);

			if  ( index ($line , " ") == -1 && index ($line , "-") == -1 ) {

				if ( $check ==1 ) {

					my $temp_start=0;

					my $size = keys %residues_1;
					#print  "the size variable 1 is $size\n";
					#print "AA in resiudes 1: ";
					#print %residues_1;
					#print "\n";

					if ( $size <= 5) {
						print "$pdb_id Motif 1 shorter than 5 \n";
						next OUTER;

					}

					my $chain1_start="";

					my $chain1_end="";

					my $get_regions1=$TMRegions{$toSearch1};
					#print "protein to search " . $toSearch1 ."\n";
					#print "Regions 1"  . " $get_regions1 " . "\n";
					@regions1=split(/,/, $get_regions1);

					foreach my $region (@regions1) {

						my @coordinates=split(/-/, $region);

						my $start_index=index($coordinates[0], "(");

						my $tm_start=substr($coordinates[0],$start_index+1);

						$tm_start=~ s/^\s+|\s+$//g;

						my $tm_end=substr($coordinates[1],0,-1);

						$tm_end=~ s/^\s+|\s+$//g;

						push @tm_regions1 , $tm_start . "-" . $tm_end;

					}

					my $get_regions2=$TMRegions{$toSearch2};
					#print "protein to search " . $toSearch2 ."\n";
					#print " Regions 2"  . " $get_regions2 " . "\n";

					@regions2=split(/,/, $get_regions2);

					foreach my $region (@regions2) {

						my @coordinates=split(/-/, $region);

						my $start_index=index($coordinates[0], "(");

						my $tm_start=substr($coordinates[0],$start_index+1);

						$tm_start=~ s/^\s+|\s+$//g;

						my $tm_end=substr($coordinates[1],0,-1);

						$tm_end=~ s/^\s+|\s+$//g;

						push @tm_regions2 , $tm_start . "-" . $tm_end;

					}

					my $fixed_res=0;

					INNER: foreach my $x ( sort {$a <=> $b} keys %residues_1) { 

						my $tm_check=0;

						foreach my $region ( @tm_regions1) {

							my @coordinates=split(/-/ , $region);

							if ( $x >= $coordinates[0] && $x<= $coordinates[1]) {

								$tm_check=1;
							}
						}

						if ($temp_start == 0){

							if ($tm_check==1) {

								$mol1=$a_a{$residues_1{$x}}; 

							}

							else {

								next INNER ;
							}

							$temp_start=$x;

						}

						else {

							for ( my $y=$temp_start+1 ; $y<$x ; $y++){

								$mol1=$mol1."."; 

							}

							if ( $tm_check ==1) {

								$mol1=$mol1.$a_a{$residues_1{$x}}; 

								$fixed_res++;

							}

							else {

								$mol1=$mol1.".";
							}

							$temp_start=$x; 

						}		

					}
					
					my $origMotif1=$mol1;
					$origMotif1=~s/(\.+)$/''/ge;
					$mol1=~s/\.+$/''/ge;
					$mol1=~s/((\.)\2+)/$2 . "{".length($1)."}"/ge;
					$mol1=~s/\.\{\d+\}$//ge;
					
					

					$temp_start=0;

					$size = keys %residues_2;
					#print  "the size variable 2 is $size\n";
					#print "AA in residues 2: ";
					#print %residues_2;
					#print "\n";

					if ( $size <= 5) {

						print "Motif 2 shorter than 5\n";
						next OUTER;

					}

					if ( $fixed_res >3) {

						#print $pdb_id . "\t" .$mol1 . "\n";
					}

					my $fixed_res1=0;

					INNER2: foreach my $x ( sort {$a <=> $b} keys %residues_2) { 

						my $tm_check=0;

						foreach my $region ( @tm_regions2) {

							my @coordinates=split(/-/ , $region);
							$coordinates[1]=~s/\)$/''/ge;
							#print $x . " ->" . $coordinates[0]. " ->". $coordinates[1]."\n";
							if ( $x >= $coordinates[0] && $x<= $coordinates[1]) {

								$tm_check=1;
							}
						}

						if ($temp_start == 0){

							if ( $tm_check ==1){

								$mol2=$a_a{$residues_2{$x}}; 

							}

							else {

								next INNER2;
							}

							$temp_start=$x;

						}

						else {

							for ( my $y=$temp_start+1 ; $y<$x ; $y++){

								$mol2=$mol2."."; 

							}

							if ($tm_check==1){

								$mol2=$mol2.$a_a{$residues_2{$x}}; 

								$fixed_res1++;

							}

							else {

								$mol2=$mol2.".";
							}

							$temp_start=$x; 

						}

					}

					my $origMotif2=$mol2;

					$origMotif2=~s/\.+$/''/ge;
					$mol2=~s/\.+$/''/ge;
					$mol2=~s/((\.)\2+)/$2 . "{".length($1)."}"/ge;
					$mol2=~s/\.\{\d+\}$//ge;

					if ( $fixed_res > 3 && $fixed_res1 >3) {



						print $pdb_id . "\t" .$mol1 . "\n";
						print $pdb_id . "\t" .$mol2 . "\n";
					

						$chain1=uc($chain1);
						$chain2=uc($chain2);
						$pdb_id=uc($pdb_id);

						my @uniprot_ids1=get_id($pdb_id , $chain1); #sub to get the uniprot id 
	    
	    				my @uniprot_ids2=get_id($pdb_id , $chain2);

	    				#print "Return get pdb id 1 $pdb_id , $chain1: @uniprot_ids1\n";
	    				#print "Return get pdb id 2 $pdb_id , $chain2: @uniprot_ids2\n";

	    				if ( scalar (@uniprot_ids1) > 0  && scalar (@uniprot_ids2) > 0 ) {

							my $sth=$dbh->prepare("SELECT count(*) from motifs WHERE motif=? AND chain=?");

							$sth->bind_param(1,$mol1);
		
							$sth->bind_param(2,$chain1);

							$sth->execute;

							my $count=$sth->fetchrow;
	
							if ( $count == 0 ) {
		
								#print "Motif 1 should be inserted!!!!!!!\n";
								$dbh->do("INSERT INTO motifs ( motif , chain) VALUES ( '$mol1' , '$chain1') ");

								print ORIG $origMotif1 . "\t" . $mol1 . "\n";
					
							}
		
							$sth=$dbh->prepare("SELECT count(*) from motifs WHERE motif=? AND chain=?");

							$sth->bind_param(1,$mol2);
		
							$sth->bind_param(2, $chain2);

							$sth->execute;

							$count=$sth->fetchrow;

							if ( $count == 0 ) {
		
								#print "Motif 2 should be inserted!!!!!!!\n";
								$dbh->do("INSERT INTO motifs ( motif , chain ) VALUES ( '$mol2', '$chain2' ) ");


								print ORIG $origMotif2 . "\t" . $mol2 . "\n";
		
							}
					
						}

						if ( scalar (@uniprot_ids1) == 0 && (@uniprot_ids2) ==0 ) {
			
						print LIST $pdb_id."\t-\t-\t".$chain1."\t". $chain2."\t". $mol1 . "\t" .$mol2."\n";
		
					}
	
					elsif ( scalar(@uniprot_ids1==1) && scalar (@uniprot_ids2 ==1 )) {
		
						print LIST $pdb_id."\t".$uniprot_ids1[0]."\t". $uniprot_ids2[0]. "\t".$chain1. "\t" . $chain2 ."\t". $mol1 . "\t" .$mol2."\n"; #we write it to the file
			
				}
				
					elsif ( scalar ( @uniprot_ids1 ==1) && scalar (@uniprot_ids2 == 0) ){
					
						print LIST $pdb_id."\t". $uniprot_ids1[0] . "\t-\t".$chain1."\t". $chain2."\t". $mol1 . "\t" .$mol2."\n";

						
					}
		
					elsif ( scalar ( @uniprot_ids1 ==0) && scalar (@uniprot_ids2 == 1) ){
					
						print LIST $pdb_id . "\t-\t" . $uniprot_ids2[0] . "\t" .$chain1."\t". $chain2."\t". $mol1 . "\t" .$mol2."\n";
					
					
						
					}
		
		
					else { 
		
						for ( my $x =0 ; $x<scalar(@uniprot_ids1) ; $x++) {
				
							for ( my $y=0; $y<scalar(@uniprot_ids2) ; $y++) {
				
								if ( $uniprot_ids1[$x] eq $uniprot_ids2[$y]) {
					
									next;
					
								}
			
								print LIST $pdb_id."\t".$uniprot_ids1[$x]."\t".$uniprot_ids2[$y]."\t".$chain1. "\t". $chain2 ."\t". $mol1 . "\t" .$mol2."\n"; #we write it to the file

							}
				
						}
			
					}	


					}

					
				}

				%residues_1=();

				%residues_2=();

				$check=0;

				$mol1_start=0;

				$mol1="";

				$mol2_start=0;

				$mol2="";

			}

			if ( $line  =~ /^ *\d/ ) {

				chomp ($line);

				if ( $bond =~ /Non-bonded contacts/) {

					$check=1;

				}

				$line =~ s/ +/_/g;

				my @elements=split("_" , $line);

				if ( substr($line , 0,1) eq "_" ) {

					my $temp=$elements[4];

					$mol1_start=$elements[5]; 

					$residues_1{$mol1_start}=$temp; 

					$temp=$elements[10];

					$mol2_start=$elements[11]; 

					$residues_2{$mol2_start}=$temp; 


				}

				else {

					my $temp=$elements[3];

					$mol1_start=$elements[4]; 

					$residues_1{$mol1_start}=$temp; 	

					$temp=$elements[9];

					$mol2_start=$elements[10]; 

					$residues_2{$mol2_start}=$temp; 

				}



			}

			if ( $line =~ /Salt bridges/ ) {

				last OUTER;

			}

			if ( ($line =~ /bond/) && ($line !~ /Number/) ) {

				$bond=$line;

				chomp ($bond);

			}

		}

	}

}

	#we get all pdbs from OPM found in our database, and for each we get its PDBsum file.

	open(FILEPDBSums , "listPDBSums.txt" ); # this file contains the list of pdb sums from OPM validated only transmembrane

	open(WRITE , ">pdbsForPBSum.txt");

	my $counter=0;

	my @pdbs=();

	while (<FILEPDBSums>) {

		chomp $_;

		#my $sth= $dbh->prepare("SELECT distinct pdb_id from pdbs where pdb_id=? and valid=1");
		my $sth= $dbh->prepare("SELECT distinct pdb_id from pdbs where pdb_id=? ");

		$sth->bind_param(1, $_);

		$sth->execute;

		while ( my @row=$sth->fetchrow_array) {

			my $line=join ( "\t" , @row);

			$counter++;

			print WRITE $line ."\n";

			push @pdbs, $row[0];

		}
	}

	#print $counter . "\n";

	close FILEPDBSums;

	close WRITE;

	open ( FILETMSegments , "TMsegments.txt");

	open ( TEST , ">testing2.txt");

	open ( MAPPING , "pdb_chain_uniprot.tsv");

	open ( WRITE , ">pdbs_with_pdbsum.txt");

	my $mappingVersion = <MAPPING>;
	my $header=<MAPPING>;

	my @arrayOfRegions=();

	while (<MAPPING>) {

		chomp $_;

		@arrayOfRegions=();

		my @line_split=split(/\t/,$_);

		my $pdb_chain=$line_split[0]." ". $line_split[1];

		my $regions=$line_split[5] . "\t" . $line_split[6] . "\t" . $line_split[7]. "\t" .$line_split[8];

		push @{$mapping{$pdb_chain}} , $regions;

	}

	my @all=<FILETMSegments>;

	my %tmNeeded=();

	my $counter2=0;

	foreach my $pdb (@pdbs) {

		
		my @matched=grep /$pdb/i, @all;

		if (scalar(@matched) == 0) {

			next;
		}

		#elsif (scalar(@matched) == 1) {
		#	print @matched;
		#	next;

		#}

		else {
		
			my @chains=();

			my @getPDB=();

			my @getChains=();

			my @getRegions=();

			foreach my $line (@matched) {
				#print $line;

				chomp $line;

				@getRegions=split(/Segments:/ , $line);

				@getChains=split(/-/ , $getRegions[0]);

				@getPDB=split(/ / , $getChains[0]);

				my $pdb_chain=  $getPDB[0]." ".  $getPDB[1];

				push @chains , $getPDB[1];

				#print "checking the pdb_chain " . $pdb_chain ."\n";
				$TMRegions{$pdb_chain}=$getRegions[1];

			}

			getPDBsumFile($getPDB[0], @chains);

			$counter2++;

			$tmNeeded{$getChains[0]}=$getRegions[1];

		}

	}

	checkCorrectPDBsum();

	getMotifs();

	sub get_id {

#this subroutine takes as paramaters the pdb id and its chain and gets the uniprot id
#that corresponds to it from our database

    my ($pdb_id , $chain) = @_;

    my @ids;
    
    my $sth=$dbh->prepare("SELECT DISTINCT p_id , chain FROM pdbs WHERE pdb_id=? ");
    
    $sth->bind_param(1, $pdb_id);
    
    $sth->execute;  
    
    SEARCH:while ( my @row=$sth->fetchrow_array) {
    
    	my @chains_split=split(/\//, $row[1]);
    	
    	for ( my $x=0; $x<scalar(@chains_split) ; $x++) {
    	
    		if ( $chain eq $chains_split[$x] ) {
    		
    			#print "GETTING THE IDS" . $row[0] ."\n";
    		
   				push @ids , $row[0];
   		
   			}
    	
    	}
        
    }

    return @ids; 
	
}
