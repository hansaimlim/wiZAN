#!/usr/bin/perl
use strict;
use warnings;

my $file='drug_infor.csv';
open my $F, '<', $file;
while (my $line=<$F>){
        chomp($line);
        my $inchi='';
        if ($line =~ m/("InChI=.*")/){
                $inchi=$1;
                $line =~ s/,("InChI=.*")//;
        }
	$line =~ s/,,/,Unknown,/g;
        my @words = split(/,/, $line);
        my $ikey = pop @words;  #last one
        my $smile= pop @words;  #may be empty
#	my $dump1= pop @words;
#	my $dump2= pop @words;
#	my $dump3= pop @words;
#	my $dump4= pop @words;
#	my $dump5= pop @words;

	
	foreach my $entry (@words){
		$entry =~ s/\s+//g;
	}
	my $num=shift @words;
	my $dname=shift @words;
	my $bname=shift @words;

	my $target='';
	while (@words){
		my $entry = shift @words;
		if ($entry =~ m/^[\"A-Z0-9\/]+$/){
			next if $entry =~ m/^[0-9]+$/;
			if ($entry =~ m/^\"/){
				my $output=$entry;
				until ($output =~ m/\"$/){
					my $add=shift @words;
					$output.="\t".$add;
				}
				$target=$output;
				last;
			} else {
				$target=$entry;
				last;
			}
		}
	}
     	if ($ikey =~ m/[A-Z]{14}-[A-Z]{10}-[A-Z]{1}/){
		print $ikey, "\t", $target, "\n";
        }
}
close $F;
