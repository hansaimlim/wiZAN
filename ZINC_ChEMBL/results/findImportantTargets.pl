#!/usr/bin/perl
use strict;
use warnings;

#important targets: GPCRs and Kinases
#input file must have at least 2 columns: Chemical	Target(Accession)

die "Usage: $0 <chem-prot file>\n" unless @ARGV == 1;
my $input = shift @ARGV;

my $gpcrIn='./ChEMBL_GPCRs.tsv';
my $gpcrOut='./GPCRs.tsv';
my $kinaseIn='./ChEMBL_Kinases.tsv';
my $kinaseOut='./Kinases.tsv';

my %gpcr;
my %kinase;
open my $G, '<', $gpcrIn;
while (my $line = <$G>){
	chomp($line);
	next if $. == 1;
	my @words = split(/\t+/, $line);
	my $acc = $words[3];	#4th column is the protein Accession
	$gpcr{$acc} = 1;
}
close $G;
open my $K, '<', $kinaseIn;
while (my $line = <$K>){
	chomp($line);
	next if $. == 1;
	my @words = split(/\t+/, $line);
	my $acc = $words[3];	#4th column is the protein Accession
	$kinase{$acc} = 1;
}
close $K;

open my $Gout, '>', $gpcrOut;
open my $Kout, '>', $kinaseOut;
open my $in, '<', $input;
while (my $line = <$in>){
	chomp($line);
	next if $. == 1;
	my @words = split(/\s+/, $line);	#white-space delimited
	my $acc = $words[1];	#2nd column is the protein Accession
	if ($gpcr{$acc}){
		print $Gout $line, "\n";
	} elsif ($kinase{$acc}){
		print $Kout $line, "\n";
	}
}
close $in;
close $Kout;
close $Gout;
