#!/usr/bin/perl
use strict;
use warnings;
use Statistics::Basic qw(:all);

die "Usapge: $0 <matlab output file>\n" unless @ARGV == 1;
my $infile = shift @ARGV;
open my $IN, '<', $infile;
my @times;
while (my $line = <$IN>){
	chomp($line);
	if ($line =~ m/Elapsed time/){
		$line =~ m/Elapsed time is (\d+\.\d*) seconds\./;
		my $time = $1;
		push @times, $time;
	} elsif ($line =~ m/>> Rank/){
		$line =~ m/>> Rank = (\d+), iteration = (\d+), rowsize = (\d+), colsize = (\d+)/;
		my $rank = $1;
		my $iter = $2;
		my $rows = $3;
		my $cols = $4;
		print "rank $rank, iter $iter, rows $rows, cols $cols\n";
	} 
}
close $IN;
my $mean = mean(@times);
my $var = var(@times);
my $sem = sqrt($var/(scalar @times));
print "mean $mean, s.e.m. $sem\n";
