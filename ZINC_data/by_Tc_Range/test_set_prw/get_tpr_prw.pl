#!/usr/bin/perl
use strict;
use warnings;

#This is NOT for 10fold cross validation tests

#This script is to calculate TPR from PRW output files.
#The PRW output files must be in .dat extension

#There must be NO other .dat files except for PRW output in the directory
#There must be wiZAN formatted test files(.txt extension) in the PRW output directory.

# # # # # # # # # #
#The test file and prw output file must be in the same order of chemicals (line number should match)
#The prw output file must contain protein names with the predicted scores
# # # # # # # # # #

die "Usage: $0 <PRW test file> <PRW result file> <output file>\n" unless @ARGV == 3;
my $testfile = shift @ARGV;
my $resultfile = shift @ARGV;
my $output_file = shift @ARGV;
chomp($output_file);

my $condition_positive = 0;	#total number of interaction pairs on test
my %rank_count;	#keys: cutoff ranks, values: true positive counts at the given rank (not cumulative)
my %linenum_protein;	#keys: line numbers, values: protein names
open my $TEST, '<', $testfile or die "Could not open test file $testfile: $!\n";
while (my $line = <$TEST>){
	chomp($line);
	my @words = split(/\t/, $line);
	my $chem = shift @words;
	my $prot = shift @words;
	$linenum_protein{$.} = $prot;
}
close $TEST;
$condition_positive += scalar(keys %linenum_protein);	#add up the number of tested pairs
print "$condition_positive Condition Positives\n";
	
open my $PRW, '<', $resultfile or die "Could not open prw file $resultfile: $!\n";
while (my $line = <$PRW>){
	#each line represents test result for one chemical-protein pair
	#looking for one single protein for each line

	chomp($line);
	my @words = split(/\s+/, $line);
	my $chem = shift @words;
	my $true_prot = $linenum_protein{$.};	#protein name for the given line(chemical)
	my $index = 0;
	my $temp_score_index;
	my $temp_score = 0;
	my $temp_rank;
	RANKING: while($index <= 398){
		#must be stopped when the true protein is found
		#get temporary rank
		if ($words[$index] eq $true_prot){
			#true protein found at the index
			$temp_score_index = $index + 1;
			$temp_score = $words[$index+1];
			last RANKING;
		} else {
			$index++;
		}
	}
	next if $temp_score == 0;	#temp score remains 0 if not found under rank 200  (index 398)
	CUMRANK: while($words[$temp_score_index] == $temp_score){
		#must count more for proteins with same prediction scores (cumulative rank)
		$temp_score_index += 2;	#look up the next score if current score is the same
	}
	my $rank = (($temp_score_index - 1)/2) + 1;	#each protein name is followed by its predicted score
	$rank_count{$rank}++;
}
close $PRW;
open my $OUT, '>', $output_file;
print $OUT "Cutoff_rank\tTPR\n";
for (my $cutoff=1; $cutoff<=200; $cutoff++){
	my $true_positives = 0;
	foreach my $rank (keys %rank_count){
		if ($rank <= $cutoff){
			$true_positives += $rank_count{$rank};
		}
	}
	print $OUT $cutoff, "\t", $true_positives/$condition_positive, "\n";
}
close $OUT;

