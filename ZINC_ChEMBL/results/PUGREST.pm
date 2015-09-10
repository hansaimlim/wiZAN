#!/usr/bin/perl

package PUGREST;
use LWP::Simple;
use IDMAP;

require Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(get_InChIKey_by_SMILES get_InChIKey_by_name get_InChIKey_by_compound get_CID_by_substance get_InChIKey_by_CID);

#--------------------------TEST AREA------------------------

#--------------------------TEST AREA------------------------

sub get_InChIKey_by_CID
{
        #input PubChem CID then output InChIKey
        #one CID gives only one InChIKey
	my $cid = shift @_;
        my $url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/$cid/property/InChIKey/TXT";
        my $inchikey = get $url;
        chomp($inchikey);
	if (defined($inchikey)){
		return $inchikey;
	} else {
		return;
	}
}
sub get_InChIKey_by_SMILES
{
	#input SMILES (non-canonical)
	#output InChIKey
	my $smiles = shift @_;
	my $url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/$smiles/property/InChIKey/TXT";
	my $inchikey = get $url;
	chomp($inchikeu);
	if (defined($inchikey)){
		return $inchikey;
	} else {
		return;
	}
}
1;
