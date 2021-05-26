
###############################################################################
# Author: Alex Di Genova 
# Laboratory: ERABLE/INRIA
# Copyright (c)
# year: 2017
###############################################################################
use Data::Dumper;
use Getopt::Std;
use strict;

sub usage {
   print "$0 usage : -a  -b  -c\n";
   print "Error in use\n";
   exit 1;
}

my %opts = ();
getopts( "a:b:c:", \%opts );
if ( !defined $opts{a}  ) {
   usage;
}


my ($alleles)=parse_xHLA_output($opts{a},$opts{b});








sub parse_xHLA_output{
	my ($file,$xhla2pvac)=@_;
	my $hash=load_xHLA2pvactools($xhla2pvac);	
	open(FILE,$file) or die "cannot open file $file\n";
	my $alleles=();
	while(my $line=<FILE>){
	    chomp $line;
            my @data=split (" ",$line);
	    #print Dumper(@data);
	    if($data[0]=~m/\*/){
		$data[0]=~s/\"//g;
		$data[0]=~s/\,//g;
		push(@{$alleles->{xHLA}},$data[0]);
		push(@{$alleles->{pvac}},$hash->{$data[0]});
		if($hash->{$data[0]} eq "NOT_PVACTOOLS"){
			print STDERR "$data[0] allele not valid for PVACTOOLS\n";
		}
		print join(" ",$data[0],$hash->{$data[0]})."\n";
	    }
	}
	return ($alleles);
}

sub load_xHLA2pvactools{
	my ($file)=@_;
	open(FILE,$file) or die "cannot open file $file\n";
	my $hash=();
	while(my $line=<FILE>){
    		chomp $line;
	        my @data=split (" ",$line);
		$hash->{$data[0]}=$data[1];
	}
	close(FILE);
	return ($hash);
}
