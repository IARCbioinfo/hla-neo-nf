
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

my ($pvac)=load_pvac($opts{a});

open(FILE,$opts{b}) or die "cannot open file $opts{b}\n";
while(my $line=<FILE>){
    chomp $line;
        #my @data=split /\t/,$line;
   if(defined $pvac->{$line}){
	print join(" ",$line, $line)."\n";
  }elsif(defined $pvac->{"HLA-".$line}){
	 print join(" ",$line, "HLA-".$line)."\n";
    }else{
	 print $line." NOT_PVACTOOLS\n";
    }
}

sub load_pvac{
	my ($file)=@_;
	my $hash=();
	open(FILE,$file) or die "cannot open file $file\n";

	while(my $line=<FILE>){
    	chomp $line;
		$hash->{$line}=1;
	}
	return ($hash);
}
