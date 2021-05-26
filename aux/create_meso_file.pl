
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

open(FILE,$opts{a}) or die "cannot open file $opts{a}\n";
my $hash=();
while(my $line=<FILE>){
    chomp $line;
	$line=~s/:##/ /g;
	$line=~s/=/ /g;
	
        my @data=split / /,$line;
	#print Dumper(@data);
	push(@{$hash->{$data[0]}->{$data[1]}},$data[2]);
}

my $bhash=load_cram($opts{b});
#print Dumper($hash);
print join("\t","id","vcf","normal_cram","normal_id","tumor_id")."\n";
foreach my $vcf (sort keys %{$hash}){
	my $s=$vcf;
	$s=~s/_filtered_PASS_norm.vcf.gz//g;
	my $normal=$hash->{$vcf}->{"normal_sample"}[0];
	print join("\t",$s,$vcf,$bhash->{$normal},$normal,join(",",@{$hash->{$vcf}->{"tumor_sample"}}))."\n";
}

#we load the meso data
sub load_cram{
	my ($file)=@_;
	open(FILE,$file) or die "cannot open file $file\n";
	my $hash=();
	while(my $line=<FILE>){
	        chomp $line;
        	my @data=split /_/,$line;
		$hash->{$data[0]}=$line;
	}
	return $hash;
}
