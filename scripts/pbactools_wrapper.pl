
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
   print "$0 usage : -a xHLA.json -b xHLA2pvac -c normal_id -d vep_vcf -t tumor_id_name -p prefix -e tools\n";
   print "Error in use\n";
   exit 1;
}

my %opts = ();
getopts( "a:b:c:d:t:p:e:", \%opts );
if ( !defined $opts{a} or !defined $opts{b} or !defined $opts{c} or !defined $opts{d} or !defined $opts{t} or !defined $opts{p}) {
   usage;
}


my ($alleles)=parse_xHLA_output($opts{a},$opts{b});
my $tools=$opts{e};
$tools=~s/,/ /g;

# we run pvacseq
my @tumors=split(",",$opts{t});
my $i=1;

foreach my $t(@tumors){
my $cmd="";
if(scalar(@tumors) == 1){
  $cmd="pvacseq run   --n-threads 2 --iedb-install-directory /opt/iedb --pass-only --normal-sample-name $opts{c} $opts{d} $t ".join(",",@{$alleles->{pvac}});
   $cmd.=" $tools $opts{p}_T_pvactools";
   print $cmd."\n";
 }else{
     $cmd="pvacseq run --n-threads 2 --iedb-install-directory /opt/iedb --pass-only --normal-sample-name $opts{c} $opts{d} $t ".join(",",@{$alleles->{pvac}});
      $cmd.=" $tools $opts{p}_T".$i."_pvactools";
      print $cmd."\n";
      $i++;
 }
 system($cmd)==0  or die "pvactools:  failed: $?";
}

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
