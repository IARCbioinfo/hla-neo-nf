
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
   print "$0 usage : -a prefix.all_epitopes.tsv   -p prefix\n";
   print "Error in use\n";
   exit 1;
}

my %opts = ();
getopts( "a:p:", \%opts );
if ( !defined $opts{a} or !defined $opts{p}) {
   usage;
}

# to fill when using only genomic data 
#$VAR30 = 'Tumor RNA Depth';
#$VAR31 = 'Tumor RNA VAF';
#$VAR34 = 'Gene Expression';
#$VAR35 = 'Transcript Expression';
my $hcols={"Tumor RNA Depth"=>"11","Tumor RNA VAF"=>"0.26","Gene Expression"=>"2.0","Transcript Expression"=>"2.0"};
=bla
--normal-cov NORMAL_COV
                        Normal Coverage Cutoff. Sites above this cutoff will
                        be considered. (default: 5)
  --tdna-cov TDNA_COV   Tumor DNA Coverage Cutoff. Sites above this cutoff
                        will be considered. (default: 10)
  --trna-cov TRNA_COV   Tumor RNA Coverage Cutoff. Sites above this cutoff
                        will be considered. (default: 10)
  --normal-vaf NORMAL_VAF
                        Normal VAF Cutoff. Sites BELOW this cutoff in normal
                        will be considered. (default: 0.02)
  --tdna-vaf TDNA_VAF   Tumor DNA VAF Cutoff. Sites above this cutoff will be
                        considered. (default: 0.25)
  --trna-vaf TRNA_VAF   Tumor RNA VAF Cutoff. Sites above this cutoff will be
                        considered. (default: 0.25)
  --expn-val EXPN_VAL   Gene and Transcript Expression cutoff. Sites above
                        this cutoff will be considered. (default: 1.0)

=cut


open(FILE,$opts{a}) or die "cannot open file $opts{a}\n";
my $h=<FILE>; #we read the header
chomp $h;
my @colnames=split("\t",$h);
my $h2cols=();
for(my $i=0; $i<scalar(@colnames); $i++){
	if(defined $hcols->{$colnames[$i]}){
		$h2cols->{$colnames[$i]}=$i;
	}
}
#print Dumper($h2cols);
open(OUT,">$opts{p}.all_epitopes_rnafill.tsv") or die "cannot open file $opts{p}.all_epitopes_rnafill.tsv\n";
print OUT join("\t",@colnames)."\n";
#we continue reading the file
while(my $line=<FILE>){
     chomp $line;
     my @data=split /\t/,$line;
     #print Dumper(@data);
     foreach my $k(keys %{$hcols}){
	#if the data is marked as missing
	if($data[$h2cols->{$k}] eq "NA"){
		$data[$h2cols->{$k}]=$hcols->{$k};
	}
     }		
     print OUT join("\t",@data)."\n";
}

#we run the pvacseq default filters which include
#The binding filter is used to remove neoantigen candidates that do not meet desired peptide:MHC binding criteria. The coverage filter is used to remove variants that do not meet desired read count and VAF criteria (in normal DNA and tumor DNA/RNA). The transcript support level filter is used to remove variant annotations based on low quality transcript annotations. The top score filter is used to select the most promising peptide candidate for each variant. Multiple candidate peptides from a single somatic variant can be caused by multiple peptide lengths, registers, HLA alleles, and transcript annotations.
my $f1="pvacseq binding_filter $opts{p}.all_epitopes_rnafill.tsv $opts{p}.all_epitopes_binding.tsv";
my $f2="pvacseq coverage_filter  --tdna-vaf 0.1 $opts{p}.all_epitopes_binding.tsv $opts{p}.all_epitopes_coverage.tsv";
my $f3="pvacseq transcript_support_level_filter $opts{p}.all_epitopes_coverage.tsv $opts{p}.all_epitopes_tsl.tsv";
my $f4="pvacseq top_score_filter $opts{p}.all_epitopes_tsl.tsv $opts{p}.filtered2.tsv";
#we agregate the filters
#my $agg="pvacseq generate_aggregated_report $opts{p}.filtered.tsv $opts{p}.filtered.aggregated.tsv";
system($f1) == 0 or die  "error excecuting : $?";
system($f2) == 0 or die "error excecuting : $?";
system($f3) == 0 or die "error excecuting : $?";
system($f4) == 0 or die "error excecuting : $?";
#system($agg) == 0 or die  "error excecuting : $?";
