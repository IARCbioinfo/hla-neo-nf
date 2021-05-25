
###############################################################################
# Author: Alex Di Genova
# Laboratory: IARC/SCG/RCG
# Copyright (c)
# year: 2020
###############################################################################
use Data::Dumper;
use Getopt::Std;
use strict;

sub usage {
   print "$0 usage : -a <hla_regions> -b <BAM/CRAM file> -r <chr6.fa> -p <prefix>\n";
   print "Error in use\n";
   exit 1;
}

my %opts = ();
getopts( "a:b:r:p:", \%opts );
if (!defined $opts{a} or !defined $opts{b} or !defined $opts{r} or !defined $opts{p}) {
   usage;
}

my $hla_regs=load_hla_reg($opts{a});

open(HLASEQS,">$opts{p}.hla.seqs.sam") or die "cannot create $opts{p}.hla.seqs.sam\n";
open(HEADER,"samtools view -H $opts{b} |") or die "Getting BAM/CRAM header failed\n";
#we add the header to the SAM file
while(my $line=<HEADER>){
	print HLASEQS $line;
}

foreach my $r(@{$hla_regs}){
	#print $r."\n";
	#we get the reads using samtools
	open(READS,"samtools view $opts{b} $r |") or die "cannot open BAM/CRAM file\n";
	while(my $line=<READS>){
		print HLASEQS $line;
	}
}

#convert bam to fastq
system("samtools view -b $opts{p}.hla.seqs.sam | sambamba sort -p -n -o - /dev/stdin | bamToFastq -i /dev/stdin -fq \"$opts{p}.hla.fwd.fq\" -fq2 \"$opts{p}.hla.rev.fq\" 2>bam2fastq.log.err");
#we run BWA-MEM on the indexed
system("bwa mem -t2 $opts{r} $opts{p}.hla.fwd.fq $opts{p}.hla.rev.fq 2>$opts{p}.bwa.log | samtools view -b - | samtools sort -T STMP -o $opts{p}.hla.chr6.bam -");
system("samtools index $opts{p}.hla.chr6.bam");
system("samtools view -o $opts{p}.mhc.bam $opts{p}.hla.chr6.bam chr6:29844528-33100696");
system("samtools index $opts{p}.mhc.bam");

sub load_hla_reg{
	my ($file)=@_;
	open(HLAR,$file) or die "cannot open file";
	my $hla_r=();
	while(my $line=<HLAR>){
		chomp $line;
		my @data=split("\t",$line);
		push(@{$hla_r},$data[0].":".$data[1]."-".$data[2]);
	}
	return $hla_r;
}
