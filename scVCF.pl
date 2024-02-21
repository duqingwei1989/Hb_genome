#!/usr/bin/perl -w
use Getopt::Long;

my $help;
my $vcf;
my $bed;

my $usage=<<USAGE;
####################################################################
Note:this script used to extract single copy vcf;       
Usage: perl $0 -v vcf_file -b sc_bed -o clean_vcf 
                                                         
        -h || --help       usage;                        
        -v || --vcf        snp.vcf, input
        -b || --bed        sc_bed, input                
        -o || --out        clean vcf, output       
							  
####################################################################
USAGE
#
#
GetOptions("help|h"=>\$help,
           "vcf|v=s"=>\$vcf,
           "bed|b=s"=>\$bed,
           "out|o=s"=>\$out);

die"$usage" if($help||(!defined $vcf)||(!defined $bed)||(!defined $out));


my $scdb=();
open(IN,$bed) or die"";
while(<IN>){
    chomp;
    my @data = split(/\s+/,$_);
    foreach my $i($data[1]..$data[2]){
        my $keys     = $data[0]."_".$i;
        $scdb{$keys}++;
    } 
}
close IN;

open(OUT,">$out") or die"";
if($vcf=~/\.gz/){
	open(IN,"gunzip -c $vcf|")or die"";
}else{
	open(IN,$vcf)or die"";
}
#my $header = <IN>;
#chomp $header;
#print OUT "$header\n";
while(<IN>){
	chomp;
	if(/##|#CHROM/){
		print OUT "$_\n";	
	}else{
	my @data =split(/\s+/,$_);
	my $l1   = length $data[3];
	my $l2   = length $data[4];
	next if($l1>1);
	next if($l2>1);
	my $posi = $data[0]."_".$data[1];
	print OUT "$_\n" if(exists($scdb{$posi}));    
	}
}
close IN;
close OUT;
