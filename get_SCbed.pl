#!/usr/bin/perl -w
die"Usage:perl $0 blast_out gff_file out_bed\n" if((!defined $ARGV[0])or(!defined $ARGV[1]));
my %scgenedb;
my %countdb;
open(IN, $ARGV[0]) or die"";
while(<IN>){
    chomp;
    my @data = split(/\s+/,$_);
    $countdb{$data[0]}++;
    
    }
close IN;
foreach my $fhgene (keys %countdb){
    next if($countdb{$fhgene}!=1);
    $scgenedb{$fhgene}++;
    }

my %single_basedb = ();
open(OUT, ">$ARGV[2]") or die"";
open(IN, "grep 'gene' $ARGV[1] |") or die"";
while(<IN>){
    chomp;
    my @data = split(/\s+/,$_);
    my $gene = $1 if(/ID=(HorSp\w+)/);  ##不同基因组改这里
    next if(!exists($scgenedb{$gene}));
    print OUT "$data[0] $data[3]    $data[4]    $gene\n";
#   foreach my $i($data[3]..$data[4]){
#       my $key = $data[0]."_".$i;
#       $single_basedb{$key} = $gene;
#       }
    
    }
close IN;
close OUT;
