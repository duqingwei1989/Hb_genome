#!/usr/bin/perl -w
use Getopt::Long;

my $help;
my $vcf;
my $sample_in;
my $sample_ex;
my $out;
my $missing_rate;

my $usage=<<USAGE;
#########################################################
Note:this script used to convert vcf to fasta;
Usage: perl $0 -v snp.vcf -o snp.fasta 

        -h || --help       usage;
        -v || --vcf        snp.vcf, input
        -o || --out        snp.fasta, output
        -s || --sample_in  sample included (options)
        -e || --sample_ex  sample excluded (options)
        -m || --missing_rate missing rate, default 0.4

#########################################################
USAGE
#
#
GetOptions("help|h"=>\$help,
           "vcf|v=s"=>\$vcf,
           "sample_in|s=s"=>\$sample_in,
           "sample_ex|e=s"=>\$sample_ex,
           "out|o=s"=>\$out);

die"$usage" if($help||(!defined $vcf)||(!defined $out));

$missing_rate = (defined $missing_rate)?$missing_rate:"0.4";
my $sam_in = (defined $sample_in)?$sample_in:"NONE";
my $sam_ex = (defined $sample_ex)?$sample_ex:"NONE";

my %indb=();
my %exdb=();
if($sam_ex ne "NONE"){
    open(IN, $sam_ex) or die"";
    while(<IN>){
        chomp;
        my $sample = (split/\s+/,$_)[0];
        $exdb{$sample}++;
        }
    close IN;
    }
    
if($sam_in ne "NONE"){
    open(IN, $sam_in) or die"";
    while(<IN>){
        chomp;
        my $sample = (split/\s+/,$_)[0];
        $indb{$sample}++;
      }
    close IN;
    }

my $infordb=();
my $namedb=();
open(IN, "grep -v '##' $vcf |") or die"";
my $head = <IN>;
my @headb = split(/\s+/,$head);
foreach my $i(9..$#headb){
    $namedb{$i} = $headb[$i];
    }
while(<IN>){
    chomp;
    my @data = split(/\s+/,$_);
    my $chrn = $data[0];
    my $posi = $data[1];
    my $refN = $data[3];
    my $altN = $data[4];
    next if((length $refN)>1 or (length $altN)>1);
    my $mr   = 0;
    my $count = 0;
    foreach my $i(9..$#data){
        my $sample = $namedb{$i};
        next if($sam_ex ne "NONE" and exists($exdb{$sample}));
        next if($sam_in ne "NONE" and !exists($indb{$sample}));
        $count++;
        $mr++ if($data[$i] eq "./.");
        }
    $mr = $mr/$count; $mr = sprintf("%.2f",$mr);
    next if($mr>$missing_rate);
    foreach my $i(9..$#data){
        my $sample = $namedb{$i};
        my $geno   = $data[$i];
        my $baseN;
        if($geno=~/0\/0/){
            $baseN = $refN;
        }elsif(($geno=~/0\/1/)or($geno =~ /1\/0/)){
             my @infodb = split(/:/,$geno);
             my $tmp    = $infodb[0];
             my ($a,$b) = split(/,/,$infodb[1]) if($tmp ne "./.");
             $baseN = ($a>=$b)?$refN:$altN;
        }elsif($geno=~/1\/1/){
            $baseN = $altN;
        }else{
            $baseN = '-';
            }
        next if($sam_ex ne "NONE" and exists($exdb{$sample}));
        next if($sam_in ne "NONE" and !exists($indb{$sample}));
        $infordb{$i}->{'sample'} = $sample;
        $infordb{$i}->{'seq'}   .= $baseN;
        }
    }
close IN;

open(OUT, "> $out") or die"";
foreach my $i (sort {$a<=>$b} keys %infordb){
    my $sample = $infordb{$i}->{'sample'};
    my $seq    = $infordb{$i}->{'seq'};
    print OUT ">$sample\n$seq\n";
    }
close OUT;
