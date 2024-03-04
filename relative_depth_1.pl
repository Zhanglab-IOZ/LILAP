#!/usr/bin/perl -w
use strict;

my %all=();
open FILE,"C01.ccs.cut19bp.fq.dm6.F904.s.depth";
while(<FILE>){
    chomp $_;
    my @line=split/\t/,$_;
    my $num=int($line[1]/100000);
    $all{$line[0]}{$num}{r}+=$line[2];
    $all{$line[0]}{$num}{c}++;
}
close FILE;

open OUT,">A4.6.9G.100000.depth1";
for my $chr(sort keys %all){
    for my $num(sort {$a<=>$b} keys %{$all{$chr}}){
        my $cov=$all{$chr}{$num}{r}/$all{$chr}{$num}{c};
        #my $cov1=$all{$chr}{$num}{r}/10000;

        print OUT "$chr\t$num\t$cov\tall\n";
    }
}
