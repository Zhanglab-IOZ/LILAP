#!/usr/bin/perl -w
use strict;

my $in=$ARGV[0];
my $out=$ARGV[1];
my $tag=$ARGV[2];

open FILE,"$in";
my @name=split/\./,$in;
my $all=0;
my %all=();
while(<FILE>){
    chomp $_;
    my @line=split/\t/,$_;
    my $num=int($line[0]*100);
    $all{$num}{gc}++;
    $all{$num}{cov}.=$line[1].",";
    $all++;
}
close FILE;

open OUT,">$out";
for my $gc(sort {$a<=>$b} keys %all){
    my @cov=split/,/,$all{$gc}{cov};
    my @median=sort {$a<=>$b} @cov;
    my $num=scalar(@median);
    my $med= int(($num-1)/2);
    print OUT "$gc\t$median[$med]\t$name[0]\t$tag\t$name[$#name-1]\n";
}

for my $gc(sort {$a<=>$b} keys %all){
    $all{$gc}{gc}=$all{$gc}{gc}/$all*5;
    print OUT "$gc\t$all{$gc}{gc}\tgc\t$tag\t$name[$#name-1]\n";
}
close OUT;