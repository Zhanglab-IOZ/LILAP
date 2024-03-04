#!/usr/bin/perl -w
use strict;

for my $file(glob "*.depth2"){
    open FILE,"$file";
    my @median=0;
    while(<FILE>){
	unless($_=~/chrX/ || $_=~/chrY/){
	    my @line=split/\t/,$_;
	    push (@median, $line[2]);
	}
    }
    close FILE;
    @median=sort {$a<=>$b} @median;
    my $num=scalar(@median);
    my $med= int(($num-1)/2);
    my $median1=$median[$med]/2;
    
    my $out=substr($file,0,length($file)-1)."3";
    open OUT,">$out";
    open FILE,"$file";
    while(<FILE>){
	chomp $_;
	my @line=split/\t/,$_;
	$line[2]=$line[2]/$median1;
	print OUT "$line[0]\t$line[1]\t$line[2]\t$line[3]\t$line[4]\n";
    }
    close OUT;
    close FILE;
}
