#!/usr/bin/perl -w
use strict;

my $in=$ARGV[0];
my $out=$ARGV[1];
my $tag=$ARGV[2];

my %all1=();
my $name='';
open FILE,"dmel.ref.fasta";
while(<FILE>){
    chomp $_;
    if($_=~/^>/){
	$name=substr($_,1,);
    }else{
	$all1{$name}.=$_;
    }
}
close FILE;

my $count=0;
for (my $i=500;$i<=1000;$i=$i*10){
    my %all=();
    open FILE,"$in";
    while(<FILE>){
	unless($_=~/mito/ || $_=~/X/ || $_=~/Y/){
	    $count++;
	    chomp $_;
	    my @line=split/\t/,$_;
	    my $num=int($count/$i);
	    $all{$line[0]}{$num}{r}+=$line[2];
	    my $base=substr($all1{$line[0]},$line[1]-1,1);
	    if($base=~/[GCgc]/){
		$all{$line[0]}{$num}{gc}++;
	    }
	    $all{$line[0]}{$num}{c}++;
	}
    }
    close FILE;
    
    my @depth=();
    my $median=0;
    my @gc=();
    for my $chr(sort keys %all){
	for my $num(sort {$a<=>$b} keys %{$all{$chr}}){
	    if($all{$chr}{$num}{c} == $i){
		my $cov=$all{$chr}{$num}{r}/$all{$chr}{$num}{c};
		if(!defined $all{$chr}{$num}{gc}){
		    $all{$chr}{$num}{gc}=0;
		}
		my $gc=$all{$chr}{$num}{gc}/$all{$chr}{$num}{c};
		push(@depth,$cov);
		push(@gc,$gc);
	    }
	}
    }
    my @median=sort {$a<=>$b} @depth;
    my $num=scalar(@median);
    my $med= int(($num-1)/2);
    $median=$median[$med];
    my $out1=$out.".".$i.".dgc";
    open OUT,">>$out1";
    for(my $j=0;$j<=$#depth;$j++){
	my $log=$depth[$j]/$median;
	print OUT "$gc[$j]\t$log\tC01\t$tag\t$i\n";
    }
}
close OUT;
