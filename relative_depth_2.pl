#!/usr/bin/perl -w
use strict;

my %all=();
my $count=0;
open FILE,"C01.ccs.cut19bp.fq.dm6.F904.s.100000.depth1";
while(<FILE>){
    chomp $_;
    my @line=split/\t/,$_;
    if($line[0] eq 'X'){
	$line[0]=23;
    }elsif($line[0] eq 'Y'){
	$line[0]=24;
    }
    $all{$line[0]}{$line[1]}=$line[2];
}
close FILE;

open OUT,">A4.6.9G.100000.depth2";
for my $chr(sort {$a<=>$b} keys %all){
    for my $pos(sort {$a <=>$b} keys %{$all{$chr}}){
	$count++;
	if(length($chr)==1){
	    print OUT "chr$chr\t";
	}else{
	    if($chr eq '24'){
		print OUT "chrY\t";
	    }elsif($chr eq '23'){
		print OUT "chrX\t";
	    }else{
		print OUT "chr$chr\t";
	    }
	    
	}
	#}elsif($chr eq '24'){
	#    print OUT "chrY\t";
	#}else{
	#    print OUT "chr$chr\t";
	#}
	print OUT "$pos\t$all{$chr}{$pos}\t$count\n";
    }
}
close OUT;
