#!/usr/bin/perl -w
use strict;

open FILE,"C01.ccs.fasta";
my %all1=();
my %all2=();
while(<FILE>){
    chomp $_;
    unless($_=~/^>/){
	if(length($_)>=150){
	my $fas1=substr($_,0,100);
	my $fas2=substr($_,-1,100);
	my $fas3=reverse $fas2;
	$fas3=~tr/ATCG/TAGC/;
	for(my $i=0;$i<length($fas1);$i++){
	    $all1{$i}.=substr($fas1,$i,1);
	}
	for(my $i=0;$i<length($fas3);$i++){
	    $all1{$i}.=substr($fas3,$i,1);
	}
	}
    }
}
close FILE;

open OUT,">insert_bias_C01";
for my $pos( sort {$a<=>$b} keys %all1){
    my $len=length($all1{$pos});
    my $g=$all1{$pos};
    $g=~s/[^G]//g;
    $g=length($g)/$len;
    
    my $c=$all1{$pos};
    $c=~s/[^C]//g;
    $c=length($c)/$len;
    
    my $a=$all1{$pos};
    $a=~s/[^A]//g;
    $a=length($a)/$len;
    
    my $t=$all1{$pos};
    $t=~s/[^T]//g;
    $t=length($t)/$len;
    
    my $pos1=$pos+1;
    print OUT "$pos1\tA\t$a\n$pos1\tT\t$t\n$pos1\tC\t$c\n$pos1\tG\t$g\n";
}

#for my $pos( sort {$a<=>$b} keys %all2){
#    my $len=length($all2{$pos});
#    my $g=$all2{$pos};
#    $g=~s/[^G]//g;
#    $g=length($g)/$len;
#    
#    my $c=$all2{$pos};
#    $c=~s/[^C]//g;
#    $c=length($c)/$len;
#    
#    my $a=$all2{$pos};
#    $a=~s/[^A]//g;
#    $a=length($a)/$len;
#    
#    my $t=$all2{$pos};
#    $t=~s/[^T]//g;
#    $t=length($t)/$len;
#    
#    my $pos1=$pos+1;
#    print OUT "$pos1\tA\t$a\n$pos1\tT\t$t\n$pos1\tC\t$c\n$pos1\tG\t$g\n";
#}
#close OUT;