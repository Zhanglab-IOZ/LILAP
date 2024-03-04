my %all=();
open FILE,"sv.C01.ccs.hifiasm.txt";
while(<FILE>){
    chomp $_;
    my @line=split/\t/,$_;
    my @so=($line[1], $line[2]);
    @so=sort {$a<=>$b} @so;
    $all{$_}{chr}=$line[0];
    $all{$_}{start}=$so[0];
    $all{$_}{end}=$so[1];
}
close FILE;

open OUT,">sv.C01.filter1";
for my $sv1(sort keys %all){
    for my $sv2(sort keys %all){
        if ($sv1 ne $sv2 && $all{$sv1}{chr} ne 'no' && $all{$sv2}{chr} ne 'no'){
            if($all{$sv1}{chr} eq $all{$sv2}{chr}){
                unless ($all{$sv1}{start}>$all{$sv2}{end} || $all{$sv1}{end}<$all{$sv2}{start}){
                    print OUT "$sv1\t$sv2\n";
                    $all{$sv1}{chr}='no';
                    $all{$sv2}{chr}='no';
                }
            }
        }
    }
    
}
close OUT;
