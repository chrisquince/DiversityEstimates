#!/usr/bin/perl

my $S = 0;
my $N = 0;
my $n1 = 0;
my $n2 = 0;

$line = <STDIN>;

while($line = <STDIN>){
    chomp($line);

    @tokens = split(/ /,$line);

    $abund = $tokens[0];
    $count = $tokens[1];

    $N += $abund*$count;
    $S += $count;

    if($abund == 1){
	$n1 = $count;
    }

    if($abund == 2){
	$n2 = $count;
    }
}


$chao = $S + 0.5*(($n1*$n1)/$n2);

print "$N $S $chao\n";
