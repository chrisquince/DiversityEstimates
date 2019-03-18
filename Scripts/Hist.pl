#!/usr/bin/perl

if($ARGV[0] eq "h"){
    $line = <STDIN>;
}

@abundance = ();

my $bin = 0;
my $binStart = 1;
my @binCounts = ();
my @start     = ();
my @end       = ();

$start[0] = 1;

while($line = <STDIN>){
    chomp($line);
    
    @tokens = split(/\s+/,$line);
    
    $abund = $tokens[0];

    $binCounts[$bin] += $tokens[1];
    $counts[$bin]++;
    $sum[$bin] += $abund*$tokens[1];

    if($binCounts[$bin] > 20){
	$end[$bin] = $abund;
	$bin++;
	$start[$bin] = $abund + 1;
    }
}

$end[$bin] = $abund;

#print "@binCounts\n";
#print "@sum\n";
#print "@counts\n";

$nBins = scalar(@binCounts);

for($i = 0; $i < $nBins; $i++){
    $density = $binCounts[$i]/($end[$i] - $start[$i] + 1);
    $mean    = $sum[$i]/$binCounts[$i];
    printf("%f %f\n",$mean, $density);
}
