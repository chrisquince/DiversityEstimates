#!/usr/bin/perl

@values = ();

$var   = shift(@ARGV);

$nBurn = shift(@ARGV);

for $file(@ARGV){

    open(FILE, $file) or die "Can't open $file\n";

    while($line = <FILE>){
	chomp($line);

	@tokens = split(/,/,$line);

	$time = $tokens[0];

	$value = $tokens[$var];
    
	if($time > $nBurn){
	    push(@values, $value);
	}
	
    }

    close(FILE);
}

printf("%.2f,", &mean(@values));
sub mean()
{
    my @array = @_;
    my $total = 0;
    
    my $size = scalar(@array);

    foreach my $val(@array){
	$total += $val;
    }

    return $total/$size;
}
