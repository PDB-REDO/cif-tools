#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;

my $input = shift;
die "Kan bestand niet vinden\n" unless -f $input;

my @isomerSets;
my %isomerIdx;

open(my $h, "<$input") or die "Kan bestand niet openen: $!\n";
while (my $line = <$h>)
{
	chomp($line);
	my ($a, $b) = split(m/ is isomer of /, $line);
	
	next if $a eq 'UNK' or $b eq 'UNK';
	
	my $ai = $isomerIdx{$a};
	my $bi = $isomerIdx{$b};
	
	if (defined $ai and defined $bi)
	{
		die "$ai != $bi voor $a en $b\n" unless $ai == $bi;
	}
	elsif (defined $ai)
	{
		push @{$isomerSets[$ai]}, $b;
		$isomerIdx{$b} = $ai;
	}
	elsif (defined $bi)
	{
		push @{$isomerSets[$bi]}, $a;
		$isomerIdx{$a} = $bi;
	}
	else
	{
		my $ix = $#isomerSets + 1;

		my @s = ( $a, $b );
		
		push @isomerSets, \@s;

		$isomerIdx{$a} = $ix;
		$isomerIdx{$b} = $ix;
	}	
}
close($h);

# print Dumper(\%isomerIdx, \@isomerSets);

print<<EOF;
<?xml version="1.0"?>
<isomers>
EOF

foreach my $set (@isomerSets)
{
	print "  <isomer-set>\n";
	foreach my $s (@{$set})
	{
		print "    <compound>$s</compound>\n";	
	}
	print "  </isomer-set>\n";
}
print "</isomers>\n";
