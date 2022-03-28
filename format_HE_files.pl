#!/usr/bin/perl -w                                                                                                                                                                     
use strict;

my $chr;
my $pos;
my $pos2;
my $strand;
my $region;
my $repeat;
my $gene;
my $distance;
my $motif;
my $line;

open(INFILE, $ARGV[0]) or die" The file is not the right one or missing or not specified\n";
while(<INFILE>)
{
    $line=$_;
    chomp $line;

if ($line=~m/\S+/)
{

($chr) = ($line =~m/(\S+)/);
($pos) = ($line =~m/\S+\s+(\S+)/);

if ($pos !~ m/Start/)
{
($strand) = ($line =~m/(\S+)$/);
($region) = ($line =~m/\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+(\S+)/);
($pos2) = ($line =~m/\S+\s+\S+\s+(\S+)/);
($gene) = ($line =~m/\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+(\S+)/);
($distance) = ($line =~m/\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+(\S+)/);
($repeat) = ($line =~m/(\S+)\s+\S+\s+\S+\s+\S+\s+\S+$/);
$repeat =~s/Name=//;
($motif) = ($line =~m/(\S+)\s+\S+$/);

if ($repeat !~ m/\./)

{
    my ($class) = ($repeat =~m/(\S+)\:\S+\:\S+/);
    my ($family) = ($repeat =~m/\S+\:(\S+)\:\S+/);
my ($repeat) = ($repeat =~m/\S+\:\S+\:(\S+)/);

#($motif) = ($line =~m/(\S+)\s+\S+$/);                                                                                                                                                 

    if ($distance =~m/NM/){
	print "$chr\t$pos\t$pos2\t$strand\t$region\t$class\t$family\t$repeat\t$gene\t\.\t$motif\n";
    }

    else {print "$chr\t$pos\t$pos2\t$strand\t$region\t$class\t$family\t$repeat\t$gene\t$distance\t$motif\n";}
}



elsif ($repeat =~ m/\./)
{

    if ($distance =~m/NM/){
        print "$chr\t$pos\t$pos2\t$strand\t$region\t\.\t\.\t\.\t$gene\t\.\t$motif\n";
    }

    else {print "$chr\t$pos\t$pos2\t$strand\t$region\t\.\t\.\t\.\t$gene\t$distance\t$motif\n";}
}


}
}
}


