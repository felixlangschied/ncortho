#!/usr/bin/perl
use strict;
use Getopt::Long;

# PROGRAMNAME: degapper.pl

# AUTHOR: INGO EBERSBERGER, ingo.ebersberger@univie.ac.at

# PROGRAM DESCRIPTION:

# DATE: Thu May 29 11:11:16 CEST 2008


# DATE LAST MODIFIED:
######################## start main #############################
my $help;
my $infile = '';
my %final;
my $limit = 0.5;
my @sum = qw();
GetOptions ("h" => \$help,
	    "in=s" => \$infile,
	    "limit=s" => \$limit);
if ($help) {
	die "degapper.pl -in=<> [-limit=<>] [-h]

DESCRIPTION:
degapper.pl processes an alignment by removing all alignment columns where less than 'limit' species are represented by a valid amino acid.
OPTIONS:
-in: provided the filename of the alignment in fasta format
-limit: determines the minimum fraction of species that have to be represented by a amino acid in an alignment column. DEFAULT: 0.5";
}

my @seqs = `less $infile |sed -e 's/>//'`;
chomp @seqs;
my %smkeep = @seqs;
my %sm = @seqs;
my $speccount = scalar(keys %sm);
## translate all aminoacids into a 1, all other characters into a 0
## and sum up
for (keys %sm) {
    $sm{$_} =~ s/[^FLIMVSPTAYHQNKDECWRG]/0/ig;
    $sm{$_} =~ s/[A-Z]/1/gi;
    my @loc = split //, $sm{$_};
    for (my $i = 0; $i < @loc; $i++) {
	$sum[$i] += $loc[$i];
    }
}

## extract the alignment
for (my $i = 0; $i < @sum; $i++) {
    if ($sum[$i]/$speccount > $limit) {
	for (keys %sm) {
	    $final{$_} .= substr($smkeep{$_}, $i, 1);
	}
    }
}
open (OUT, ">$infile.proc") or die "could not open infile\n";
for (my $i = 0; $i < @seqs; $i+= 2) {
    print OUT ">$seqs[$i]\n$final{$seqs[$i]}\n";
}
close OUT;

