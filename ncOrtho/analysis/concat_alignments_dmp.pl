#!/usr/bin/perl -w 
use strict;
use Getopt::Long;

# PROGRAMNAME: concat_alignments.pl

# AUTHOR: INGO EBERSBERGER, ingo.ebersberger@univie.ac.at

# PROGRAM DESCRIPTION:

# DATE: Thu Feb 22 18:51:40 CET 2007


# DATE LAST MODIFIED:
######################## start main #############################
my $help;
my $indir = '';
my $error;
my $maxlength = 0;
my @track;
my $refid;
my @files;
my $path = '.';
my $checklength = 0; # will be the sum of the individual alignment lengths
my $outfile;
my $concatenate; #hashref for the concatenated alignment
GetOptions ("h" => \$help,
	    "in=s" => \$indir,
	    "out=s"=> \$outfile);
#### usage and errors
my $usage = "USAGE: concat_alignments.pl -in=<> [-out=<>] [-h]
OPTIONS:
-in: specify the directory where the alignment files are located
-out: this gives you the option to specify the name of the output file
-h: Displays this help message
";

if ($help) {
    die "$usage\n";
}
if (!(defined($indir))) {
    die "Please provide the directory where the alignment files are located!\n$usage";
}
## the user has provided an indir
$indir =~ s/\/$//;
if ($indir =~ /\//) {
    ($path, $indir) = $indir =~ /(.*)\/(.*)/;
}
if (!(-e "$path/$indir")) {
    $error = "The specified directory $indir does not exist!\n";
}
else {
    opendir (INDIR, "$path/$indir") or die "could not open $indir\n";
    @files = grep !/^\./, readdir(INDIR);
    closedir(INDIR) or warn "could not close $indir. Not a fatal error!\n";
    if (@files == 0) {
	$error = "The specified directory is empty!\n";
    }
    else {
	print "Read in " . scalar(@files) . " potential alignment files that will be processed!\n";
    }
}
## something was wrong with the indir
if (defined $error) {
    die "$error\n\n$usage\n";
}
## the output
if (!(defined($outfile))) {
    $outfile = 'concat_alignments.fa';
}

## now do the concatenation
### read the file in
for (my $i = 0; $i < @files; $i++) {
    open(IN, "$path/$indir/$files[$i]") or die "could not open $indir/$files[$i]\n";
    my @content = <IN>;
    close (IN);
    chomp (@content);
    if (!(@content) or $content[0] !~ /^>/) {
	print "$files[$i] seems not to be in fasta-format. Skipped!\n";
	next;
    }
    else {
	&Concatenate($files[$i], @content);
    }
}
## generate the output
## print out the alignment
## check that sum of alignment lengths is the total length
open (OUT, ">$path/$outfile") or die "could not open outfile!\n";
for (keys %$concatenate) {
    if ($checklength = $concatenate->{$_}) {
	print OUT ">$_\n$concatenate->{$_}\n";
    }
    else {
	die "there is a problem. Generated alignment is longer than the sum of the individual alignments...\n";
    }
}
close OUT;
## print out the positional information
open (OUT, ">$path/subalignment_positions.txt") or die "could not open outfile\n";
print OUT join "\n", @track;
close OUT;


###########
sub Concatenate {
    my ($filename, @seq) = @_;
    my %localref;
    my $id;
    my $check; 
    if (defined($concatenate)) {
	$check = 1;
    }
    else {
	$check = 0;
    }
    ## concatenate the sequences in the file and store in hash
    ($id, %localref) = toHash(@seq);

    ## now walk through the hash.
    my $length = length($localref{$id}); 
    $checklength += $length;
    ## keep track of alignment start and end
    my $substart = $maxlength + 1;
    my $subend = $maxlength + $length;
    push @track, "$filename\t$substart\t$subend\t$length";
    ## now add the sequences to the concatenated alignment
    &addToConcat($check, $length, $filename, %localref);
 
    ## now make sure, that all sequences in the alignment have the same length. This guarantees
    ## that taxa present in the hashref, but not present in a subsequent processed alignment get
    ## represented by a string of 'X'
    for (keys %$concatenate) {
	if (length($concatenate->{$_}) < $maxlength) {
	    my $diff = $maxlength - length($concatenate->{$_});
	    $concatenate->{$_} .= 'X' x $diff;
	}
    }
}
##########
sub toHash {
    my @seq = @_;
    my %localref;
    my $id;
    for (my $i = 0; $i < @seq; $i++) {
	if ($seq[$i] =~ /^>/) {
	    if ($seq[$i] =~ /\s/) {
		($id) = $seq[$i] =~ />(.*?)\s.*/;
	    }
	    else {
		($id) = $seq[$i] =~ />(.*)$/;
	    }
	    if (!(defined $refid)) {
		$refid = $id;
	    }
	    $localref{$id} = '';
	}
	else {
	    chomp $seq[$i];
	    $localref{$id} .= $seq[$i];
	}
    }
    return ($id,%localref);
}
#########
sub addToConcat {
    my ($check, $length, $filename, %localref) = @_;
    my $oldlength;
    ## determine the length of the alignment up to now. This will be
    ## needed for adding taxa that were not represented in the so far
    ## concatenated alignments.
    if (defined $concatenate->{$refid}) {
	$oldlength = length($concatenate->{$refid});
    }
    else {
	$oldlength = 0;
    }
    ## now add the alignment sequences to the concatenated alignment
    ## we have to distinguish different cases
    for (keys %localref) {
	## first: check if the alignment lengths are the same
	if (length($localref{$_}) != $length) {
	    die "problem in $filename. Sequences have different lenghts\n";
	}
	## second: test if the hashref already contains any sequences. Note
	## that the test-variable gets its value when the subroutine is called.
	if ($check == 1) {
	    ## first test passed, now test if the taxon is represented
	    ## case: NO
	    if (!(defined($concatenate->{$_}))) {
		$concatenate->{$_} = 'X' x $oldlength;
		$concatenate->{$_} .= $localref{$_};
		$maxlength = length($concatenate->{$_});
	    }
	    ## case: YES
	    else {
		## if first test passed, test if the taxon is already represented
		## in the hashref representing the concatenated alignment. If so,
		## add the sequence
		$concatenate->{$_} .= $localref{$_};
		$maxlength = length($concatenate->{$_});
	    }
	}
	## the hashref contains no sequences, yet
	else {
	    $concatenate->{$_} = $localref{$_};
	    $maxlength = length($concatenate->{$_});
	}
    }
}
