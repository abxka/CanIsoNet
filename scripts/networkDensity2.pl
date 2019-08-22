#!/usr/bin/perl -w
# Author: Abdullah Kahraman
# Date: 26.10.2016

###############################################################################
###############################################################################
### Calculates for every node in STRING an interaction density score.       ###
###############################################################################
###############################################################################

use strict;
use warnings;
use Getopt::Long;
use List::Util qw[min max];
$|=1;

my (
    # variable for parameters which are read in from commandline
    $help,
    $stringFile,
    $minStringScore,
    $maxNeighbourhood,
    $verbose,

    %string,
    %int,
    %nInt,
    $count,
   );

my $interactionSign = "<=>";
##############################################################################
### read all needed parameters from commandline ##############################

&GetOptions(
    "help!"    => \$help,               # print this help.
    "in=s"     => \$stringFile,         # STRING DB e.g. protein.links.v10.human.txt.gz
    "min:i"    => \$minStringScore,     # minimum confidence score of STRING interaction, default is 400.
    "max:i"    => \$maxNeighbourhood,   # maximum neighbourhood size to search for the interaction partner, default is 1.
    "verbose!" => \$verbose,            # output additional information on computational progess.
) or die "\nTry \"$0 -h\" for a complete list of options\n\n";

##############################################################################

# help
if ($help) {printHelp(); exit}

$minStringScore ||= 400;
$maxNeighbourhood ||= 1;
$count = 0;
##############################################################################
### SETTINGS #################################################################
##############################################################################

##############################################################################
### SUBROUTINES ##############################################################
############################################################################## 

###############################################################################
sub printHelp {
###############################################################################
    # prints a help about the using and parameters of this scripts 
    # (execute if user types commandline parameter -h)
    # param:  no paramaters
    # return: no return value

    my (
	$usage,
	$sourceCode,
	@rows,
	$row,
	$option,
	$scriptInfo,
	$example,
       );

    $usage = "$0\n";


    print "\nUsage: " .  $usage . "\n";

    print "Valid options are:\n\n";
    open(MYSELF, "$0") or
      die "Cannot read source code file $0: $!\n";
    $sourceCode .= join "", <MYSELF>;
    close MYSELF;
    $sourceCode =~ s/^.+?\&GetOptions\(\n//s;
    $sourceCode =~ s/\n\).+$//s;
    @rows = split /\n/, $sourceCode;
    foreach $row (@rows){
        $option = $row;
	$option =~ s/\s+\"//g;
	$option =~ s/\"\s.+\#/\t\#/g;
	$option =~ s/=./\t<value> [required]/;
	$option =~ s/:./\t<value> [optional]/;
	$option =~ s/!/\t<non value> [optional]/;

	$row =~ s/^.*//;
	print "\t";
	printf("%-1s%-30s%-30s\n", "-",$option,$row);
    } # end of foreach $row (@rows)
    print "\n";
    print "Options may be abreviated, e.g. -h for --help\n\n";

    $example  = "$0";
}
############################################################################## 
sub readStringFile {
    open(F4, "zcat $stringFile|") or die "ERROR: Failed to open $stringFile\n";
    while(<F4>) {
	next if(! /^\d+/);
	chomp($_);
	my @a=split(/\s+/);
#	next if($a[0]!~/^$taxId\./ or $a[1]!~/^$taxId\./);
	$a[0] =~ s/^\d+\.//;
	$a[1] =~ s/^\d+\.//;
	$a[0] =~ s/(ENSP\d+)\..*/$1/;
	$a[1] =~ s/(ENSP\d+)\..*/$1/;
	if($a[-1] >= $minStringScore) {
	    $int{$a[0]}++ if(!exists $string{$a[0]}->{$a[1]}); 
	    $int{$a[1]}++ if(!exists $string{$a[1]}->{$a[0]}); 
	    $string{$a[0]}->{$a[1]} = $a[-1]/1000;
	    $string{$a[1]}->{$a[0]} = $a[-1]/1000;
	}
    }
    close(F4);
}
##############################################################################
sub getInteractors {
    (my $exclude, my $p1) = @_;

    my $interactors = $string{$p1};
    my @interactors;
    foreach my $p2 (sort {$interactors->{$b}<=>$interactors->{$a}} keys %$interactors) {
	next if($exclude =~ /\b$p2\b/);
	push(@interactors, $p2);
    }
    return \@interactors;
}

##############################################################################
sub breadth1stSearch {
    my ($proteins, $root) = @_;

    $count++;

    my @proteins = @$proteins;
    my @neighbours;
    my @neighboursScore;

    print STDERR "\tChecking $count-th degree neighbourhood. Need to check ".($#proteins+1)." proteins ... " if($verbose);
    for(my $i=0; $i<=$#proteins; $i++) {
	my  $p1 = $proteins[$i];
	my $lastProtein = $p1;
	$lastProtein =~ s/.*$interactionSign//;
	my @interactors = @{&getInteractors($p1, $lastProtein)};
	if($#interactors >= 0){
	    for(my $j=0; $j<=$#interactors; $j++){
		my $p2 = $interactors[$j];
		# use negative exponential as a weighting factor
		$nInt{$root} += (1/(2**$count)) * $int{$p2} * $string{$lastProtein}->{$p2};
		$interactors[$j] = "$p1$interactionSign$p2";
	    }
	    push(@neighbours, @interactors);
	}
    }
    print STDERR "Done\n" if($verbose);

    if($count == $maxNeighbourhood) {
	return;
    }
    &breadth1stSearch(\@neighbours, $root);
}
##############################################################################
sub searchNeighbourhood {
    (my $prot) = @_;
    my @proteins = ($prot);
    $count = 0;
    %nInt = ();
    &breadth1stSearch(\@proteins, $prot);
}
############################################################################## 
### END OF SUBROUTINES########################################################
############################################################################## 


############
### MAIN ###
############
if(!defined $stringFile ) {
    print "\nPlease provide the path to an STRING interaction file. Try \"$0 -h\" for a complete list of options.\n\n";
} else {
    print STDERR "Reading in STRING interactions from $stringFile ... " if(defined $verbose);
    &readStringFile();
    print STDERR "Done\n" if(defined $verbose);

    my %density = ();
    print "#Protein1\tProtein2\tDensityScore1\tDensityScore2\tInteractionScore\n";
    foreach my $ensp1 (sort keys %string) {
	if(!exists $density{$ensp1}) {
	    &searchNeighbourhood($ensp1);
	    $density{$ensp1} = $int{$ensp1} + $nInt{$ensp1};
	}
	foreach my $ensp2 (sort keys %{$string{$ensp1}}) {	    
	    if(!exists $density{$ensp2}) {
		&searchNeighbourhood($ensp2);
		$density{$ensp2} = $int{$ensp2} + $nInt{$ensp2};
	    }
	    printf("%s\t%s\t%.4f\t%.4f\t%.4f\t%i\t%i\n",
		   $ensp1, $ensp2, $density{$ensp1}, $density{$ensp2}, $string{$ensp1}->{$ensp2}, $int{$ensp1}, $int{$ensp2});
	}
    }
}
