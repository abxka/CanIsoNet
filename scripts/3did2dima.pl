#!/usr/bin/perl -w
# Author: Abdullah Kahraman
# Date: 01.01.2000

###############################################################################
###############################################################################
### bla bla bla bla bla bla bla bla bla bla bla bla bla bla bla bla         ###
###############################################################################
###############################################################################

use strict;
use warnings;
use Getopt::Long;
no autovivification;

my (
    # variable for parameters which are read in from commandline
    $help,
   );

##############################################################################
### read all needed parameters from commandline ##############################

&GetOptions(
    "help!"      => \$help,     # print this help
) or die "\nTry \"$0 -h\" for a complete list of options\n\n";

##############################################################################

# help
if ($help) {printHelp(); exit}

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
### END OF SUBROUTINES########################################################
############################################################################## 


############
### MAIN ###
############

if(@ARGV < 1) {
    print STDERR "\n\tUsage: $0 3did_pfamInteractions_1689.tsv.gz\n\n";
    exit;
}

open(F1, "zcat $ARGV[0] |") or die "\n\tERROR: Failed to open 3did file $ARGV[0]: $!\n\n";
my %ipfam = ();
while(my $l = <F1>) {
    chomp($l);
    my ($pfam1, $pfam2) = split(/\t/, $l);
    $ipfam{$pfam1}->{$pfam2} = 1;
}
close(F1);

foreach my $pfam1 (sort keys %ipfam) {
    foreach my $pfam2 (sort keys %{$ipfam{$pfam1}}) {
	print "$pfam1\t$pfam2\t1\t\t1.0\t\t\t\t\t\n";
    }
}
