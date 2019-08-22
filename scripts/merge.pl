#!/usr/bin/perl -w
# Author: Abdullah Kahraman
# Date: 04.03.2008

###############################################################################
###############################################################################
### Merges the two files.                                                   ###
###############################################################################
###############################################################################

use strict;
use Getopt::Long;

my (
    # variable for parameters which are read in from commandline
    $help,
    $infile,
    $in2,
    $del,
    $col1,
    $col2,
    $app,
    $remove,
    $uniq,
    $verbose,
    $skip,
    $ignoreCase,
   );

##############################################################################
### read all needed parameters from commandline ##############################

&GetOptions(
 "help!"        => \$help,        # print this help
 "in1=s"        => \$infile,      # path of file1, which will be used as hash.
 "in2:s"        => \$in2,         # path of file2.
 "del:s"        => \$del,         # delimiter (default is tab)
 "col1:s"       => \$col1,        # columns in common in both files (first column is 0)(several columns seperated by colon)
 "col2:s"       => \$col2,        # columns in common in both files (first column is 0)(several columns seperated by colon)
 "uniq!"        => \$uniq,        # ignore multiple occurences of values in file1. Uses only first occuring value.
 "append!"      => \$app,         # just append each line of file2 to each line of file1
 "remove:s"     => \$remove,      # removes this regular expression from columns prior to matching.
 "skip!"        => \$skip,        # does not output an empty line for entries that are not found in in1.
 "verbose!"     => \$verbose,     # show additional information during calculation
 "ignore!"      => \$ignoreCase,  # ignore case
) 
or die "\nTry \"$0 -h\" for a complete list of options\n\n";

##############################################################################
$| = 1;
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

    $usage = "$0 -in=<path> -field=[string] -stop=[integer] " .
             "-output [stdout|database] -help -list -verbose\n";


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

    $example  = "$0 -in=../../databases/locuslink/locuslink_031007.txt ";
}


##############################################################################
### END OF SUBROUTINES########################################################
##############################################################################

# help
if ($help) {printHelp(); exit}

##############################################################################

$del = "\t" unless(defined $del);

############
### MAIN ###
############
if (defined $infile and defined $in2 and defined $app) {
    my @fileContent;
    open(F2, "$in2") or die "Can't open file \"$in2\"\n";
    while (my $lineBuffer = <F2>){
	chomp($lineBuffer);
	push(@fileContent, $lineBuffer);
    }
    close(F2);
    open(F1, "$infile") or die "Can't open file \"$infile\"\n";
    my $i=0;
    while (my $lineBuffer = <F1>){
	chomp($lineBuffer);
	print $lineBuffer."\t".$fileContent[$i++]."\n";
    }
    close(F1);
    exit;
}

###########################################################################

if (defined $infile and defined $col1 and defined $col2) {
    $col1 = ":".$col1 if($col1 !~ /^:/);
    $col1 = $col1.":" if($col1 !~ /:$/);

    if($infile =~ /\.gz$/) {
	open(F1, "zcat $infile |") or die "\nERROR: Can't open gzipped reference file \"$infile\": $!\n\n";
    } else {
	open(F1, $infile) or die "\nERROR: Can't open txt reference file \"$infile\": $!\n\n";
    }

    print STDERR "Reading in reference file $infile ...\n" if($verbose);

    my %h1;
    my @cols1 = split(/:/,$col1);
    my $empty = "";
    my $i = 0;
    my %warningMsg = ();
    while(my $l = <F1>){
	$i++;
	if($l !~ /^#/){
	    chomp($l);
	    my @a = split(/$del/,$l);
	    my $key="";
	    for(my $i=1; $i<=$#cols1; $i++){
		$a[$cols1[$i]] =~ s/$remove// if(defined $remove);
		$key .= $a[$cols1[$i]]."\t";
	    }
	    chop($key);
	    my $value = "";
	    $empty = "";
	    for(my $i=0; $i<=$#a; $i++){
		if($col1!~/:$i:/){
		    $value .= $a[$i]."\t";
		    $empty .= "-\t";
		}
	    }
	    if(exists $h1{$key} and !defined $uniq) {
		my @a = @{$h1{$key}};
		push(@a, $value);
		$h1{$key} = \@a;
		$h1{lc($key)} = \@a if(defined $ignoreCase);
		$warningMsg{$key} = "WARNING: Multiple occurances of value $key in file $infile\n";
	    } else {
		my @a = ($value);
		$h1{$key} = \@a;
		$h1{lc($key)} = \@a if(defined $ignoreCase);
	    }
	}
    }
    close(F1);
    foreach my $key (sort keys %warningMsg) {
	print STDERR $warningMsg{$key} if(defined $verbose);
    }

    print STDERR "Reading in reference file $infile ... DONE\n" if($verbose);
    print STDERR "Reading in ".(defined $in2 ? "input file $in2" : "STDIN")."...\n" if($verbose);
    
    $col2 = ":".$col2 if($col2 !~ /^:/);
    $col2 = $col2.":" if($col2 !~ /:$/);

    if(defined $in2) {
	if($in2 =~ /\.gz$/) {
	    open(F2, "zcat $in2 |") or die "\nERROR: Can't open gzipped file \"$in2\"\n\n";
	} else {
	    open(F2, $in2) or die "\nERROR: Can't open txt file \"$in2\"\n\n";
	}
    }
    my @cols2 = split(/:/,$col2);
    while(my $l = (defined $in2 ? <F2> : <>)){
	if($l !~ /^#/){
	    chomp($l);
	    my @a = split(/$del/,$l);
	    my $key="";
	    for(my $i=1; $i<=$#cols2; $i++){
		$a[$cols2[$i]] =~ s/$remove// if(defined $remove);
		$key .= $a[$cols2[$i]]."\t";
	    }
	    chop($key);

	    my $value = "";
	    for(my $i=0; $i<=$#a; $i++){
		if($col2!~/:$i:/){
		    $value .= $a[$i]."\t";
		}
	    }
	    $value =~ s/\t$//;
	    if(exists $h1{$key}) {
		foreach my $hashValue (@{$h1{$key}}) {
		    print "$key\t$hashValue$value\n";
		}
	    } elsif(defined $ignoreCase and exists $h1{lc($key)}){
		foreach my $hashValue (@{$h1{lc($key)}}) {
		    print "$key\t$hashValue$value\n";
		}
	    } else{
		print STDERR "#WARNING: No entry to \"$key\" found in -in1\n" if($verbose);
		print "$key\t$empty$value\n" if(!defined $skip);
	    }
	}
    }
    close(F2);    
    print STDERR "Reading ".(defined $in2 ? "input file $in2" : "STDIN")." DONE\n" if($verbose);
}
else {
    die "\nTry \"$0 -h\" for a complete list of options\n\n";
}
