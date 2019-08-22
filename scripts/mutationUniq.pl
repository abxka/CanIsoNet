#!/usr/bin/perl -w
# Author: Abdullah Kahraman
# Date: 26.06.2017

###############################################################################
###############################################################################
### Counts mutations shared between regions only once. Prioritization is as ###
### follows: cds->ss->intron->promCore->promoter->5'UTR->promDomain->3'UTR  ###
### ->enhancer.                                                             ###           
###############################################################################
###############################################################################

use strict;
use warnings;
use Getopt::Long;
no autovivification;

my (
    # variable for parameters which are read in from commandline
    $help,
    $mutCol,
    $eqtlCol,
    $promDomain,
    $formatV2,
   );

##############################################################################
### read all needed parameters from commandline ##############################

&GetOptions(
    "help!"       => \$help,       # print this help
    "mutCol=i"    => \$mutCol,     # e.g. 50. Starts from 0.
    "eqtlCol=i"   => \$eqtlCol,    # e.g. 51. Starts from 0.
    "promDomain!" => \$promDomain, # includes promoter domain mutations
    "v2!"         => \$formatV2,   # format v2.
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

if(!defined $mutCol or !defined $eqtlCol) {
    print STDERR "\n\tPlease provide a column number for mutational and eQTL information\n\n";
    exit;
}

while(my $l = <>) {
    if($l =~ /^#/) {
	print $l;
	next;
    }
    chomp($l);
    my @a = split(/\t/, $l);    
    print $a[0];
    for(my $i=1; $i < @a-2; $i++) {
	print "\t$a[$i]";
    }

    if(defined $mutCol) {
	my @m = split(/,/, $a[$mutCol]);
	my %uniqMut = ();
	my %m = ();
	foreach my $m (@m){
	    my ($regionType, $chr, $rStartEnd, $varType, $mStartEnd, $tumorSampleBarcode) = split(/\:/, $m);
	    ($chr, my $strand, $rStartEnd, $regionType, $mStartEnd, $varType) = split(/\:/, $m) if(defined $formatV2);
	    next if($m eq "-");
	    next if($regionType eq "gc19_pc.promDomain" and !defined $promDomain);
	    if(exists $uniqMut{"$chr:$mStartEnd"}) {
		if($uniqMut{"$chr:$mStartEnd"} eq "gc19_pc.cds") {
		    next;
		} elsif($regionType eq "gc19_pc.cds") {
		    $uniqMut{"$chr:$mStartEnd"} = $regionType;
		    $m{"$chr:$mStartEnd"} = $m;
		    next;
		} elsif($uniqMut{"$chr:$mStartEnd"} eq "gc19_pc.ss") {
		    next;
		} elsif($regionType eq "gc19_pc.ss") {
		    $uniqMut{"$chr:$mStartEnd"} = $regionType;
		    $m{"$chr:$mStartEnd"} = $m;
		    next;
		} elsif($uniqMut{"$chr:$mStartEnd"} eq "gc19_pc.promCore") {
		    next;
		} elsif($regionType eq "gc19_pc.promCore") {
		    $uniqMut{"$chr:$mStartEnd"} = $regionType;
		    $m{"$chr:$mStartEnd"} = $m;
		    next;
		} elsif($uniqMut{"$chr:$mStartEnd"} eq "gc19_pc.5utr") {
		    next;
		} elsif($regionType eq "gc19_pc.5utr") {
		    $uniqMut{"$chr:$mStartEnd"} = $regionType;
		    $m{"$chr:$mStartEnd"} = $m;
		    next;
		} elsif($uniqMut{"$chr:$mStartEnd"} eq "gc19_pc.3utr") {
		    next;
		} elsif($regionType eq "gc19_pc.3utr") {
		    $m{"$chr:$mStartEnd"} = $m;
		    $uniqMut{"$chr:$mStartEnd"} = $regionType;
		    next;
		} elsif($uniqMut{"$chr:$mStartEnd"} eq "gc19_pc.promDomain") {
		    next;
		} elsif($regionType eq "gc19_pc.promDomain") {
		    $m{"$chr:$mStartEnd"} = $m;
		    $uniqMut{"$chr:$mStartEnd"} = $regionType;
		    next;
		} else {
		    $m{"$chr:$mStartEnd"} = $m;
		    $uniqMut{"$chr:$mStartEnd"} = $regionType;
		}
	    } else {
		$m{"$chr:$mStartEnd"} = $m;
		$uniqMut{"$chr:$mStartEnd"} = $regionType;
	    }
	}
    
	my $t = "";
	foreach my $m (sort keys %m) {
	    $t .= "$m{$m},";
	}
	$t = "-" if($t eq ""); 
	$t =~ s/\,$//;
	print "\t$t";
    } 
    if(defined $eqtlCol) {
	my @e = split(/,/, $a[$eqtlCol]);

	my %uniqEqtl = ();
	my %e = ();
	foreach my $e (@e){
	    my ($eelement, $eType, $chr, $eStartEnd, $gpval, $beta, $mStartEnd, $tumorSampleBarcode) = split(/:/, $e);
	    ($chr, $eStartEnd, $eelement, $eType, $gpval, $beta, $mStartEnd, my $varType) = split(/\:/, $e) if(defined $formatV2);
	    next if($e eq "-");
	    if(exists $uniqEqtl{"$chr:$mStartEnd"}){
		if($uniqEqtl{"$chr:$mStartEnd"} eq "intron") {
		    next;
		} elsif($eelement eq "intron") {
		    $uniqEqtl{"$chr:$mStartEnd"} = $eelement;
		    $e{"$chr:$mStartEnd"} = $e;
		} elsif($uniqEqtl{"$chr:$mStartEnd"} eq "promoter") {
		    next;
		} elsif($eelement eq "promoter") {
		    $uniqEqtl{"$chr:$mStartEnd"} = $eelement;
		    $e{"$chr:$mStartEnd"} = $e;
		} elsif($uniqEqtl{"$chr:$mStartEnd"} eq "enhancer") {
		    next;
		} elsif($eelement eq "enhancer") {
		    $uniqEqtl{"$chr:$mStartEnd"} = $eelement;
		    $e{"$chr:$mStartEnd"} = $e;
		} else { 
		    $e{"$chr:$mStartEnd"} = $e;
		    $uniqEqtl{"$chr:$mStartEnd"} = $eelement;
		}
	    } else {
		$e{"$chr:$mStartEnd"} = $e;
		$uniqEqtl{"$chr:$mStartEnd"} = $eelement;
	    }
	}
	
	my $t = "";
	foreach my $e (sort keys %e) {
	    $t .= "$e{$e},";
	}
	$t = "-" if($t eq ""); 
	print "\t$t";
    }
    print "\n";
}

