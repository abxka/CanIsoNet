#!/usr/bin/perl -w
# Author: Abdullah Kahraman
# Date: 10.11.2016
# version: 0.1
###############################################################################
###############################################################################
### Adds various annotations to a missingInteractionsInDTU_filt_int.tsv.gz  ###
### file.                                                                   ###
###############################################################################
###############################################################################

use strict;
use warnings;
use Getopt::Long;
no autovivification;

my (
    # variable for parameters which are read in from commandline
    $help,
    $ensgFile,
    $stringDensFile,
    $cosmicFile,
    $keggCancerFile,
    $step,
    $keggFile,
    $expressFile,
    $isoIntFile, # no need anymore. Random should be selected from STRING interactions
    $ensp1col,
    $ensp2col,
    $ensg1col,
    $minExpress,
    $verbose,
   );

##############################################################################
### read all needed parameters from commandline ##############################

&GetOptions(
    "help!"             => \$help,           # print this help
    "ensgFile=s"        => \$ensgFile,       # e.g. ensg_ensp_enst_ense_geneName_v75.tsv.gz
    "densityFile=s"     => \$stringDensFile, # e.g. string_v10_density_900_3.tsv.gz
    "cosmicFile=s"      => \$cosmicFile,     # e.g. cosmic_census_20160527_network.tsv.gz
    "keggCancerFile=s"  => \$keggCancerFile, # e.g. keggCancerPathways_ensg.tsv.gz
    "keggFile=s"        => \$keggFile,       # e.g. kegg_benchmarking.CONN_maps_in_human.tsv
    "expressFile=s"     => \$expressFile,    # e.g. isoformExpressionInCancer.tsv.gz
    "isoIntFile:s"      => \$isoIntFile,     # e.g. interactionsInIsoforms_0_2.tsv.gz
    "step:s"            => \$step,           # restrict calculations on steps. Default 1234
    "ensp1col:i"        => \$ensp1col,       # column number of 1st STRING ENSP. Starts with 0.
    "ensp2col=i"        => \$ensp2col,       # column number of 2nd STRING ENSP. Starts with 0.
    "ensg1col:i"        => \$ensg1col,       # column number of 1st ENSG. Starts with 0.
    "minExpress:s"      => \$minExpress,     # Minimum expression for random background determination. Default is 2.
    "verbose!"          => \$verbose,        # print out additional information.
) or die "\nTry \"$0 -h\" for a complete list of options\n\n";

##############################################################################

# help
if ($help) {printHelp(); exit}

##############################################################################
### SETTINGS #################################################################
##############################################################################
$step = "1234" if(!defined $step);
my $randMaxRep = 5000;
my $expressWindow = 0.1; # +-10%;
$minExpress = 2 if(!defined $minExpress);
srand(23);
$| = 1;
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
sub readEnsgFile {
    my %p = ();
    my %g = ();
    open(F1,"zcat $ensgFile|") or die "ERROR: Failed to open ENSG file $ensgFile: $!\n";
    my %d = ();
    while(<F1>){
	next if(!/^ENSG/);
	chomp($_);
	my ($ensg, $ensp,$enst,$ense, $gene, $db) = split(/\t/);
	next if($ensp !~ /ENSP/);

	$p{$ensp} = $ensg;
	if(exists $g{$ensg}) {
	    if($d{$ensg} =~ /^HGNC/){
		$g{$ensg} = $gene;
		$g{$ensp} = $gene;
		$d{$ensg} = $db;
	    }
	} else {
	    $g{$ensg} = $gene;
	    $g{$ensp} = $gene;
	    $d{$ensg} = $db;
	}
    }
    close(F1);
    die "$ensgFile has wrong format. Should have in 1st column ENSG, 2nd column ENSP IDs\n"
	if(keys %p == 0);
    return (\%p, \%g);
}
##############################################################################
sub readStringDensityFile {
    my %d = ();
    my %n = ();
    my %r = ();
    open(F2,"zcat $stringDensFile|") or die "ERROR: Failed to open STRING density file $stringDensFile: $!\n";
    my $i = 0;
    while(<F2>){
	chomp($_);
	next if(/^#/);
	my @a = split(/\t/);
	$d{$a[0]} = $a[2];
	$d{$a[1]} = $a[3];
	$n{$a[0]} = $a[5];
	$n{$a[1]} = $a[6];
    }
    close(F2);
  
    my $n = keys %d;  
    # calculate rank based score
    my $j = 0;
    foreach my $ensp (sort {$d{$a}<=>$d{$b}} keys %d) {
	$j++;
	$r{$ensp} = $j/$n;
    }

    return (\%n, \%d, \%r);
}
##############################################################################
sub readCosmicFile {
    my %h = ();
    my %n = ();
    open(F3, "zcat $cosmicFile|") or die "ERROR: Failed to open COSMIC file $cosmicFile: $!\n";
    while(<F3>){
	next if(/^#/);
	chomp($_);
	my @a = split(/\t/);
	if (!exists $h{$a[2]}->{$a[0]}){
	    $h{$a[2]}->{$a[0]} = $a[5];
	    $n{$a[2]}->{$a[0]} = $a[1];
	} else{
	    if($a[5] > $h{$a[2]}->{$a[0]}) {
		$h{$a[2]}->{$a[0]} = $a[5];
		$n{$a[2]}->{$a[0]} = $a[1];
	    }
	}
    }
    close(F3);
    return(\%h, \%n);
}
##############################################################################
sub readKeggCancerFile {
    my %h = ();
    open(F4,"zcat $keggCancerFile|") or die "ERROR: Failed to open KEGG cancer pathway file $keggCancerFile: $!\n";
    while(<F4>){
	chomp($_);
	next if(/^#/);
	my @a = split(/\t/);
	if(exists $h{$a[5]}) {
	    $h{$a[5]} .= ",$a[0]" if($h{$a[5]} !~ /\b$a[0]\b/); 
	} else {
	    $h{$a[5]} = "$a[0]";
	}
	if(!exists $h{$a[6]}) {
	    $h{$a[6]} = "$a[0]";
	} else {
	    $h{$a[6]} .= ",$a[0]" if($h{$a[6]} !~ /\b$a[0]\b/); 
	}
    }
    close(F4);
    return \%h;
}
##############################################################################
sub readKeggFile {
    my %h = ();
    open(F5, $keggFile) or die "ERROR: Failed to open KEGG pathway file $keggFile: $!\n";
    while(<F5>){
	chomp($_);
	next if(/^#/);
	my @a = split(/\t/);
	my @ensp = split(/\s/, $a[3]);
	foreach my $ensp (@ensp) {
	    if(!exists $h{$ensp}) {
		$h{$ensp} = $a[1];
	    } else {
		$h{$ensp} .= ",$a[1]" if($h{$ensp}!~/$a[1]/);
	    }
	}
    }
    close(F5);
    return \%h;
}
###############################################################################
sub readIsoExpressFile {
###############################################################################
    my %express = ();
    my %stringExpress = ();
    open(F6, "zcat $expressFile|") or die "\nERROR: Failed to open $expressFile: $!\n\n";
    while(my $l = <F6>) {
	chomp($l);
	my ($cancerType, $aliquotId, $ensg, $stringEnsp, $ensp, $enst,
	    $totalGeneExpression, $isoformExpression, $ratio) = split(/\s/, $l);

	$express{$ensp}->{$ensg} = $totalGeneExpression;
	$stringExpress{$ensg}->{$ensp} = $stringEnsp;
    }
    close(F6);
    return (\%express, \%stringExpress);
}
###############################################################################
sub readIsoIntFile {
###############################################################################
    my %isoInt = ();
    open(F7, "zcat $isoIntFile|") or die "\nERROR: Failed to open $isoIntFile: $!\n\n";
    while(my $l = <F7>) {
	next if($l =~ /^#/);
	chomp($l);
	my ($ensg, $ensp, $enst,
	    $stringEnsp, $stringIntN, $existIntN, $missIntN, $relMissIntN,
	    $existInts, $missInts) = split(/\t/, $l);
	$isoInt{$stringEnsp} = 1 if($existIntN > 0 or $missIntN > 0);
    }
    close(F7);
    return \%isoInt;
}
##############################################################################
### END OF SUBROUTINES########################################################
############################################################################## 


############
### MAIN ###
############

my $error = 0;
if(!defined $ensgFile) {
    print STDERR "Please provide an ENSG file.\n";
    $error = 1;
}
if(!defined $stringDensFile and $step =~ /1/) {
    print STDERR "Please provide a STRING density file.\n";
    $error = 1;
}
if(!defined $cosmicFile and $step =~ /2/) {
    print STDERR "Please provide a Cosmic file.\n";
    $error = 1;
}
if(!defined $keggCancerFile and $step =~ /3/) {
    print STDERR "Please provide a KEGG cancer file.\n";
    $error = 1;
}
if(!defined $keggFile and $step =~ /4/) {
    print STDERR "Please provide a KEGG file.\n";
    $error = 1;
}
if(!defined $expressFile) {
    print STDERR "Please provide an isoform expression file file.\n";
    $error = 1;
}
if(!defined $ensp1col) {
    print STDERR "Please provide the column number of the 1st STRING ENSP ID.\n";
    $error = 1;
}
if(!defined $ensp2col) {
    print STDERR "Please provide the column number of the 2nd STRING ENSP ID.\n";
    $error = 1;
}
if(!defined $ensg1col) {
    print STDERR "Please provide the column number of the 1st ENSG ID.\n";
    $error = 1;
}
if($error) {
    print "\n\tTry $0 -help to get more information.\n\n";
    exit 0;
}

my ($string, $int, $intDens, $intRank);
my ($cosmicMinScore, $cosmicIntProt);
my $keggCancer;
my $kegg;

print STDERR "Reading in $expressFile ..." if($verbose);
my ($express, $stringExpress) = &readIsoExpressFile();
print STDERR "done\n" if($verbose);

print STDERR "Reading in $ensgFile ..." if($verbose);
my ($ensp2ensg, $ensg2gene) = &readEnsgFile();
print STDERR "done\n" if($verbose);

my $isoInt = ();
if(defined $isoIntFile) {
    print STDERR "Reading in $isoIntFile ..." if($verbose);
    $isoInt = &readIsoIntFile();
    print STDERR "done\n" if($verbose);
}

if($step =~ /1/) {
    print STDERR "Reading in $stringDensFile ..." if($verbose);
   ($int, $intDens, $intRank) = &readStringDensityFile();
    print STDERR "done\n" if($verbose);
}
if($step =~ /2/) {
    print STDERR "Reading in $cosmicFile ..." if($verbose);
    ($cosmicMinScore, $cosmicIntProt) = &readCosmicFile();
    print STDERR "done\n" if($verbose);
}
if($step =~ /3/) {
    print STDERR "Reading in $keggCancerFile ..." if($verbose);
    $keggCancer = &readKeggCancerFile();
    print STDERR "done\n" if($verbose);
}
if($step =~ /4/) {
    print STDERR "Reading in $keggFile ..." if($verbose);
    $kegg = &readKeggFile();
    print STDERR "done\n" if($verbose);
}

# for random selection, choose a transcript that is expressed and has multiple isoforms
my @ensp = ();
foreach my $ensg (sort keys %$stringExpress) {
    if(keys %{$stringExpress->{$ensg}} > 1) {
	foreach my $ensp (sort keys %{$stringExpress->{$ensg}}) {
	    if($express->{$ensp}->{$ensg} >= $minExpress) {
		my $stringEnsp = $stringExpress->{$ensg}->{$ensp};
		next if(defined $isoIntFile and !exists $isoInt->{$stringEnsp});
		push(@ensp, $stringEnsp);
	    }
	}
    }
}
my %random = ();
my %selRandom = ();

print STDERR "Starting mapping ...\n" if($verbose);
my $nline = 0;
while(<>) {
    print STDERR "\t$nline entries mapped ...\n" if($verbose and $nline++ % 100 == 0);
    chomp($_);
    if(/^#/) {
	print $_;
	print "\tENSG2\tGeneName1\tGeneName2";
	print "\tSTRINGintN1\tSTRINGdensityScore1\tSTRINGdensityRank1\t".
	      "STRINGintN2\tSTRINGdensityScore2\tSTRINGdensityRank2" if($step =~ /1/);
	print "\tENSPcanonCosmic\tENSPcanonCosmicDist\t".
 	        "ENSPcanonCosmicMinScore\tENSPcanon2cosmic\t".
		"ENSPcanon2cosmicDist\tENSPcanon2cosmicMinScore\t".
		"ENSPrandom\tENSPrandomCosmic\tENSPrandomCosmicDist\tENSPrandomCosmicMinScore\t".
		"MeanRandomDist\tMeanRandomMinScore" if($step =~ /2/);
	print "\tKeggCancerPathway1\tKeggCancerPathway2" if($step =~ /3/);
	print "\tKeggPathway1\tKeggPathway2" if($step =~ /4/);
	print "\n";
	next;
    }
    my @a = split(/\t/);
    print "$_";
    my $ensg1 = $a[$ensg1col];
    my $ensp1 = $a[$ensp1col];
    my $ensp2 = $a[$ensp2col];
    my $ensg2 = "-";
    if(exists $ensp2ensg->{$ensp2}) {
	$ensg2 = $ensp2ensg->{$ensp2};
	print "\t$ensg2"; 
	if(exists $ensg2gene->{$ensg1}) {
	    print "\t$ensg2gene->{$ensg1}\t$ensg2gene->{$ensg2}";
	} else {
	    print "\t-\t$ensg2gene->{$ensg2}";
	}
    } else {
	print "\t-";
	print STDERR "STRING ID $ensp2 unknown in ENSEMBL. Assigning \"-\" !\n" if($ensp2 ne "-");
	if(exists $ensg2gene->{$ensg1}) {
	    print "\t$ensg2gene->{$ensg1}\t-";
	} else {
	    print "\t-\t-";
	}
    }
    if($step =~ /1/) {
	if(exists $int->{$ensp1}) {
	    print "\t$int->{$ensp1}\t$intDens->{$ensp1}\t$intRank->{$ensp1}";
	} else {
	    print "\t0\t0\t0";
	    print STDERR "No interactions for $ensp1. Assigning \"0\t0\t0\" !\n" if($ensp1 ne "-");
	}
	if(exists $int->{$ensp2}) {
	    print "\t$int->{$ensp2}\t$intDens->{$ensp2}\t$intRank->{$ensp2}";
	} else {
	    print "\t0\t0\t0";
	    print STDERR "No interactions for $ensp2. Assigning \"0\t0\t0\" !\n" if($ensp2 ne "-");
	}
    }
    if($step =~ /2/) {
	if(exists $cosmicMinScore->{$ensp1}){
	    foreach my $dist (sort {$a<=>$b} keys %{$cosmicMinScore->{$ensp1}}) {
		print "\t$cosmicIntProt->{$ensp1}->{$dist}\t$dist\t$cosmicMinScore->{$ensp1}->{$dist}";
		last;
	    }
	} else {
	    print "\t-\t4\t0";
	}
	if(exists $cosmicMinScore->{$ensp2}){
	    foreach my $dist (sort {$a<=>$b} keys %{$cosmicMinScore->{$ensp2}}) {
		print "\t$cosmicIntProt->{$ensp2}->{$dist}\t$dist\t$cosmicMinScore->{$ensp2}->{$dist}";
		last;
	    }
	} else {
	    print "\t-\t4\t0";
	}
	
	# select random protein that has expression similar to DTU gene expression +- 10%. Try maximum 1000 times.
	my $randomEnsp = "";

	# if you have already determined a random protein for ENSP1, reuse it.
	if(!exists $random{$ensp1}) {
	    my $n = 0;
	    while(1) {
		$randomEnsp = $ensp[int(rand($#ensp + 1))];

		$n++;
		if($n >= $randMaxRep) {
		    print STDERR "WARNING: $n attempts done. Reached maximum random repetition cylcle $randMaxRep. Using $randomEnsp as random protein. Changing ID to _"."$randomEnsp\n";
		    $randomEnsp = "_".$randomEnsp;
		    last;
		}
		
		next if(exists $selRandom{$randomEnsp}); # randomEnsp should be selected only once, as there is also only one ensp in the dataset.
		last;
	    }
	    $random{$ensp1} = $randomEnsp;
	    my $rEnsp = $randomEnsp;
	    $rEnsp =~ s/^\_//g;
	    $selRandom{$rEnsp} = 1;
	} else {
	    $randomEnsp = $random{$ensp1};
	}

	my $rEnsp = $randomEnsp;
	$rEnsp =~ s/^\_//g;

	if(exists $cosmicMinScore->{$rEnsp}){
	    foreach my $dist (sort {$a<=>$b} keys %{$cosmicMinScore->{$rEnsp}}) {
		print "\t$randomEnsp\t$cosmicIntProt->{$rEnsp}->{$dist}\t$dist\t$cosmicMinScore->{$rEnsp}->{$dist}";
		last;
	    }
	} else {
	    print "\t$randomEnsp\t-\t4\t0";
	}
	
	my $distSum = 0;
	my $scoreSum = 0;
	my $repeats = 10;
	for(my $i = 0; $i<$repeats; $i++) {
	    my $randomEnsp = $ensp[int(rand($#ensp + 1))];
	    if(exists $cosmicMinScore->{$randomEnsp}){
		foreach my $dist (sort {$a<=>$b} keys %{$cosmicMinScore->{$randomEnsp}}) {
		    $distSum += $dist;
		    $scoreSum += $cosmicMinScore->{$randomEnsp}->{$dist};
		    last;
		}
	    } else {
		$distSum += 4;
		$scoreSum += 0;
	    }
	}
	printf("\t%.0f\t%.0f", $distSum/$repeats, $scoreSum/$repeats);
    }

    if($step =~ /3/) {
	if(exists $keggCancer->{$ensg1}) {
	    print "\t$keggCancer->{$ensg1}";
	} else {
	    print "\t-";
	}
	if(exists $keggCancer->{$ensg2}) {
	    print "\t$keggCancer->{$ensg2}";
	} else {
	    print "\t-";
	}
    }

    if($step =~ /4/) {
	if(exists $kegg->{$ensp1}) {
	    print "\t$kegg->{$ensp1}";
	} else {
	    print "\t-";
	}
	if(exists $kegg->{$ensp2}) {
	    print "\t$kegg->{$ensp2}";
	} else {
	    print "\t-";
	}
    }
    print "\n";
}
print STDERR "ALL done\n" if($verbose);
 
