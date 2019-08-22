#!/usr/bin/perl -w
# Author: Abdullah Kahraman
# Date: 06.05.2016

###############################################################################
###############################################################################
### Outputs relative transcript expression per gene.                        ###
###############################################################################
###############################################################################

use strict;
use warnings;
use Getopt::Long;
no autovivification;

my (
    # variable for parameters which are read in from commandline
    $help,
    $metaFile,
    $expressFile,
    $stringSeqFile,
    $seqFile,
    $gffFile,
    $org,
    $out,
    $verbose,
   );

##############################################################################
### read all needed parameters from commandline ##############################

&GetOptions(
    "help!"                => \$help,               # print this help
    "metaFile=s"           => \$metaFile,           # e.g. rnaseq_metadata.tsv.gz
    "expressFile=s"        => \$expressFile,        # e.g. freeze3_v2.kallisto.lib.trans.fpkm.wl.aliquot_id.tsv.gz
    "stringSeqFile=s"      => \$stringSeqFile,      # e.g. 9606.protein.sequences.v10.fa.gz
    "seqFile=s"            => \$seqFile,            # e.g. Homo_sapiens.GRCh37.75.pep.all.fa.gz
    "gffFile:s"            => \$gffFile,            # e.g. Homo_sapiens.GRCh37.75.gff   
    "org:i"                => \$org,                # e.g. default is 9606 
    "out!"                 => \$out,                # outputs statistics on single canonical isoforms rather than summary statistics on cancer types
    "verbose!"             => \$verbose,            # print out additional information on calculation progress plus warning messages
) or die "\nTry \"$0 -h\" for a complete list of options\n\n";

##############################################################################

# help
if ($help) {printHelp(); exit}

##############################################################################
### SETTINGS #################################################################
##############################################################################
$| = 1;
$org = 9606 if(!defined $org);

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
    $sourceCode =~ s/.*\&GetOptions\(//s;
    $sourceCode =~ s/\).+//s;
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

###############################################################################
sub readStringSeqFile {
###############################################################################
    my %ensp2seq = ();
    open(F1, "zcat $stringSeqFile |") or die "\nERROR: Failed to open $stringSeqFile: $!\n\n";
    my $ensp = "";
    my $read = 0;

    while(my $l = <F1>) {
	next if($l =~ /^#/);
	chomp($l);
	if($l =~ /^>/) {
	    $ensp = $l;
	    if($ensp =~ /^>$org\./) {
		$ensp =~ s/^>$org\.//;
		$read = 1;
	    } else {
		$read = 0;
	    }
	} else {
	    if($read == 1) {
		$ensp2seq{$ensp} .= $l;
	    }
	}
    }
    close(F1);

    return \%ensp2seq;
}
###############################################################################
sub readMetaDataFile {
    my %sample2cancerType = (); 
    open(F3, "zcat $metaFile|") or die "\nERROR: Failed to open $metaFile: $!\n";
    while(my $l = <F3>) {
	next if($l =~ /^#/);
	chomp($l);
	my @a = split(/\t/, $l);
	$sample2cancerType{$a[13]} = $a[1];
    }
    close(F3);
    return \%sample2cancerType;
}
###############################################################################
sub readExpressionFile {
    (my $enst2ensg, my $gff) = @_;
    my $sample2cancerType = &readMetaDataFile();

    my %express = ();
    my %totalExpress = ();
    open(F4, "zcat $expressFile|") or die "\nERROR: Failed to open $expressFile: $!\n";
    my $n = 0;
    my @samples = ();
    while(my $l = <F4>) {
	chomp($l);
	if($n == 0) {
	    @samples = split(/\t/, $l);
	    $n++;
	    next;
	}
	my @a = split(/\t/, $l);
	$a[0] =~ s/\..*//;
	$a[0] =~ s/ENSTR/ENST0/;
	my %ensts = ($a[0], 1);
	if(defined $gff) {
	    if(exists $gff->{$a[0]}) {
		%ensts = %{$gff->{$a[0]}};
	    } else {
		print STDERR "WARNING: $a[0] does not exist in GFF file. Please check in $gffFile!\n";
	    }
	}
	# had to use a enst hash to allow for multiple ensts from GFF file
	foreach my $enst (sort keys %ensts) {
	    for(my $i = 1; $i < @a; $i++) {
		if(exists $sample2cancerType->{$samples[$i]}) {
		    my $cancerType = $sample2cancerType->{$samples[$i]};
		    my $expression = $a[$i];
		    $expression = 0 if($expression eq "NA");
		    if($expression > 0 and exists $enst2ensg->{$enst}) {
			if(exists $express{$cancerType}->{$samples[$i]}->{$enst}) {
			    $express{$cancerType}->{$samples[$i]}->{$enst} += $expression;
			} else {
			    $express{$cancerType}->{$samples[$i]}->{$enst} = $expression;
			}
			if(exists $totalExpress{$cancerType}->{$samples[$i]}->{$enst2ensg->{$enst}}) {
			    $totalExpress{$cancerType}->{$samples[$i]}->{$enst2ensg->{$enst}} += $expression;
			} else {
			    $totalExpress{$cancerType}->{$samples[$i]}->{$enst2ensg->{$enst}} = $expression;
			}
		    }
		} else {
		    print STDERR "WARNING: Sample $samples[$i] could not be found in $metaFile\n" if($verbose);
		}
	    }
	}
    }
    close(F4);
    return (\%express, \%totalExpress);
}
###############################################################################
sub readSeqFile {
###############################################################################
    my %ensp2seq = ();
    my %ensp2enst = ();
    my %enst2ensp = ();
    my %ensp2ensg = ();
    my %ensg2ensp = ();
    my %enst2ensg = ();
    my %ensg2enst = ();
    open(F6, "zcat $seqFile |") or die "\nERROR: Failed to open $seqFile: $!\n\n";
    my $ensp = "";
    my $ensg = "";
    my $enst = "";
    while(my $l = <F6>) {
	next if($l =~ /^#/);
	chomp($l);
	if($l =~ /^>/) {
	    $ensp = $l;
	    $ensg = $l;
	    $enst = $l;

	    $ensp =~ s/^>(.*?) .*/$1/;
	    $ensg =~ s/.* gene:(.*?) .*/$1/;
	    $enst =~ s/.* transcript:(.*?) .*/$1/;

	    $ensp2ensg{$ensp} = $ensg;
	    $ensp2enst{$ensp} = $enst;

	    $enst2ensp{$enst} = $ensp;
	    $enst2ensg{$enst} = $ensg;

	    $ensg2ensp{$ensg}->{$ensp} = $enst;
	    $ensg2enst{$ensg}->{$enst} = $ensp;
	} else {
	    $ensp2seq{$ensp} .= $l;
	}
    }
    close(F6);

    return (\%ensp2ensg, \%ensg2ensp,
	    \%enst2ensg, \%ensg2enst,
	    \%ensp2enst, \%enst2ensp,
	    \%ensp2seq);
}
###############################################################################
sub readGffFile {
###############################################################################
    my %gff = ();
    open(F1, $gffFile) or die "\nERROR: Failed to read $gffFile: $!\n\n";
    
    while(my $l = <F1>) {
	next if($l =~ /^#/);
	chomp($l);
	my @a = split(/\t/, $l);
	my $anno = $a[8];
	if($anno =~ /^transcripts/) {
	    $anno =~ s/\"//g;
	    $anno =~ s/\;//g;
	    my ($label, $ensts, $exonLabel, $exon, $geneLabel, $ensg) = split(/\s/, $anno);
	    my @ensts = split(/\+/, $ensts);
	    foreach my $enst(@ensts) {
		$gff{"$ensg:$exon"}->{$enst} = "";
	    }
	}
    }
    return \%gff;
}

###############################################################################
sub printHeader {
###############################################################################
    print "#CancerType\t".
	"AliquotID\t".
	"ENSG\t".
	"StringEnsp\t".
	"ENSP\t".
	"ENST\t".
	"TotalGeneExpression\t".
	"IsoformExpression\t".
	"IsoformExpressionInPercent\n";
}
##############################################################################
### END OF SUBROUTINES########################################################
############################################################################## 


############
### MAIN ###
############

if(!defined $metaFile) {
    print STDERR "\n\tPlease provide an META data file for the expression file. Try $0 -help to get more information\n\n";
    exit;
}
if(!defined $expressFile) {
    print STDERR "\n\tPlease provide a expression table file. Try $0 -help to get more information\n\n";
    exit;
}
if(!defined $stringSeqFile) {
    print STDERR "\nPlease specify a STRING sequence file. See $0 -help for more information.\n\n";
    exit;
}
if(!defined $seqFile) {
    print STDERR "\nPlease specify a sequence file. See $0 -help for more information.\n\n";
    exit;
}

my $gff = ();
$gff = &readGffFile() if(defined $gffFile);
my ($ensp2ensg, $ensg2ensp, $enst2ensg, $ensg2enst, $ensp2enst, $enst2ensp, $ensp2seq) = &readSeqFile();
my ($express, $totalExpress) = &readExpressionFile($enst2ensg, $gff);
my $stringEnsp2seq = &readStringSeqFile();

&printHeader();

foreach my $cancerType (sort keys %$express) {
    foreach my $sample (sort keys %{$express->{$cancerType}}) {

	foreach my $stringEnsp (sort keys %$stringEnsp2seq) {
	    # skip this STRING protein if it does not exist in ENSEMBL 
	    next if(!exists $ensp2ensg->{$stringEnsp});
	    next if($stringEnsp2seq->{$stringEnsp} ne $ensp2seq->{$stringEnsp});

	    my $ensg = $ensp2ensg->{$stringEnsp};	    
	    my $totalGeneExpression = 0;
	    $totalGeneExpression = $totalExpress->{$cancerType}->{$sample}->{$ensg} if(exists $totalExpress->{$cancerType}->{$sample}->{$ensg});
	    next if($totalGeneExpression == 0);

	    foreach my $ensp (sort keys %{$ensg2ensp->{$ensg}}) {
		if(exists $ensp2enst->{$ensp}) {
		    my $enst = $ensp2enst->{$ensp};

		    my $isoformExpression = 0;
		    $isoformExpression = $express->{$cancerType}->{$sample}->{$enst} if(exists $express->{$cancerType}->{$sample}->{$enst});
		    printf("%s\t%s\t%s\t%s\t%s\t%s\t%.3f\t%.3f\t%.3f\n",
			   $cancerType,
			   $sample,
			   $ensg,
			   $stringEnsp,
			   $ensp,
			   $enst,
			   $totalGeneExpression,
			   $isoformExpression,
			   $isoformExpression/$totalGeneExpression);
		}
	    }
	}
    }
}
