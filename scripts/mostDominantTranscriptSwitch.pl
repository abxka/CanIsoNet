#!/usr/bin/perl -w
# Author: Abdullah Kahraman
# Date: 06.02.2018

###############################################################################
###############################################################################
### Determines MDI based on 2-fold expression difference to 2nd MDI and     ###
### subsequently compares the relative expression difference to expression  ###
### differences in GTEx using a Weibull distribution based P-value.         ### 
###############################################################################
###############################################################################

use strict;
use warnings;
use Getopt::Long;
no autovivification;
use Statistics::R;

my (
    # variable for parameters which are read in from commandline
    $help,
    $enst2ensgFile,
    $expressPcawgFile,
    $expressGtexFile,
    $isoformIntFile,
    $enstColumnNo,
    $minStringScore,
    $minEnrichment,
    $maxQvalue,
    $minExpress,
    $canonFile,
    $verbose,
    $pvalueFile,
    $sequenceFile,
   );

##############################################################################
### read all needed parameters from commandline ##############################

&GetOptions(
    "help!"                 => \$help,             # print this help
    "pcawg=s"               => \$expressPcawgFile, # e.g. databases/pcawg/expression/pcawg.rnaseq.transcript.expr.tpm.tsv.gz
    "gtex=s"                => \$expressGtexFile,  # e.g. databases/pcawg/gtex/GTEX_v4.pcawg.transcripts.tpm.tsv.gz
    "ensgFile=s"            => \$enst2ensgFile,    # e.g. ensg_ensp_enst_ense_geneName_v75.tsv.gz
    "canonFile=s"           => \$canonFile,        # e.g. /mnt/mnemo2/abdullah/databases/string/v10.0/canonical_isoforms_ensg_ensp_enst_geneName_v75.tsv.gz
    "isoformIntFile=s"      => \$isoformIntFile,   # isoform interaction file, e.g. interactionsInIsoforms_900_2.tsv.gz
    "sequenceFile=s"        => \$sequenceFile,     # file with sequences of ENSPs e.g. ensp_ensg_enst_sequence.tsv.gz
    "minStringScore=i"      => \$minStringScore,   # e.g. 900
    "enstColumnNuo=i"       => \$enstColumnNo,     # e.g. 0. Column number of enst in pcawg and gtex file. Starts with 0.
    "minEnrichment=f"       => \$minEnrichment,    # e.g. 2
    "maxQvalue=f"           => \$maxQvalue,        # e.g. 0.05
    "minExpress=f"          => \$minExpress,       # e.g. 1. Minimum Expression below which isoform is considered to be not expressed.
    "pvalueFile:s"          => \$pvalueFile,       # e.g. pvalues.tsv
    "verbose!"              => \$verbose,          # print out additional information on calculation progress plus warning messages
) or die "\nTry \"$0 -h\" for a complete list of options\n\n";

##############################################################################

# help
if ($help) {printHelp(); exit}

##############################################################################
### SETTINGS #################################################################
##############################################################################
$| = 1;
my $multFactor = 1000;
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
sub readEnst2ensgFile {
    my %enst2ensg = (); 
    open(F1, "zcat $enst2ensgFile|") or die "\nERROR: Failed to open $enst2ensgFile: $!\n";
    while(my $l = <F1>) {
	next if($l !~ /^ENS/);
	chomp($l);
	my @a = split(/\t/, $l);
	$enst2ensg{$a[2]} = $a[0];
    }
    close(F1);
    return \%enst2ensg;
}
###############################################################################
sub readExpressionFile {
    (my $file) = @_;
    
    my $enst2ensg = &readEnst2ensgFile();

    my %enstExpress = ();
    my %ensgExpress = ();
    open(F2, "zcat $file |") or die "\nERROR: Failed to open $file: $!\n";
    my $n = 0;
    my @samples = ();
    while(my $l = <F2>) {
	chomp($l);
	# read header with sample names
	if($n == 0) {
	    @samples = split(/\t/, $l);
	    $n++;
	    next;
	}

	my @a = split(/\t/, $l);
	my $enst = $a[$enstColumnNo];
	next if($enst !~ /ENST0/);
	$enst =~ s/\..*//;

	if(!exists $enst2ensg->{$enst}) {
	    print STDERR "WARNING: $enst does not exist in $enst2ensgFile. Skipping transcript.\n";
	    next;
	}

	my $ensg = $enst2ensg->{$enst};
	
	for(my $i = $enstColumnNo+1; $i < @a; $i++) {
	    my $expression = $a[$i];
	    $expression = 0 if($expression eq "NA");
	    $enstExpress{$ensg}->{$samples[$i]}->{$enst} = $expression;

	    if(exists $ensgExpress{$ensg}->{$samples[$i]}) {
		$ensgExpress{$ensg}->{$samples[$i]} += $expression;
	    } else {
		$ensgExpress{$ensg}->{$samples[$i]} = $expression;
	    }
	}
    }
    close(F2);

    my %relEnstExpress = ();
    my %relEnstExpressString = ();
    # calculate relative transcript expression
    foreach my $ensg (sort keys %enstExpress) {
	foreach my $sample (sort keys %{$enstExpress{$ensg}}) {
	    my $geneExpress = $ensgExpress{$ensg}->{$sample};
	    foreach my $enst (sort keys %{$enstExpress{$ensg}->{$sample}}) {
		my $transExpress = $enstExpress{$ensg}->{$sample}->{$enst};

		my $relExpress = 0.0;
		$relExpress = $transExpress / $geneExpress * $multFactor if($geneExpress > 0);

		$relEnstExpress{$ensg}->{$enst}->{$sample} = $relExpress;
		if(exists $relEnstExpressString{$ensg}->{$enst}) {
		    $relEnstExpressString{$ensg}->{$enst} .= sprintf(" %.3f", $relExpress);
		} else {
		    $relEnstExpressString{$ensg}->{$enst} = sprintf("%.3f", $relExpress); 
		}
	    }
	}
    }

    return (\%ensgExpress, \%enstExpress, \%relEnstExpress, \%relEnstExpressString);
}
##############################################################################
sub enrichment {

    my ($express, $minimumExpress) = @_;

    my %enrichment = ();
    foreach my $ensg (sort keys %$express) {
	foreach my $sample (sort keys %{$express->{$ensg}}) {
	    my $mostDominantEnst = "-";
	    my $mostDominantEnstExpress = 0;
	    my $enrichment = -1;
	    my $noEnst = keys %{$express->{$ensg}->{$sample}};
	    
	    # sort isoforms by highest expression
	    my $n = 0;
	    foreach my $enst (sort {$express->{$ensg}->{$sample}->{$b} <=> $express->{$ensg}->{$sample}->{$a}}
			      keys %{$express->{$ensg}->{$sample}}) {

		my $enstExpress = $express->{$ensg}->{$sample}->{$enst};

		# the first transcript is the one with the hightest expression, i.e. the one which is potentially the MDI. 
		if(++$n == 1) {
		    $mostDominantEnst = $enst;
		    $mostDominantEnstExpress = $express->{$ensg}->{$sample}->{$enst};
		    next;
		}
		my $enrichment = $mostDominantEnstExpress;
		# ignore transcripts with insignifcant expression
		if($mostDominantEnstExpress >= $minimumExpress){
		    $enrichment = $mostDominantEnstExpress/$enstExpress if($enstExpress > 0);
		    $enrichment{$ensg}->{$mostDominantEnst}->{$sample} = $enrichment if($enrichment >= $minEnrichment);
		}
		last if($n == 2);
	    }
	}
    }
    return \%enrichment;
}
###############################################################################
sub readCanonFile {
    my %canon = (); 
    open(F3, "zcat $canonFile|") or die "\nERROR: Failed to open $canonFile: $!\n";
    while(my $l = <F3>) {
	next if($l !~ /^ENS/);
	chomp($l);
	my ($canonEnsp, $ensg, $ensp, $enst, $geneName, $source) = split(/\t/, $l);
	$canon{$enst} = $canonEnsp;
	$canon{$ensg} = $canonEnsp;
	$canon{$ensp} = $ensg;
    }
    close(F3);
    return \%canon;
}
###############################################################################
sub readIsoformIntFile {
    my %missInt4isoform = ();
    my %intN = ();
    open(F4, "zcat $isoformIntFile |") or die "\nERROR: Failed to open $isoformIntFile: $!\n\n";
    while(my $l = <F4>) {
	next if($l !~ /^EN/);
	chomp($l);
	my ($ensg1, $ensp1, $enst1, $stringEnsp1, $stringIntN, $existIntN, $missIntN, $relMissIntN, $existInt, $missInt) = split(/\t/, $l);
	$intN{$enst1} = 0;
	
	my @int = split(/\,/, $missInt);
	foreach my $int (@int) {
	    if($int =~ /\:/) {
		my ($pfam1, $domain1, $range1, $pfam2, $domain2, $ensp2, $stringScore, $dimaScore) = split(/\:/, $int);
		next if($stringScore < $minStringScore);
		$missInt4isoform{$enst1}->{$int} = 1;
		$intN{$enst1}++;
	    }
	}
	@int = split(/\,/, $existInt);
	foreach my $int (@int) {
	    if($int =~ /\:/) {
		my ($pfam1, $domain1, $range1, $pfam2, $domain2, $ensp2, $stringScore, $dimaScore) = split(/\:/, $int);
		next if($stringScore < $minStringScore);
		$intN{$enst1}++;
	    }
	}
    }
    close(F4);
    return (\%missInt4isoform, \%intN);
}
###############################################################################
sub readSequenceFile {
    my %sequences = ();
    open(F5, "zcat $sequenceFile |") or die "\nERROR: Failed to open $sequenceFile: $!\n\n";
    while(my $l = <F5>) {
	next if($l !~ /^EN/);
	chomp($l);
	my ($ensp, $ensg, $enst, $sequence) = split(/\t/, $l);
	$sequences{$enst} = $sequence;
    }
    close(F5);
    print STDERR "Found ".(keys %sequences)." sequences\n" if(defined $verbose);
    return \%sequences;
}
##############################################################################
sub median {
    my @values = sort{$a<=>$b} @_;

    my $med = 0;
    # if odd number of elements take middle element. Array size is given in size-1
    if(@values%2==1){
	$med = $values[int(@values/2)];
    }
    # if even number of elements take average of middle two elements
    else{
	my $n= @values;	
	$med = ($values[int(@values/2)-1] + $values[int(@values/2)]) / 2;
    }
    return $med;
}
##############################################################################
### END OF SUBROUTINES########################################################
############################################################################## 


############
### MAIN ###
############

my $exit = 0;
if(!defined $expressPcawgFile) {
    print STDERR "\tPlease provide an expression table file for PCAWG samples. Try $0 -help to get more information\n";
    $exit = 1;
}
if(!defined $expressGtexFile) {
    print STDERR "\tPlease provide an expression table file for GTEx samples. Try $0 -help to get more information\n";
    $exit = 1;
}
if(!defined $enst2ensgFile) {
    print STDERR "\tPlease provide an Ensembl ID mapping data file. Try $0 -help to get more information\n";
    $exit = 1;
}
if(!defined $canonFile) {
    print STDERR "\tPlease provide a canonical isoform list file. Try $0 -help to get more information\n";
    $exit = 1;
}
if(!defined $isoformIntFile) {
    print STDERR "\tPlease provide an isoform interaction file. Try $0 -help to get more information\n";
    $exit = 1;
}
if(!defined $sequenceFile) {
    print STDERR "\tPlease provide a sequence file. Try $0 -help to get more information\n";
    $exit = 1;
}
if(!defined $enstColumnNo){
    print STDERR "\tPlease provide a column number for ENST IDs. Try $0 -help to get more information\n";
    $exit = 1;
}
if(!defined $minStringScore) {
    print STDERR "\tPlease provide a minimum STRING score. Try $0 -help to get more information\n";
    $exit = 1;
}
if(!defined $minEnrichment) {
    print STDERR "\tPlease provide a minimum fold difference value. Try $0 -help to get more information\n";
    $exit = 1;
}
if(!defined $maxQvalue) {
    print STDERR "\tPlease provide a maximum Q-value/P-value. Try $0 -help to get more information\n";
    $exit = 1;
}
if(!defined $minExpress) {
    print STDERR "\tPlease provide a minimum expression value. Try $0 -help to get more information\n";
    $exit = 1;
}
if($exit) {
    print "\n";
    exit;
}

print STDERR "Reading in $expressPcawgFile ... " if($verbose);
my ($ensgExpressPcawg, $enstExpressPcawg, $relEnstExpressPcawg, $relEnstExpressStringPcawg) = &readExpressionFile($expressPcawgFile);
print STDERR "done\n" if($verbose);

print STDERR "Reading in $expressGtexFile ... " if($verbose);
my ($ensgExpressGtex,  $enstExpressGtex, $relEnstExpressGtex, $relEnstExpressStringGtex)  = &readExpressionFile($expressGtexFile);
print STDERR "done\n" if($verbose);

print STDERR "Reading $canonFile ... " if($verbose);
my $canonical = &readCanonFile();
print STDERR "done\n" if($verbose);

print STDERR "Reading $isoformIntFile ... " if($verbose);
my ($missIsoInt, $intN) = &readIsoformIntFile();
print STDERR "done\n" if($verbose);

print STDERR "Reading $sequenceFile ... " if($verbose);
my ($sequences) = &readSequenceFile();
print STDERR "done\n" if($verbose);

print STDERR "Calculating Enrichment for PCAWG ... " if($verbose);
my $enrichmentPcawg = &enrichment($enstExpressPcawg, $minExpress);
print STDERR "done\n" if($verbose);

print STDERR "Calculating Enrichment for GTEx ... " if($verbose);
my $enrichmentGtex = &enrichment($enstExpressGtex, $minExpress*0.1);
print STDERR "done\n" if($verbose);

# print file header
print "#ENSG\tENSPcanon\tDomCancerTrans\tCancerSampleId\t".
      "GTExMDIs\ttotalRelevantGtexSamples\tMedianEnrichmentGTExMDIs\tPvalue\tQvalue\tNumberOfGtexMDIs\tEnrichmentMDIinCancer\t".
      "MedianRelMDIexpressionInNormal\tRelMDIexpressionInCancer\t".
      "UniqMissedInteractionsOfDomCancerTrans\tTotalNumberOfStringInt\t".
      "NumberOfUniqMissedInteractionsOfDomCancerTrans\tNumberOfCommonMissedInteractionsOfDomTrans\n";

my $R = Statistics::R->new();
$R->startR();

if(defined $pvalueFile) {
    open(O, ">$pvalueFile");
    print O "#ENSG\tDomCancerTrans\tCancerSampleId\tMedianRelMDIexpressionInNormal\tRelMDIexpressionInCancer\t".
	  "PCAWGhigherExpressionN\tPCAWGlowerExpressionN\tPvalue\n";
}


my @outputPart1 = ();
my @outputPart2 = ();
my $pvalues = "";
foreach my $ensg (sort keys %$enrichmentPcawg) {
    if(exists $relEnstExpressStringGtex->{$ensg}) {
	foreach my $enst (sort keys %{$enrichmentPcawg->{$ensg}}) {
	    if(exists $relEnstExpressStringGtex->{$ensg}->{$enst}) {
		foreach my $sample (sort keys %{$enrichmentPcawg->{$ensg}->{$enst}}) {
		    my $enrichmentPcawg = $enrichmentPcawg->{$ensg}->{$enst}->{$sample};
		    my $enstRelExpressPcawg = $relEnstExpressPcawg->{$ensg}->{$enst}->{$sample};
		    my $enstRelExpressesGtex = $relEnstExpressStringGtex->{$ensg}->{$enst};
		    # only do analysis for MDI, i.e. those transcripts having a two-fold expression surplus + their frequency is 
		    # significantly different from GTEx and they are specific to PCAWG + they are expressed more than minExpress
		    if(!exists $enrichmentGtex->{$ensg}->{$enst} and (exists $enrichmentGtex->{$ensg} and keys %{$enrichmentGtex->{$ensg}} > 0)) {

			my @enstRelExpressesGtex = split(/\s/, $enstRelExpressesGtex);
			
			my $relExpressionGtexMedian = 0;
			$relExpressionGtexMedian = &median(@enstRelExpressesGtex) if(@enstRelExpressesGtex > 0);

			my $pos = 0;
			my $neg = 0;
			foreach my $enstRelExpressGtex (@enstRelExpressesGtex) {
			    $pos++ if($enstRelExpressPcawg > $enstRelExpressGtex);
			    $neg++ if($enstRelExpressPcawg < $enstRelExpressGtex);
			}
			
			$R->send(qq`pcawg = $enstRelExpressPcawg`);
			$R->send(qq`test = binom.test(c($pos, $neg), p=0.5, alternative="two.sided", conf.level=1-$maxQvalue)`);
			$R->send(q`cat(test$p.value)`);
			my $pvalue = $R->read();

			print O "$ensg\t$enst\t$sample".sprintf("\t%.3f\t", $relExpressionGtexMedian/$multFactor).
			        sprintf("%.3f\t", $enstRelExpressPcawg/$multFactor)."$pos\t$neg".
				sprintf("\t%.3e\n", $pvalue) if(defined $pvalueFile);

			# So P-value must be smaller than maxQvalue
			if($pvalue < $maxQvalue) {
			    $pvalues .= $pvalue." ";

			    # in cancer, cancer driver isoforms should have missed interaction.
			    my $missedInt = "";
			    my $stringIntN = -1;
			    $stringIntN = $intN->{$enst} if(exists $intN->{$enst});
			    my $missedCancerIntN = -1;
			    my $missedCommonIntN = -1;
			    my %gtexMDIs = ();
			    my @mdiEnrichmentsGtex = (); 

			    # check how many GTEx MDIs exist for the gene. Only count a GTEx MDI if it has a different sequence. If it hasn't 
			    if(exists $enrichmentGtex->{$ensg}) {
				foreach my $nEnst (sort keys %{$enrichmentGtex->{$ensg}}) {
				    foreach my $nSample (sort keys %{$enrichmentGtex->{$ensg}->{$nEnst}}) {
					my $enrichmentGtex = $enrichmentGtex->{$ensg}->{$nEnst}->{$nSample};
					if($enrichmentGtex >= $minEnrichment and exists $sequences->{$nEnst} and exists $sequences->{$enst} and $sequences->{$nEnst} ne $sequences->{$enst}) {
					    $gtexMDIs{$nEnst}++;
					    push(@mdiEnrichmentsGtex, $enrichmentGtex);
					}
				    }
				}			
			    }
			    my $gtexMDIs = "";
			    # skip if GTEx doesn't have MDI for 50% of its cases.
			    my $skip = 1;
			    foreach my $g (sort{$gtexMDIs{$b}<=>$gtexMDIs{$a}} keys %gtexMDIs) {
				# check that sequence of MDI is different to sequence
				$gtexMDIs .= "$g:$gtexMDIs{$g},";
				$skip = 0 if($gtexMDIs{$g} >= ($#enstRelExpressesGtex+1)*0.5);
			    }
			    next if($skip);
			    $gtexMDIs =~ s/,$//g;
  			    my $medianMDIenrichmentGtex = "-";
			    $medianMDIenrichmentGtex = sprintf("%.3f", &median(@mdiEnrichmentsGtex)) if(@mdiEnrichmentsGtex > 0);

			    # check whether interactions are disrupted. If there are, check that interaction partner
			    # is expressed and that interactions are not already disrupted by GTEx MDIs. 
			    if(exists $missIsoInt->{$enst}) {
				# set values to 0 to indicate that interaction data is available for transcript
				$missedCancerIntN = 0;
				$missedCommonIntN = 0;
				
				foreach my $cInt (sort keys %{$missIsoInt->{$enst}}) {
				    # interaction partners must be expressed too in same cancer sample!
				    my ($pfam1, $domain1, $range1, $pfam2, $domain2, $ensp2, $stringScore, $dimaScore) = split(/\:/, $cInt);
				    if(exists $ensgExpressPcawg->{$canonical->{$ensp2}}->{$sample}) {
					my $found = 0;
					if(exists $enrichmentGtex->{$ensg}) {
					    foreach my $nEnst (sort keys %{$enrichmentGtex->{$ensg}}) {
						foreach my $nSample (sort keys %{$enrichmentGtex->{$ensg}->{$nEnst}}) {
						    my $enrichmentGtex = $enrichmentGtex->{$ensg}->{$nEnst}->{$nSample};
						    if($enrichmentGtex >= $minEnrichment) {
							$found = 1 if(exists $missIsoInt->{$nEnst}->{$cInt});
						    }
						}
					    }			
					}
					if($found) {
					    $missedCommonIntN++;
					} else {
					    # Only list disrupted interactions in cancer that are not disrupted in normal
					    $missedInt .= "$cInt,";
					    $missedCancerIntN++;
					}
				    }
				}
			    }
			    $missedInt =~ s/\,$//g;
			    $missedInt = "-" if($missedInt eq "");
			    				
			    push(@outputPart1,
				 sprintf("%s\t%s\t%s\t%s\t%s\t%i\t%s",
					 $ensg, (exists $canonical->{$enst} ? $canonical->{$enst} : "-"), $enst, $sample,
					 ($gtexMDIs eq "" ? "-" : $gtexMDIs), $#enstRelExpressesGtex+1, $medianMDIenrichmentGtex));
			    
			    push(@outputPart2,
				 sprintf("%i\t%.3f\t".
					 "%.3f\t%.3f\t".
					 "%s\t%i\t".
					 "%i\t%i",
					 $#mdiEnrichmentsGtex+1, $enrichmentPcawg,
					 $relExpressionGtexMedian/$multFactor, $enstRelExpressPcawg/$multFactor,
					 $missedInt, $stringIntN,
					 $missedCancerIntN, $missedCommonIntN));
			}
		    }
		}
	    } else {
		print STDERR "WARNING: $ensg - $enst does not exists in $expressGtexFile. Skipping it.\n";
		next;
	    }
	}
    } else {
	print STDERR "WARNING: $ensg does not exists in $expressGtexFile. Skipping $ensg.\n";
	    next;
    }
}
$pvalues =~ s/ $//g;
$R->send(qq`pvalues = unlist(strsplit("$pvalues", split=" "))`);
$R->send(qq`pvalues.corrected = p.adjust(pvalues, method = "fdr")`);
$R->send(qq`cat(pvalues.corrected)`);
my $qvalues = $R->read();
$R->stopR();
my @pvalues = split(/\s/, $pvalues);
my @qvalues = split(/\s/, $qvalues);

# TODO Q-value ausrechnen.
for(my $i=0; $i<@outputPart1; $i++){
    printf("%s\t%.3e\t%.3e\t%s\n", $outputPart1[$i], $pvalues[$i], $qvalues[$i], $outputPart2[$i]) if($qvalues[$i] < $maxQvalue);
}

close(O) if(defined $pvalueFile);
