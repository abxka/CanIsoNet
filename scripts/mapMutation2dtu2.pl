#!/usr/bin/perl -w
# Author: Abdullah Kahraman
# Date: 16.01.2017

###############################################################################
###############################################################################
### Maps mutations to missingInteractionsInDTU_filt_int_anno.tsv.gz files,  ###
### where a mutation is found in a functional region that is associated     ###
### with a gene in missingInteractionsInDTU_filt_int.tsv.gz.                ###
###############################################################################
###############################################################################

use strict;
use warnings;
use Getopt::Long;
no autovivification;

my (
    # variable for parameters which are read in from commandline
    $help,
    $dtuFile,
    $regionFile,
    $gtfFile,
    $mutFile,
    $metaFile,
    $annoFile,
    $eqtlFile,
    $enhancerFile,
    $geneNameCol,
    $maxStat,
    $projectCode,
    $verbose,
   );

##############################################################################
### read all needed parameters from commandline ##############################

&GetOptions(
    "help!"          => \$help,         # print this help
    "dtuFile=s"      => \$dtuFile,      # e.g. interactionDisruptionInDominantTranscripts_min2_filt_int_anno.tsv.gz
    "regionFile=s"   => \$regionFile,   # e.g. all_gc19_pc.bed
    "gtfFile:s"      => \$regionFile,   # to identify transcript specific UTR, promoter regions e.g. Homo_sapiens.GRCh37.75.gtf
    "mutFile=s"      => \$mutFile,      # e.g. October_2016_whitelist_2583.snv_mnv_indel.maf
    "metaFile:s"     => \$metaFile,     # relevant if input file is MDT based e.g. pcawg.rnaseq.extended.metadata.aliquot_id.V4.tsv.gz
    "annoFile=s"     => \$annoFile,     # e.g. promoter.annotation.tsv
    "eqtlFile=s"     => \$eqtlFile,     # e.g. promoter.cis.eQTL.all.tsv
    "enhancerFile:s" => \$enhancerFile, # e.g. map.enhancer.gene_max10.gz
    "geneCol=i"      => \$geneNameCol,  # column number of HGNC gene name. Starts with 0. 
    "maxStat=f"      => \$maxStat,      # maximum global sample size corrected p-value above which eQTLs are ignored.
    "project=s"      => \$projectCode,  # project name. Only mutations from this project will be mapped.
    "verbose!"       => \$verbose,      # output additional information on calculations
) or die "\nTry \"$0 -h\" for a complete list of options\n\n";

##############################################################################

# help
if ($help) {printHelp(); exit}

##############################################################################
### SETTINGS #################################################################
##############################################################################

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
sub readRegionFile {
    (my $genes) = @_;
     
    my $n= 0;
    my %region = ();
    open(F1, $regionFile) or die "\nERROR: Failed to open $regionFile: $!\n\n";
    my $header = <F1>;
    while(my $l = <F1>) {
	chomp($l);
	my ($chr, $start, $end, $name, $score, $strand, $thickStart, $thickEnd, $itemRgb,
	    $exonCount, $exonSizes, $exonStarts) = split(/\t/, $l);
	my @exonSizes = split(/,/, $exonSizes);
	my @exonStarts = split(/,/, $exonStarts);

	for(my $i=0; $i<$exonCount; $i++) {
	    my ($type, $db, $geneName, $ensg) = split(/\:+/, $name);
	    $ensg =~ s/\..*//;

	    $chr =~ s/chr//;
	    next if($ensg !~ /^ENS/ or !exists $genes->{$ensg});

	    my $rStart = $start+$exonStarts[$i];
	    my $rEnd = $start+$exonStarts[$i]+$exonSizes[$i];
	    for(my $j = $rStart; $j <= $rEnd; $j++) {
		$region{$ensg}->{$chr}->{$j}->{"$rStart-$rEnd"}->{$type} = 1;
	    }
	}
    }
    close(F1);
    return \%region;
}
##############################################################################
sub readMutFile {

    my %mut = ();
    open(F2, $mutFile) or die "\nERROR: Failed to open $mutFile: $!\n\n";
    my $header = <F2>;
    while(my $l = <F2>) {
	chomp($l);
	my @a = split(/\t/, $l);
	next if(@a < 43);
	my ($geneName, $chr, $start, $end, $strand, $varClass, $varType,
	    $RefAllele, $tumorSeqAllele1, $tumorSeqAllele2, $dbSNPrs, $dbSNPvalStatus,
	    $tumorSampleBarcode, $matchedNormSampleBarcode, $genomeChange, $refContext, $gcContent,
	    $i1000genomesAF, $i1000genomesID, $iCallers, $iGERM1000G, $iGERMOVLP, $iLOWSUPPORT, $iNORMALPANEL,
	    $iNumCallers, $iOXOGFAIL, $iREMAPFAIL, $iSEXF, $iVAF, $iBpcr, $iBseq, $iQual, $iRepeatMasker,
	    $iSignatureN3, $iSignatureR1, $iSignatureR2, $iSnvNearIndel, $tAltCount, $tRefCount, $iModelScore,
	    $iNvaf, $projectCode, $donorID) = split(/\t/, $l);

	for(my $i=$start; $i<=$end; $i++) {
	    $mut{$chr}->{$i}->{"$start-$end"}->{$projectCode}->{$tumorSampleBarcode}->{$varType} = 1;
	}
    }
    close(F2);
    return \%mut;
}
##############################################################################
sub readMetaFile {
    my %meta = ();
    open(F3, "zcat $metaFile |") or die "\nERROR: Failed to open $metaFile: $!\n\n";
    my $header = <F3>;
    while(my $l = <F3>) {
	chomp($l);
	my @a = split(/\t/, $l);
	my $rnaAliquotId = $a[0];
	my $wgsAliquotId = $a[-1];
	$meta{$rnaAliquotId} = $wgsAliquotId;
    }
    close(F3);
    return \%meta;
}
##############################################################################
sub readAnnoFile {
    my %anno = ();
    open(F1, $annoFile) or die "\nERROR: Failed to open $annoFile: $!\n\n";
    my $header = <F1>;
    while(my $l = <F1>) {
	chomp($l);
	my ($element, $chr, $start, $end, $pos) = split(/\t/, $l);
	$anno{$chr}->{$pos}->{$start}->{$end}->{$element} = 1;
    }
    close(F1);
    return \%anno;
}
##############################################################################
sub readEqtlFile {
    
    my ($genes) = @_;
    my $anno = &readAnnoFile();    
    my %eqtl = ();
    open(F2, $eqtlFile) or die "\nERROR: Failed to open $eqtlFile: $!\n\n";
    my $header = <F2>;
    while(my $l = <F2>) {
	chomp($l);
	my ($element, $type, $ensg, $chr, $pos, $pval, $lqval, $gpval, $beta, $diff)  = split(/\t/, $l);
	$ensg=~ s/\..*//;

	next if($gpval > $maxStat or !exists $genes->{$ensg});

	if(exists $anno->{$chr}->{$pos}) {
	    foreach my $start (sort keys %{$anno->{$chr}->{$pos}}) {
		foreach my $end (sort keys %{$anno->{$chr}->{$pos}->{$start}}) {
		    foreach my $elem (sort keys %{$anno->{$chr}->{$pos}->{$start}->{$end}}) {
			for(my $i=$start; $i<=$end; $i++) {
			    if($element eq $elem) {
				$eqtl{$ensg}->{$chr}->{$i}->{"$start-$end"}->{$elem}->{$type}->{$gpval}->{$beta} = 1;
			    } else {
				print STDERR "eQTL element $element differs from the one in the annotation file, $elem. Ignoring annotation file: \"$l\"\n";
				$eqtl{$ensg}->{$chr}->{$i}->{"$start-$end"}->{$element}->{$type}->{$gpval}->{$beta} = 1;
			    }
			}
		    }
		}
	    }
	} else {
	    print STDERR "WARNING $chr:$pos does not exist in $annoFile. Skipping eQTL \"$l\"\n";
	}
    }
    close(F2);
    return \%eqtl;
}
##############################################################################
sub readEnhancerFile {
    (my $regions) = @_;
    open(F6, "zcat $enhancerFile |") or die "\nERROR: Failed to open $enhancerFile: $!\n\n";
    while(my $l = <F6>) {
	chomp($l);
	my ($location, $genes) = split(/\t/, $l);
	my ($chr, $startEnd) = split(/:/, $location);
	$chr =~ s/chr//;
	my @ensg = split(/\;/, $genes);
	my ($start, $end) = split(/\-/, $startEnd);
	foreach my $ensg (@ensg) {
	    for(my $i=$start; $i<=$end; $i++){
		$regions->{$ensg}->{$chr}->{$i}->{$startEnd}->{"enhancer"} = 1;
	    }
	}
    }
    close(F6);
}
##############################################################################
sub getGenesFromDtuFile {
    my %genes = ();
    if($dtuFile =~ /.gz$/) {
	open(F4, "zcat $dtuFile |") or die "\nERROR: Failed to open $dtuFile: $!\n\n";
    } else {
	open(F4, $dtuFile) or die "\nERROR: Failed to open $dtuFile: $!\n\n";
    }
    while(my $l = <F4>) {
	next if($l =~ /^#/);
	chomp($l);
	
	my ($ensg, $enspCanon, $enstCanon, $enstIso, @a) = split(/\t/, $l);
	my $geneName = $a[$geneNameCol-4];
	$genes{$geneName} = 1;
	$genes{$ensg} = 1;
    }
    close(F4);
    return \%genes;
}

##############################################################################
### END OF SUBROUTINES########################################################
############################################################################## 


############
### MAIN ###
############

if(!defined $regionFile) {
    print STDERR "\nPlease provide a path to a BED file. Try $0 -help to get more information\n\n";
    exit;
}
if(!defined $mutFile) {
    print STDERR "\nPlease provide a path to a MAF file. Try $0 -help to get more information\n\n";
    exit;
}
if(!defined $dtuFile) {
    print STDERR "\nPlease provide a path to a MDI/DTU file. Try $0 -help to get more information\n\n";
    exit;
}
if(!defined $annoFile) {
    print STDERR "\nPlease provide a path to a annotation file. Try $0 -help to get more information\n\n";
    exit;
}
if(!defined $eqtlFile) {
    print STDERR "\nPlease provide a path to eQTL file. Try $0 -help to get more information\n\n";
    exit;
}
if(!defined $geneNameCol) {
    print STDERR "\nPlease provide a column number for HGNC gene names. Try $0 -help to get more information\n\n";
    exit;
}
if(!defined $projectCode) {
    print STDERR "\nPlease provide a project code. Try $0 -help to get more information\n\n";
    exit;
}
if(!defined $maxStat) {
    print STDERR "\nPlease provide a maximum threshold for eQTL p-value. Try $0 -help to get more information\n\n";
    exit;
}

print STDERR "6. Reading in DTU/MDI file ...\n" if($verbose);
my $genes = &getGenesFromDtuFile();
print STDERR "done\n" if($verbose);

print STDERR "4. Reading in region file ...\n" if($verbose);
my $regions = &readRegionFile($genes);
print STDERR "done\n" if($verbose);
if(defined $enhancerFile) {
    print STDERR "1. Adding enhancer regions to set of regions ...\n" if($verbose);
    &readEnhancerFile($regions); 
    print STDERR "done\n" if($verbose);
}

print STDERR "5. Reading in mutation file ...\n" if($verbose);
my $mutations = &readMutFile(); 
print STDERR "done\n" if($verbose);
print STDERR "3. Reading in eQTL file ...\n" if($verbose);
my $eqtl = &readEqtlFile($genes); 
print STDERR "done\n" if($verbose);
my $meta = ();
if(defined $metaFile) {
    print STDERR "2. Reading in meta data file ...\n" if($verbose);
    $meta = &readMetaFile(); 
    print STDERR "done\n" if($verbose);
}

print STDERR "1. Mapping mutation information ...\n" if($verbose);

# identify mutations only once for a sample and its gene
my %doneMuts = ();
my %doneEqtls = ();

if($dtuFile =~ /.gz$/) {
    open(F5, "zcat $dtuFile |") or die "\nERROR: Failed to open $dtuFile: $!\n\n";
} else {
    open(F5, $dtuFile) or die "\nERROR: Failed to open $dtuFile: $!\n\n";
}
while(my $l = <F5>) {
    chomp($l);
    if($l =~ /^#/) {
	print "$l\tMutationInfo\teQTLinfo\n";
	next;
    } 
    my ($ensg, $enspCanon, $enstCanon, $enstIso, @a) = split(/\t/, $l);
    my $geneName = $a[$geneNameCol-4];

    # 1st test functional regions for mutations
    my $muts = "";
    if(!exists $doneMuts{"$ensg-$enstIso"} ){
	if(exists $regions->{$ensg}) {
	    foreach my $chr (sort keys %{$regions->{$ensg}}) {
		my %done = ();
		foreach my $i (sort {$a<=>$b} keys %{$regions->{$ensg}->{$chr}}) {

		    if(exists $mutations->{$chr}->{$i}) {
			
			foreach my $rStartEnd (sort keys %{$regions->{$ensg}->{$chr}->{$i}}){
			    foreach my $mStartEnd (sort keys %{$mutations->{$chr}->{$i}}) {
				
				next if(exists $done{"$chr-$rStartEnd-$mStartEnd"});

				foreach my $pCode (sort keys %{$mutations->{$chr}->{$i}->{$mStartEnd}}) {
				    foreach my $tumorSampleBarcode (sort keys %{$mutations->{$chr}->{$i}->{$mStartEnd}->{$pCode}}) {

					next if(defined $metaFile and $tumorSampleBarcode ne $meta->{$enstIso});

					foreach my $varType (sort keys %{$mutations->{$chr}->{$i}->{$mStartEnd}->{$pCode}->{$tumorSampleBarcode}}) {
					    foreach my $regionType (sort keys %{$regions->{$ensg}->{$chr}->{$i}->{$rStartEnd}}){
						my $region = "$chr:$rStartEnd";
						my $mRegion = "$mStartEnd";
						my $mut = "$regionType:$region:$varType:$mRegion:$tumorSampleBarcode";
						$muts .= "$mut," if($muts !~ /\b$mut\b/);
						$done{"$chr-$rStartEnd-$mStartEnd"} = 1;
					    }
					}
				    }
				}
			    }
			}
		    }
		}
	    }
	}
    } else {
	$muts = $doneMuts{"$ensg-$enstIso"};
    }

    my $eqtls = "";
    # now test eQTL regions for mutations
    if(!exists $doneEqtls{"$ensg-$enstIso"}) {
	if(exists $eqtl->{$ensg}) {
	    foreach my $chr (sort keys %{$eqtl->{$ensg}}) {
		my %done = ();
		foreach my $i (sort {$a<=>$b} keys %{$eqtl->{$ensg}->{$chr}}) {
		    if(exists $mutations->{$chr}->{$i}) {
			
#			print "#2. $chr\t$i\n";
			foreach my $eStartEnd (sort keys %{$eqtl->{$ensg}->{$chr}->{$i}}){
			    foreach my $mStartEnd (sort keys %{$mutations->{$chr}->{$i}}) {
				
				next if(exists $done{"$chr-$eStartEnd-$mStartEnd"});

				foreach my $pCode (sort keys %{$mutations->{$chr}->{$i}->{$mStartEnd}}) {
				    foreach my $tumorSampleBarcode (sort keys %{$mutations->{$chr}->{$i}->{$mStartEnd}->{$pCode}}) {

#					print "#3. $i\t$eStartEnd\t$tumorSampleBarcode\n";
					next if(defined $metaFile and $tumorSampleBarcode ne $meta->{$enstIso});

					foreach my $varType (sort keys %{$mutations->{$chr}->{$i}->{$mStartEnd}->{$pCode}->{$tumorSampleBarcode}}) {
					    foreach my $eelement (sort keys %{$eqtl->{$ensg}->{$chr}->{$i}->{$eStartEnd}}){
						foreach my $eType (sort keys %{$eqtl->{$ensg}->{$chr}->{$i}->{$eStartEnd}->{$eelement}}){
						    foreach my $gpval (sort keys %{$eqtl->{$ensg}->{$chr}->{$i}->{$eStartEnd}->{$eelement}->{$eType}}){
							foreach my $beta (sort keys %{$eqtl->{$ensg}->{$chr}->{$i}->{$eStartEnd}->{$eelement}->{$eType}->{$gpval}}){
							    my $eRegion = "$chr:$eStartEnd";
							    my $mRegion = "$mStartEnd";
							    my $e = "$eelement:$eType:$eRegion:$gpval:$beta:$mRegion:$tumorSampleBarcode";
							    $eqtls .= "$e," if($eqtls !~ /\b$e\b/);
							    $done{"$chr-$eStartEnd-$mStartEnd"} = 1;
							}
						    }
						}
					    }
					}
				    }
				}
			    }
			}
		    }
		}
	    }
	}
    } else {
	$eqtls = $doneEqtls{"$ensg-$enstIso"};
    }
       
    $muts  =~ s/,$//g;
    $muts = "-" if($muts eq "");
    $eqtls  =~ s/,$//g;
    $eqtls = "-" if($eqtls eq "");
    $doneMuts{"$ensg-$enstIso"} = $muts;
    $doneEqtls{"$ensg-$enstIso"} = $eqtls;
    print "$l\t$muts\t$eqtls\n";
}
close(F5);

print STDERR "ALL DONE\n" if($verbose);
