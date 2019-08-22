#!/usr/bin/perl -w
# Author: Abdullah Kahraman
# Date: 02.11.2016

###############################################################################
###############################################################################
### List for all ENSEMBL isoforms existing and missing interactions based   ###
### on STRING, DIMA information.                                            ###
###############################################################################
###############################################################################

use strict;
use warnings;
use Getopt::Long;
no autovivification;

my (
    # variable for parameters which are read in from commandline
    $help,
    $idFile,
    $domainFile,
    $stringSeqFile,
    $isoSeqFile,
    $stringIntFile,
    $domainIntFile,
    $pfamFile,
    $pfamFile2,
    $minIntScore,
    $minDomainIntScore,
    $verbose,
    );

##############################################################################
### read all needed parameters from commandline ##############################

&GetOptions(
    "help!"              => \$help,              # print this help
    "stringIdFile=s"     => \$idFile,            # e.g. items.proteins.v10.after_haploduplications.dump.gz
    "stringDomainFile=s" => \$domainFile,        # e.g. items.proteins_smartlinkouts.tsv.gz
    "stringIntFile=s"    => \$stringIntFile,     # e.g. 9606.protein.links.detailed.v10.txt.gz
    "stringSeqFile=s"    => \$stringSeqFile,     # e.g. 9606.protein.sequences.v10.fa.gz
    "ensemblSeqFile=s"   => \$isoSeqFile,        # e.g. Homo_sapiens.GRCh37.75.pep.all.fa.gz
    "dimaFile=s"         => \$domainIntFile,     # e.g. dima_network_all.tbl.gz
    "pfamFile=s"         => \$pfamFile,          # e.g. 9606.tsv.gz
    "pfamFile2=s"        => \$pfamFile2,         # e.g. ensp_pfam_v78.txt.gz 
    "minScore:i"         => \$minIntScore,       # minimum STRING combined score. Default is 500.
    "minDimaIntScore:i"  => \$minDomainIntScore, # minimum domain interaction score. Default is 2.
    "verbose!"           => \$verbose,           # print out additional information on calculation progress plus warning messages    
) or die "\nTry \"$0 -h\" for a complete list of options\n\n";

##############################################################################

# help
if ($help) {printHelp(); exit}

##############################################################################
### SETTINGS #################################################################
##############################################################################
$| = 1;
$minIntScore = 400 if(!defined $minIntScore);
my $org = 9606;
# miminum number of non-structural interaction evidence in DIMA 
$minDomainIntScore = 2 if(!defined $minDomainIntScore);
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

    $usage = "$0 "."-stringIdFile \$databases/string/v10.0/items.proteins.v10.after_haploduplications.dump.gz ".
	           "-stringDomainFile \$databases/string/v10.0/items.proteins_smartlinkouts.tsv.gz ".
		   "-stringSeqFile \$databases/string/v10.0/9606.protein.sequences.v10.fa.gz ".
		   "-stringIntFile \$databases/string/v10.0/9606.protein.links.detailed.v10.txt.gz ".
		   "-ensemblSeqFile \$databases/ensembl/v75/Homo_sapiens.GRCh37.75.pep.all.fa.gz ".
		   "-dimaFile \$databases/dima/dima_network_all.tbl.gz ".
		   "-pfamFile \$databases/pfam/9606_v29.0.tsv.gz ".
		   "-minScore 400 ".
		   "-v";

    print "\nUsage: ".$usage."\n\n";

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
###############################################################################
sub readIdFile {
###############################################################################
    my %ids = ();
    open(F2, "zcat $idFile |") or die "\nERROR: Failed to open $idFile: $!\n\n";
    while(my $l = <F2>) {
	next unless($l =~ /^\d/);
	chomp($l);
	my ($id, $stringIdP, $orgId, $protChecksum, $protSize, $annotation, $name, $annotationWords) = split(/\t/, $l);

	# focus only on humans
	if($orgId == $org) {
	    my $ensp = $stringIdP;
	    $ensp =~ s/^\d+\.//;
	    $ids{$id} = $ensp;
	}
    }
    close(F2);
    return \%ids;
}
###############################################################################
sub readPfamFile {
###############################################################################
    my ($file) = @_; 
 
    my %pfamName2acc = ();
    my %pfam = ();
    open(F9, "zcat $file |") or die "\nERROR: Failed to open $file: $!\n\n";
    while(my $l = <F9>) {
	next if($l =~ /^#/);
	chomp($l);
	
	my ($seqId, $alignmentStart, $alignmentEnd, $envelopeStart, $envelopeEnd, $hmmAcc, $hmmName, $type, $hmmStart, $hmmEnd, $hmmLength, $bitScore, $eValue, $significance, $clan)
	    = split(/\t/, $l);
	$pfamName2acc{$hmmName} = $hmmAcc;
	$hmmAcc =~ s/\..*//;
	$pfam{$seqId}->{"$hmmAcc:$hmmName"}->{"$alignmentStart-$alignmentEnd"} = 1;
    }
    close(F9);
    return (\%pfamName2acc, \%pfam);
}
###############################################################################
sub readDomainFile {
###############################################################################
    my $ids = &readIdFile();
    my ($pfamName2acc, $pfam) = &readPfamFile($pfamFile);

    my %domain = ();
    open(F5, "zcat $domainFile |") or die "\nERROR: Failed to open $domainFile: $!\n\n";
    while(my $l = <F5>) {
	next unless($l =~ /^\d/);
	chomp($l);
	my ($id, $protSize, $smartURL) = split(/\t/, $l);
	if(exists $ids->{$id}) {
	    my $ensp = $ids->{$id};

	    my $smart = $smartURL;
	    $smart =~ s/http\:\/\/smart\.embl\.de\/smart\/.*.cgi\?smart\=//;
	    my ($protSize2, $annotations) = split(/\:/, $smart);
	    my @annotations = split(/\+/, $annotations);
	    foreach my $annotation (@annotations) {
		if($annotation =~ /(.*)\((\d+)\|(\d+)\)/) {
		    my $domainName = $1;
		    my $start = $2;
		    my $end = $3;

		    if($domainName =~ /^Pfam_(.*)/) {
			my $pfamName = $1;
			if(exists $pfamName2acc->{$pfamName}) {
			    $domainName = $pfamName2acc->{$pfamName}.":".$pfamName;
			}
		    }

#		    print "$ensp\t$domainName\t$start-$end\n";
		    $domain{$ensp}->{$domainName}->{"$start-$end"} = 1;
		} else {
		    print STDERR "WARNING: Domain annotation format is unknown of $annotation\n"
			if($verbose);
		}
	    }
# ID file is parsed only for HUMAN proteins. This file holds many non-human information, so many warning messages will be printed.
#	} else {
#	    print STDERR "WARNING: $id doesn't exist in STRING ID file $idFile!\n"
#		if($verbose);
	}
    }
    my ($pfamName2acc2, $pfam2) = &readPfamFile($pfamFile2);
    foreach my $ensp (sort keys %$pfam2) {
	foreach my $pfamName (sort keys %{$pfam2->{$ensp}}) {
	    foreach my $startEnd (sort keys %{$pfam2->{$ensp}->{$pfamName}}) {
		$domain{$ensp}->{$pfamName}->{$startEnd} = 1 if(!exists $domain{$ensp}->{$pfamName}->{$startEnd});
	    }
	}
    }

    close(F5);
    return \%domain;
}
###############################################################################
sub readStringSeqFile {
###############################################################################
    my %seq = ();
    open(F3, "zcat $stringSeqFile |") or die "\nERROR: Failed to open $stringSeqFile: $!\n\n";
    my $stringIdP = "";
    my $read = 0;
    while(my $l = <F3>) {
	next if($l =~ /^#/);
	chomp($l);
	if($l =~ /^>(.*)/) {
	    ($stringIdP, my $annot) = split(/\s+/, $1);
	    if($stringIdP =~ /^$org\./) {
		$read = 1;
		$stringIdP =~ s/^\d+\.//;
	    } else {
		$read = 0;
	    }
	    next;
	} else {
	    if($read == 1) {
		$seq{$stringIdP} .= $l;
# 		print "$stringIdP\t$l\n";
	    }
	}
    }
    close(F3);
    return \%seq;
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
    open(F6, "zcat $isoSeqFile |") or die "\nERROR: Failed to open $isoSeqFile: $!\n\n";
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
#	    print "$ensp\t$enst\t$ensg\t$l\n";
	}
    }
    close(F6);

    return (\%ensp2ensg, \%ensg2ensp,
	    \%enst2ensg, \%ensg2enst,
	    \%ensp2enst, \%enst2ensp,
	    \%ensp2seq);
}
###############################################################################
sub readStringIntFile {
###############################################################################
    my %int = ();

    open(F5, "zcat $stringIntFile |") or die "\nERROR: Failed to open $stringIntFile: $!\n\n";
    while(my $l = <F5>) {
	next if($l !~ /^\d/);
	chomp($l);
	my @a = split(/\s/, $l);
	my $ensp1 = $a[0];
	my $ensp2 = $a[1];
	my $score = $a[-1];
	$ensp1 =~ s/^\d+\.//;
	$ensp2 =~ s/^\d+\.//;

#	print "$ensp1\t$ensp2\t$score\n";
	$int{$ensp1}->{$ensp2} = $score if($score >= $minIntScore);
    }
    close(F5);
    return \%int;
}
###############################################################################
sub readDomainIntFile {
###############################################################################
    my %domainInt = ();

    open(F8, "zcat $domainIntFile |") or die "\nERROR: Failed to open $domainIntFile: $!\n\n";
    while(my $l = <F8>) {
	next if($l !~ /^PF/);
	chomp($l);
	my ($dom1, $dom2, $score, $ipfam, $did3, $dprofm, $dpea, $dpeaString, $cmm, $dipd) = split(/\t/, $l);

	next if($score < $minDomainIntScore and $ipfam eq "" and $did3 eq "");

	$domainInt{$dom1}->{$dom2} = $score;
	$domainInt{$dom2}->{$dom1} = $score;
#	print "$dom1\t$dom2\t$score\n";
    }
    close(F8);
    return \%domainInt;
}
###############################################################################
sub whichDomainsInteract {
###############################################################################
    my ($domains1, $domains2, $domainInt) = @_;

    # found out which domains are used for the interaction
    my $domain1 = "";
    my $domain2 = "";
    my $maxScore = 0;
    if(defined $domains1 and defined $domains1) {

	foreach my $domainInfo1 (sort keys %$domains1) {
	    next if($domainInfo1 !~ /^PF\d+/);
	    my ($domainId1, $domainName1) = split(/\:/, $domainInfo1);
	    foreach my $domainInfo2 (sort keys %$domains2) {
		next if($domainInfo2 !~ /^PF\d+/);
		my ($domainId2, $domainName2) = split(/\:/, $domainInfo2);
		if(exists $domainInt->{$domainId1}->{$domainId2} and $domainInt->{$domainId1}->{$domainId2} > $maxScore) {
		    $domain1 = $domainInfo1;
		    $domain2 = $domainInfo2;
		    $maxScore = $domainInt->{$domainId1}->{$domainId2};
		}
	    }
	}
    }
    return ($domain1, $domain2, $maxScore);
}

##############################################################################
### END OF SUBROUTINES########################################################
############################################################################## 


############
### MAIN ###
############

if(!defined $idFile) {
    print STDERR "\nPlease specify a STRING ID file. See $0 -help for more information.\n\n";
    exit;
}
if(!defined $domainFile) {
    print STDERR "\nPlease specify a domain file. See $0 -help for more information.\n\n";
    exit;
}
if(!defined $stringSeqFile) {
    print STDERR "\nPlease specify a STRING sequence file. See $0 -help for more information.\n\n";
    exit;
}
if(!defined $isoSeqFile) {
    print STDERR "\nPlease specify an ENSEMBL isoform sequence file. See $0 -help for more information.\n\n";
    exit;
}
if(!defined $stringIntFile) {
    print STDERR "\nPlease specify a STRING interaction file. See $0 -help for more information.\n\n";
    exit;
}
if(!defined $domainIntFile) {
    print STDERR "\nPlease specify a DIMA domain interaction file. See $0 -help for more information.\n\n";
    exit;
}
if(!defined $pfamFile) {
    print STDERR "\nPlease specify a PFAM domain DB file. See $0 -help for more information.\n\n";
    exit;
}
if(!defined $pfamFile2) {
    print STDERR "\nPlease specify a self created PFAM domain DB file. See $0 -help for more information.\n\n";
    exit;
}

print STDERR "Reading $isoSeqFile ... " if($verbose);
my ($ensp2ensg, $ensg2ensp, $enst2ensg, $ensg2enst, $ensp2enst, $enst2ensp, $ensp2seq) = &readSeqFile();
print STDERR "done\n" if($verbose);
print STDERR "Reading $stringSeqFile ... " if($verbose);
my $stringEnsp2seq = &readStringSeqFile();
print STDERR "done\n" if($verbose);
print STDERR "Reading $stringIntFile ... " if($verbose);
my $int = &readStringIntFile();
print STDERR "done\n" if($verbose);
print STDERR "Reading domain file ... " if($verbose);
my $domains = &readDomainFile();
print STDERR "done\n" if($verbose);
print STDERR "Reading $domainIntFile ... " if($verbose);
my $domainInt = &readDomainIntFile();
print STDERR "done\n" if($verbose);
print STDERR "Starting to analyse interactions in isoforms ... \n";



print "#ENSG\tENSP\tENST\tSTRINGensp\tSTRINGintN\tExistIntN\tMissIntN\tRelMissIntN\tExistInts\tMissInts\n";
foreach my $ensg1 (sort keys %$ensg2ensp) {
    foreach my $ensp1 (sort keys %{$ensg2ensp->{$ensg1}}) {
	my $existInt = "";
	my $missInt  = "";
	my $stringIntN = 0;
	my $existIntN = 0;
	my $missIntN = 0;

	my $enst1 = $ensp2enst->{$ensp1};
	my $stringEnsp1;
	foreach my $ensp (sort keys %{$ensg2ensp->{$ensg1}}) {
	    if(exists $stringEnsp2seq->{$ensp}) {
		$stringEnsp1 = $ensp;
		last;
	    }
	}
	if(!defined $stringEnsp1) {
	    print STDERR "WARNING: No canonical ENSP Id found for $ensg1:$ensp1 in $stringSeqFile. Skipping protein.\n";
	    next;
	}
	# skip this STRING protein if it does not exist in ENSEMBL or if its sequence is different in ENSEMBL
	if($stringEnsp2seq->{$stringEnsp1} ne $ensp2seq->{$stringEnsp1}) {
	    print STDERR "WARNING: 1 $stringEnsp1 has different sequence in $stringSeqFile and $isoSeqFile. Skipping protein.\n";
	    next;
	}

	# interaction of ensp1 is defined via its canonical isoform in STRING. 
	foreach my $stringEnsp2 (sort keys %{$int->{$stringEnsp1}}) {
	    next if(!exists $ensp2seq->{$stringEnsp2});
	    if($stringEnsp2seq->{$stringEnsp2} ne $ensp2seq->{$stringEnsp2}) {
		print STDERR "WARNING: 2 $stringEnsp2 has different sequence in $stringSeqFile and $isoSeqFile. Skipping protein.\n";
		next;
	    }

	    $stringIntN++;
	    my ($domain1, $domain2, $maxScore) = &whichDomainsInteract($domains->{$stringEnsp1}, $domains->{$stringEnsp2}, $domainInt);
	    if($domain1 ne "") {
		# fantastic ... interacting domains are known for this interaction. Check now whether domains exist in transcripts
		my ($domainId1, $domainName1) = split(/\:/, $domain1);
		my $found = 0;
		my $anno = "";
		foreach my $startEnd1 (sort keys %{$domains->{$stringEnsp1}->{$domain1}}) {
		    (my $start, my $end) = split(/-/, $startEnd1);
		    my $domainSeq = substr($ensp2seq->{$stringEnsp1}, $start-1, $end-$start+1);
		    $anno = "$domain1:$start-$end:$domain2:$stringEnsp2:$int->{$stringEnsp1}->{$stringEnsp2}:$maxScore";
		    if($ensp2seq->{$ensp1} =~ /$domainSeq/) {
			$found = 1;
			last;
		    }		    
		}
		if($found) {
		    $existInt .= ",$anno";
		    $existIntN++;
		} else {
		    $missInt .= ",$anno";
		    $missIntN++;
		}
	    }
	}
	if($existInt eq "") {
	    $existInt = "-";
	} else {
	    $existInt =~ s/^\,//g;
	}
	if($missInt eq "") {
	    $missInt = "-";
	} else {
	    $missInt =~ s/^\,//g;
	}
	print "$ensg1\t$ensp1\t$enst1\t$stringEnsp1\t$stringIntN\t$existIntN\t$missIntN\t".($existIntN+$missIntN > 0 ? sprintf("%.3f", $missIntN/($existIntN+$missIntN)) : "-")."\t$existInt\t$missInt\n";
    }
}
print STDERR "DONE\n";
