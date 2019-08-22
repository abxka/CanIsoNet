#!/usr/bin/perl -w
# Author: Abdullah Kahraman
# Date: 14.02.2014

###############################################################################
###############################################################################
### finds the connecting neighbours for two proteins in STRING.             ###
###############################################################################
###############################################################################

use strict;
use warnings;
use Getopt::Long;
use List::Util qw[min max sum];
use Parallel::ForkManager;

$|=1;

my (
    # variable for parameters which are read in from commandline
    $help,
    $infile,
    $col,
    $stringFile,
    $ensemblAliasesFile,
    $minStringScore,
    $maxNeighbourhood,
    $cpu,
    $stringId,
    $verbose,

    # global variables
    %string,
   );


##############################################################################
### read all needed parameters from commandline ##############################

&GetOptions(
    "help!"        => \$help,               # print this help.
    "id=s"         => \$stringId,           # protein STRING id 1, e.g. ENSP00000288602
    "infile=s"     => \$infile,             # file with STRING ids
    "col:i"        => \$col,                # column number in tab delimted file with ENSP information.
    "stringFile=s" => \$stringFile,         # STRING DB e.g. human_protein.links.detailed.v9.1.txt.gz
    "alias=s"      => \$ensemblAliasesFile, # aliases for ENSEMBL ids, downloaded from STRING e.g. human_protein.aliases.v9.1.txt.gz
    "min:i"        => \$minStringScore,     # minimum confidence score of STRING interaction, default is 400.
    "max:i"        => \$maxNeighbourhood,   # maximum neighbourhood size to search for the interaction partner, default is 1.
    "cpu:i"        => \$cpu,                # number of CPU cores to run calculations on. Default is 8.
    "verbose!"     => \$verbose,            # output additional information on computational progess.
) or die "\nTry \"$0 -h\" for a complete list of options\n\n";

##############################################################################

# help
if ($help) {printHelp(); exit}

$minStringScore ||= 400;
$maxNeighbourhood ||= 1;
$cpu ||= 8;

print "#minStringScore\tmaxNeighbourhood\tcpu\n".
      "#$minStringScore\t$maxNeighbourhood\t$cpu\n" if($verbose);
##############################################################################
### SETTINGS #################################################################
##############################################################################
my $interactionSign = "<=>";
my $pm = Parallel::ForkManager->new($cpu);
$pm->{auto_cleanup} = 1;

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
sub readIdFile {
    (my $path) = @_;
    open(F, $path) or die "ERROR: Failed to open $path!: $!\n";

    my @ids = ();
    while(<F>){
	chomp($_);
	my @a = split(/\t/, $_);
	push(@ids, $a[$col]);
    }
    close(F);
    return \@ids;
}
############################################################################## 
sub readStringFile {
    open(F4, "zcat $stringFile|") or die "ERROR: Failed to open $stringFile\n";
    while(<F4>) {
	next if(/^#/);
	chomp($_);
	my @a=split(/\s+/);
	$a[0] =~ s/^\d+\.//;
	$a[1] =~ s/^\d+\.//;
	$a[0] =~ s/(ENS[A-Z0-9]+\d+)\..*/$1/;
	$a[1] =~ s/(ENS[A-Z0-9]+\d+)\..*/$1/;
	if(exists $string{$a[0]} and defined $string{$a[0]}->{$a[1]}){
	    print STDERR "WARNING: STRING interaction $a[0]$interactionSign$a[1] is NOT unique\n";
	} else {
	    if($a[-1] >= $minStringScore) {
		$string{$a[0]}->{$a[1]} = $a[-1];
	    }
	}
    }
    close(F4);
}
##############################################################################
sub readAliasesFile {
    open(F5, "zcat $ensemblAliasesFile|") or die "ERROR: Failed to open $ensemblAliasesFile\n";

    my %geneNames = ();
    my %proteinFunction = ();
    my %uniprotId = ();
    while(<F5>){
	next if(/^#/);
	chomp($_);
	my @a = split(/\t/);
	$a[0] =~ s/.*\.//; 
	$proteinFunction{$a[0]} = $a[1] if(!exists $proteinFunction{$a[0]} and $a[2]=~ /\sEnsembl_UniProt_DE/i);
	$geneNames{$a[0]} = $a[1] if(!exists $geneNames{$a[0]} and $a[2] =~ /\bEnsembl_UniProt_GN\b/i);
	$uniprotId{$a[0]} = $a[1] if(!exists $uniprotId{$a[0]} and $a[2] =~ /\sEnsembl_UniProt_ID$/i);
    }
    close(F5);
    return (\%geneNames, \%proteinFunction, \%uniprotId);
}
##############################################################################
sub getInteractors {
    (my $exclude, my $p1) = @_;

    my $interactors = $string{$p1};
    my %interactors;
    foreach my $p2 (sort {$interactors->{$b}<=>$interactors->{$a}} keys %$interactors) {
	next if($exclude =~ /\b$p2\b/);
	$interactors{$p2}->{int} = $interactors->{$p2};
	$interactors{$p2}->{min} = 0;
	$interactors{$p2}->{sum} = 0;
	if(exists $string{$p2}) {
	    $interactors{$p2}->{min} = keys %{$string{$p2}};
	    $interactors{$p2}->{sum} = keys %{$string{$p2}};
	}
    }
    return \%interactors;
}

##############################################################################
sub neighbourhoodSearch {
    my ($neighbourhood, $count) = @_;

    my $isDeadEnd = 1;

    $count++;
    print STDERR "\tChecking $count-th degree neighbourhood. Need to check ".(keys %$neighbourhood)." proteins ... " if($verbose);
    foreach my $protein (sort keys %$neighbourhood){
	my $lastProtein = $protein;
	$lastProtein =~ s/.*$interactionSign//;
	my %interactors = %{&getInteractors($protein, $lastProtein)};

	if (keys %interactors >= 1){
	    foreach my $interactor (sort keys %interactors) {
		$neighbourhood->{"$protein$interactionSign$interactor"}->{int} = min($neighbourhood->{$protein}->{int}, $interactors{$interactor}->{int});
		$neighbourhood->{"$protein$interactionSign$interactor"}->{min} = min($neighbourhood->{$protein}->{min}, $interactors{$interactor}->{min});
		$neighbourhood->{"$protein$interactionSign$interactor"}->{sum} = sum($neighbourhood->{$protein}->{sum}, $interactors{$interactor}->{sum});
	    }
	    $isDeadEnd = 0;
	}
    }
    print STDERR "Done\n" if($verbose);

    if($isDeadEnd) {
	my %a = ("", "-1");
	return \%a;
    }
    if($count == $maxNeighbourhood) {
	return $neighbourhood;
    }
    &neighbourhoodSearch($neighbourhood, $count);
}

##############################################################################
sub neighbourhood {
    (my $alias) = @_;
    my %h = ();
    $h{$alias}->{int} = 1000;
    $h{$alias}->{min} = keys %{$string{$alias}};
    $h{$alias}->{sum} = keys %{$string{$alias}};
    return &neighbourhoodSearch(\%h, 0);
}

############################################################################## 
sub outputHeader {
    print "#Distance\tProteinId\tNeighbour\tPath\tProteinNumberOfInteractions\tPathMinScore\tPathNumberOfInteractionsSum\tPathNumberOfInteractionsMin\t".
	  "ProteinGeneName\tNeighbourGeneName\tProteinUniprotId\tNeighbourUniprotId";
#    print "\tProteinFunction\tNeighbourFunction";
    print "\n";
}
############################################################################## 
sub output {
    my ($protein, $neighbourhood, $geneName, $proteinFunction, $uniprotId) = @_;

    foreach my $interactor (sort {
	split(/$interactionSign/, $b) <=> split(/$interactionSign/, $a)
			    } keys %$neighbourhood) {
	my @neighbourhood = split(/$interactionSign/, $interactor);
	my $depth = @neighbourhood - 1;

	my $gn = "-";
	my $pf = "-";
	my $id = "-";
	$gn = $geneName->{$protein} if(exists $geneName->{$protein});
	$pf = $proteinFunction->{$protein} if(exists $proteinFunction->{$protein});
	$id = $uniprotId->{$protein} if(exists $uniprotId->{$protein});

	print "$depth\t$protein\t$neighbourhood[-1]\t$interactor\t".
	    (keys %{$string{$protein}})."\t".
	    "$neighbourhood->{$interactor}->{int}\t".
	    "$neighbourhood->{$interactor}->{sum}\t".
	    "$neighbourhood->{$interactor}->{min}\t";
	print "$gn\t";
	print exists $geneName->{$neighbourhood[-1]} ? $geneName->{$neighbourhood[-1]}."\t" : "-\t";
	print "$id\t";
	print exists $uniprotId->{$neighbourhood[-1]} ? $uniprotId->{$neighbourhood[-1]} : "-";
#	print "$pf\t";
#	print exists $proteinFunction->{$neighbourhood[-1]} ? $proteinFunction->{$neighbourhood[-1]}."\n" : "-\n";
	print "\n";
    }
}

############################################################################## 

$pm->run_on_finish(
    sub { my ($pid, $exit_code,  $ident, $exit_signal, $core_dump, $ng) = @_;
	  print STDERR "$pid has finished:  $exit_code\n";
	  my $alias = $ng->[0];
	  my $neighbourhood = $ng->[1];
	  my $geneName = $ng->[2];
	  my $proteinFunction = $ng->[3];
	  my $uniprotId = $ng->[4];

	  if(defined $neighbourhood) {
	      &output($alias, $neighbourhood, $geneName, $proteinFunction, $uniprotId);
	  }
    }
    );

$pm->run_on_start(
    sub { my ($pid)=@_; print STDERR "$pid starts\n"; }
    );

$pm->run_on_wait(
    sub {print STDERR "Waiting 1 sec ...\n"}, 1
    );


############################################################################## 
### END OF SUBROUTINES########################################################
############################################################################## 


############
### MAIN ###
############
if( ((!defined $stringId) and (!defined $infile or !defined $col)) or !defined $stringFile or !defined $ensemblAliasesFile) {
    print "\nMissing parameters. Try \"$0 -h\" for a complete list of options.\n\n";
} else {
    my $n = 0;
    my @aliases = ();
    ###############################################################################################
    if(defined $infile) {
	print STDERR "Reading in STRING IDs from $infile ..." if(defined $verbose);
	@aliases = @{&readIdFile($infile)};
	$n = @aliases;
	print STDERR "Done\n" if($verbose);
    } elsif(defined $stringId) {
	@aliases = ($stringId);
	$n = 1;
    }
    $n += 2;
    ###############################################################################################
    print STDERR $n--.". Reading STRING interactions from $stringFile ... " if(defined $verbose);
    &readStringFile();
    print STDERR "Done\n" if(defined $verbose);
    ###############################################################################################
    print STDERR $n--.". Reading alias files $ensemblAliasesFile ... " if(defined $verbose);
    my ($geneNames, $proteinFunction, $uniprotId) = &readAliasesFile();
    print STDERR "Done\n" if(defined $verbose);

    ###############################################################################################
    &outputHeader();
    ###############################################################################################
    foreach my $alias (@aliases) {
	my $pid = $pm->start and next; 
	if(!exists $string{$alias}) {
	    print STDERR "WARNING: ".$n--." Skipping $alias as it does not exist in the STRING database\n"
		if($verbose);
	    $pm->finish;
	    next;
	}

	print STDERR $n--.". Searching neighbourhood of $alias ...\n" if(defined $verbose);
	$pm->finish(0, [$alias, &neighbourhood($alias), $geneNames, $proteinFunction, $uniprotId]);
    }
    ###############################################################################################
    $pm->wait_all_children;
}
