# CanIsoNet
In the following, we are describing how an isoform-specific interaction network can be created using various PERL scripts and BASH commands to assess the pathological impact of alternatively spliced isoforms. 

The whole methodology follows following procedure:
![Overview of CanIsoNet methodology](https://)

Top of the figure with grey background shows the steps and filters for MDT switch detection based on the relative expression value analysis of transcripts in PCAWG and GTEx. Bottom of figure with violet background describes the methods and databases we used to develop CanIsoNet which combines functional interactions from STRING, with physical domain-domain interactions from 3did on all known transcripts from Ensembl. For the assessment of the functional and pathogenic impact, CanIsoNet counts the relative number of disrupted interactions, collects network density information from the STRING database and proximity information to genes to the COSMIC Cancer Gene Census within the STRING database. The middle section with white background indicates that the impact assessment of alternatively spliced isoforms is done by combining MDT switch information with data from CanIsoNet. 


Directory structure:

- bin
	- PfamScan
- databases
	- 3did/v2018_04 
	- ensembl/v75
	- string/v10
	- pcawg
	- pfam/v32
- scripts

## Isoform specific Protein-Protein Interactions
Required files are:

* [3did interaction file](https://3did.irbbarcelona.org/download/2018_04/3did_flat.gz). Save in directory databases/3did/v2018_04.
* [Ensembl transcript FASTA file](http://ftp.ensembl.org/Homo_sapiens.GRCh37.75.pep.all.fa.gz). Save in directory databases/ensembl/v75.
* Ensembl PFAM domain file. [Download from here](http://feb2014.archive.ensembl.org/biomart/martview). Save in directory databases/ensembl/v75.
* [PFAM domain file](ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam32.0/proteomes/9606.tsv.gz). Save in directory databases/pfam/v32.
* STRING haploduplications. Find in directory databases/string/v10.
* STRING SMART Linkouts. Fin in directory databases/string/v10.
* [STRING sequence file](http://version10.string-db.org/download/protein.sequences.v10/9606.protein.sequences.v10.fa.gz). Save in directory databases/string/v10.
* [STRING interaction file](http://version10.string-db.org/download/protein.links.detailed.v10/9606.protein.links.detailed.v10.txt.gz). Save in directory databases/string/v10.

<hr>
#### 
1. Prepare 3did domain interaction file
#### <a name="3did"></a>

	```
	# Filter out PFAM domain-domain interaction information 
	less databases/3did/v2018_04/3did_flat.gz | grep "^#=ID" | cut -f4,5 | perl -ane '$F[0]=~ s/.*(PF\d+).*/$1/; $F[1]=~ s/.*(PF\d+).*/$1/; print "$F[0]\t$F[1]\n$F[1]\t$F[0]\n"' | sort -u | gzip > 3did_pfamInteractions_1804.tsv.gz
	
	# Reformat to DIMA format. This is a legacy step 
	scripts/3did2dima.pl 3did_pfamInteractions_1804.tsv.gz | gzip > 3did_1804.tbl.gz
	```

2. Prepare PFAM domain identification
#### <a name="pfam"></a>

	```
	# Create directory structure for Ensembl protein FASTA files
	mkdir -p fasta/ensembl_v75/human
	
	# Create ENSEMBL protein FASTA files
	zless databases/ensembl/v75/Homo_sapiens.GRCh37.75.pep.all.fa.gz | perl -ne 'if(/^>(ENSP\d+) .*/){$i=$1; $h{$i}=$_; next} else{$h{$i}.=$_} END{foreach $i(sort keys %h){open(O, ">fasta/ensembl_v75/human/$i.fasta"); print O $h{$i}; close(O) }}'
	
	# Run pfam_scan.pl. HMM files of Pfam-A have been downloaded to 
	for s in fasta/ensembl_v75/human/*; do 	bin/PfamScan/pfam_scan.pl -fasta $s -dir bin/PfamScan | grep -v ^#; done | gzip > ensp_pfam_v78.txt.gz
	```

3. Create isoform-specific protein-protein interaction network.
#### <a name="isoformNetwork"></a>

	```
	# Get list of shared and missed interactions for each isoform
	scripts/missingDomainAndInteractions.pl -stringIdFile databases/string/v10/items.proteins.v10.after_haploduplications.dump.gz -stringDomainFile databases/string/v10/items.proteins_smartlinkouts.tsv.gz -stringSeqFile databases/string/v10/9606.protein.sequences.v10.fa.gz -stringIntFile databases/string/v10/9606.protein.links.detailed.v10.txt.gz -ensemblSeqFile databases/ensembl/v75/Homo_sapiens.GRCh37.75.pep.all.fa.gz -dimaFile 3did_1804.tbl.gz -pfamFile databases/pfam/v32/9606_v32.0.tsv.gz -minScore 900 -minDimaIntScore 2 -v -pfamFile2 ensp_pfam_v78.txt.gz >& interactionsInIsoforms_900_2_status.tsv && gzip interactionsInIsoforms_900_2_status.tsv
	
	# Extract interaction list from output file
	zgrep "^[#\|E]" interactionsInIsoforms_900_2_status.tsv.gz | gzip > interactionsInIsoforms_900_2.tsv.gz
	```

### Network density
#### <a name="networkDensity"></a>

```
scripts/networkDensity2.pl -in databases/string/v10/9606.protein.links.detailed.v10.txt.gz -min 900 -max 3 | gzip > string_v10_density_900_3.tsv.gz
```

### Breadth-first search for COSMIC Cancer Gene Census Genes
#### <a name="networkDensity"></a>

```
# Get all ENSP IDs for COSMIC CGC IDs
scripts/merge.pl -in1 databases/ensembl/v75/ensg_ensp_enst_geneName_v75.tsv.gz -col1 3 -in2 databases/cosmic/GRCh37/cancer_gene_census-grch37_v89.tsv -col2 0 -ignore | sort -u | grep ENSP > cosmic_census_grch37_v89_idMap_ensp.tsv

# Search in STRING for nearest neighbours in COSMIC CGC genes  
scripts/stringNeighbourhood.pl -in cosmic_census_grch37_v89_idMap_ensp.tsv -col 1 -string databases/string/v10/human.protein.network.string.v10.tsv.gz -min 900 -max 3 -alias databases/string/v10/human_protein.aliases.v10.txt.gz -cpu 20 | gzip > string_v10_cosmic_grch37_v89_neighbourhood_min900_shell3.tsv.gz
```
	

## MDT Detection
#### <a name="requiredFilesMDT"></a>
Required files are:

* [RNAseq metadata file](http://synapse.org/#!Synapse:syn7416381). Save in directory databases/pcawg/histo.
* [GTEx metadata file](http://synapse.org/#!Synapse:syn7596611). Save in directory databases/pcawg/histo.
* [PCAWG Kallisto TPM file from PCAWG](http://synapse.org/#!Synapse:syn7536587)
* [GTEx Kallisto TPM file from PCAWG](http://synapse.org/#!Synapse:syn7596599)
* Ensembl Biomart ID mapping. [Download from here](http://feb2014.archive.ensembl.org/biomart/martview)
* [isoform-specific-interaction file][#isoform]

<hr>

#### <a name="sampleTable"></a>
1. Determine for each cancer type all samples and create directory structure with GTEx tissue names and cancer types as subdirectories:

	```
	# Add PCAWG sample IDs
	zcat databases/pcawg/pcawg.rnaseq.extended.metadata.aliquot_id.V4.tsv.gz | tail -n +2 | perl -ne 'chomp($_); @a=split(/\t/); next if($a[34] eq "no" or $a[37] ne "Whitelist"); $a[57]=~s/\s//g; $a[75]=~s/\s//g; $h{$a[75]}->{$a[57]}.="$a[0]\tCA\n"; END{foreach $o(sort keys %h){mkdir $o if(!-d $o); foreach $k(sort keys %{$h{$o}}){mkdir "$o/$k" if(!-d "$o/$k"); open(O,">$o/$k/sampleTable.tsv"); print O "#sample\tcondition\n$h{$o}->{$k}"; close(O)}}}'
	
	# Add GTEx sample IDs
	for s in */*/sampleTable.tsv; do t=`echo $s | sed "s/\/.*//" | sed "s/Uteri/ Uteri/" | sed "s/Gland/ Gland/"`; for n in `zgrep "$t$" databases/pcawg/GTEX_v4.metadata.tsv.gz | cut -f2`; do echo -e "$n\tNO"; done >> `dirname $s`/sampleTable.tsv; done
	```

2. Extract for each cancer type the PCAWG expression values:
#### <a name="kallistoPcawg"></a>
	
	```
	for s in */*/sampleTable.tsv; do echo $s; less databases/pcawg/pcawg.rnaseq.transcript.expr.tpm.tsv.gz | perl -ane 'BEGIN{open(F,"'$s'"); while(<F>){chomp($_); @a=split(/\t/); $h{$a[0]}=$a[1]}} chomp($_); @a=split(/\t/); if($i++ == 0){@h=@a} print $a[0]; for($i=1; $i<@a; $i++){print "\t$a[$i]" if(exists $h{$h[$i]} and $h{$h[$i]} eq "CA")}; print "\n"' | gzip > `dirname $s`/pcawg.rnaseq.transcript.expr.tpm.tsv.gz & done

	# Add ENSG ID to each line	
	for s in */*/pcawg.rnaseq.transcript.expr.tpm.tsv.gz; do d=`dirname $s`; echo $d; less $s | perl -ane 'BEGIN{open(F,"zcat databases/ensembl/v75/ensg_ensp_enst_ense_geneName_v75.tsv.gz |"); while(<F>){@a=split(/\t/); $h{$a[2]}=$a[0]}}; $F[0]=~s/\..*//; print "$h{$F[0]}\t$_"' | gzip > $d/loschen.gz && /bin/mv $d/loschen.gz $d/pcawg.rnaseq.transcript.expr.tpm.tsv.gz & done
	```

3. Extract for each tissue type the GTEx expression values:
#### <a name="kallistoGtex"></a>

	```
	for s in */*/sampleTable.tsv; do echo $s; less databases/pcawg/GTEX_v4.pcawg.transcripts.tpm.tsv.gz | perl -ane 'BEGIN{open(F,"'$s'"); while(<F>){chomp($_); @a=split(/\t/); $h{$a[0]}=$a[1]}} chomp($_); @a=split(/\t/); if($i++ == 0){@h=@a} print $a[0]; for($i=1; $i<@a; $i++){print "\t$a[$i]" if(exists $h{$h[$i]} and $h{$h[$i]} eq "NO")}; print "\n"' | gzip > `dirname $s`/GTEX_v4.pcawg.transcripts.tpm.tsv.gz & done

	# Add ENSG ID to each line
	for s in */*/GTEX_v4.pcawg.transcripts.tpm.tsv.gz; do d=`dirname $s`; echo $d; less $s | perl -ane 'BEGIN{open(F,"zcat databases/ensembl/v75/ensg_ensp_enst_ense_geneName_v75.tsv.gz |"); while(<F>){@a=split(/\t/); $h{$a[2]}=$a[0]}}; $F[0]=~s/\..*//; print "$h{$F[0]}\t$_"' | gzip > $d/loschen.gz && /bin/mv $d/loschen.gz $d/GTEX_v4.pcawg.transcripts.tpm.tsv.gz & done
	```

4. Compute relative expression values for each transcript

	```
	# Determine relative expression values for each transcript
	scripts/isoformExpression.pl -meta $databases/pcawg/histo/pcawg.rnaseq.extended.metadata.aliquot_id.tsv.gz -expressFile databases/pcawg/expression/pcawg.rnaseq.transcript.expr.tpm.tsv.gz -stringSeq databases/string/v10/9606.protein.sequences.v10.fa.gz -seqFile databases/ensembl/v75/Homo_sapiens.GRCh37.75.pep.all.fa.gz | gzip >& isoformExpressionInCancer.tsv.gz
	
	# Split isoformExpressionInCancer.tsv.gz file for each cancer type
for s in */*/sampleTable.tsv; do echo $s; d=`dirname $s`; zcat isoformExpressionInCancer.tsv.gz | perl -ane 'BEGIN{open(F,"'$s'"); while(<F>){@a=split(/\t/); $h{$a[0]}=$a[1]}} @a=split(/\t/); if(exists $h{$F[1]}){print $_}' | gzip > $d/isoformExpressionInCancer.tsv.gz & done
	```


5. Identify for each cancer type cancer-specific MDT switches
#### <a name="mdtSwitch"></a>

	```
	for s in ../pcawgCAvsOnlyGtexFromCAtissue/*/*pcawg.rnaseq.transcript.expr.tpm.tsv.gz; do d=`dirname $s`; m=`echo $d | sed "s/.*tissue\///"`; mkdir -p $m; echo $m; scripts/mostDominantTranscriptSwitch.pl -pcawg $s -gtex $d/GTEX_v4.pcawg.transcripts.tpm.tsv.gz -ensg databases/ensembl/v75/ensg_ensp_enst_ense_geneName_v75.tsv.gz -canon databases/string/v10/canonEnsp_ensg_ensp_enst_geneName_v75.tsv.gz -iso interactionsInIsoforms_900_2.tsv.gz -minString 900 -minEnr 2 -maxQ 0.01 -minExp 1 -enst 1 | gzip > $m/interactionDisruptionInDominantTranscripts_min2.tsv.gz; done
	```

6. List for each MDT switch all disrupted interactions
#### <a name="mdtSwitchIntDisrupt"></a>

	```
	f=interactionDisruptionInDominantTranscripts_min2; for s in */*/$f.tsv.gz; do d=`dirname $s`; zgrep ^# $s | cut -f1-13,15-40 | tr -d '\n' > $d/$f"_int.tsv"; echo -e "\tPfam1\tDomain1\tRegion1\tPfam2\tDomain2\tENSPcanon2\tStringScore\tDIMAscore" >> $d/$f"_int.tsv"; zgrep -v ^# $s | perl -ane '@a=split(/\,/, $F[13]); foreach $a(@a){if($a=~/:/){$a=~s/\:/\t/g}else{$a="-\t-\t-\t-\t-\t-\t-\t-"}; for($i=0; $i<@F; $i++){print "$F[$i]\t" if($i!=13)}; print "$a\n"}' >> $d/$f"_int.tsv"; gzip -f $d/$f"_int.tsv" & done
	```

7. Add annotations to each disrupted interaction
#### <a name="annotate"></a>

	```
	for s in */*/interactionDisruptionInDominantTranscripts_min2_int.tsv.gz; do d=`dirname $s`; m=${s%.tsv.gz}; zcat $s | scripts/addAnnotation.pl -ensgFile databases/ensembl/v75/ensg_ensp_enst_ense_geneName_v75.tsv.gz -express $d/isoformExpressionInCancer.tsv.gz -density string_v10_density_900_3.tsv.gz -cosmic string_v10_cosmic_grch37_v89_neighbourhood_min900_shell3.tsv.gz -ensp1 1 -ensg1col 0 -ensp2 21 -v -step 12 -isoInt interactionsInIsoforms_900_2.tsv.gz > $m"_anno.tsv" & done
	```

8. Map mutations
#### <a name="mapMutations"></a>

	```
	# Extract mutations for each cancer type from main PCAWG MAF file
	cat databases/pcawg/October_2016_whitelist_2583.snv_mnv_indel.maf | perl -ane 'BEGIN{open(F, "zcat databases/pcawg/pcawg.rnaseq.extended.metadata.aliquot_id.V4.tsv.gz | tail -n +2|"); while(<F>){chomp($_); @a=split(/\t/); $h{$a[-1]}="$a[75]/$a[57]"}} chomp($_); if($i++==0){$header=$_}; @a=split(/\t/); if(exists $h{$a[12]}){$m{$h{$a[12]}}.="$_\n"} END{foreach $d(sort keys %m){open(O,">$d/October_2016_whitelist_2583.snv_mnv_indel.maf"); print O "$header\n$m{$d}\n"; close(O)}}'
	
	# Map mutations to each MDT interaction.
	for s in */*/interactionDisruptionInDominantTranscripts_min2_int_anno.tsv; do m=`basename $s | sed "s/.tsv//"`; d=`dirname $s`; scripts/mapMutation2dtu2.pl -dtu $s -regionFile databases/pcawg/all_gc19_pc.bed -mutFile $d/October_2016_whitelist_2583.snv_mnv_indel.maf -geneCol 25 -project `basename $d` -meta databases/pcawg/pcawg.rnaseq.extended.metadata.aliquot_id.V4.tsv.gz -anno databases/pcawg/all_annotation.tsv -eqtl databases/pcawg/all_eQTL.tsv -max 0.05 -enhancer databases/pcawg/map.enhancer.gene_max5.gz -v > $d/$m"_mutationsWithEnhancer.tsv"; done >& mapMutation_status.txt &
	```
	
9. Merge all results
### <a name="mergeAll"></a>

	```
	f=interactionDisruptionInDominantTranscripts_min2_int_anno_mutationsWithEnhancer; zcat Bladder/Bladder-TCC/$f.tsv.gz | head -1 | perl -ne '$_=~s/^#//; print "#Tissue\t$_"' > $f"_all.tsv"; for s in */*/$f.tsv.gz; do d=`dirname $s`; for t in `less $s | grep -v "^#"`; do echo -e "$d\t$t"; done; done >> $f"_all.tsv"; gzip $f"_all.tsv"
	```

10. Clean up mutational annotations. Allow a mutation to be either first as cds > ss > promCore > 5utr > 3utr > promDomain
### <a name="uniqMuts"></a>

	```
	less interactionDisruptionInDominantTranscripts_min2_int_anno_mutationsWithEnhancer_all.tsv.gz | scripts/mutationUniq.pl -mutCol 50 -eqtlCol 51 | gzip > interactionDisruptionInDominantTranscripts_min2_int_anno_mutationsWithEnhancerUniq_all.tsv.gz

