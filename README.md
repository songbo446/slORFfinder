####################################################################################

slORFfinder: a tool to detect trans-ORFs resulted from trans-splicing of spliced leader sequences

Bo Song (songbo446@yeah.net)

AGIS

####################################################################################

0. DEPENDENCIES

0.1 R packages:
	
	multitaper
	pheatmap
	
0.2 Python packages:

	math
	pickle
	scipy
	multiprocessing
	functools
	
0.3 Other tools:

	bedtools
	samtools

1. INTRODUCTION

	slORFfinder is coded in Python3 and is developed for Linux users. 
	Make sure Python3 has been installed.
	
2. INSTALLATION

	Download slORFfinder.tar.gz and extract the files in YOURPATH
	> tar xzf slORFfinder.tar.gz
	> cd slORFfinder
	> chmod +x slORFfinder
	
	Now you can excute directly by typing "YOURPATH/slORFfinder" in command line
	Alternatively, you can also run the program by typing "python3 YOURPATH/slORFfinder".
	If it works, th help message will appear as follows:

        	slORFfinder [options]* --genome <reference> --gtf <annotation> --SL <spliced-leader sequence> --rna <alignment of RNA-Seq reads> --ribo <alignment of Ribo-Seq reads>

        	<reference>                     Sequence of reference genome (fasta)
        	<annotation>                    Annotation of reference genome (gtf)
        	<spliced-leader sequence>       Spliced-leader sequence in captical letters (string) (Multiple SLs separated by comma; Ambiguous nucleotide is not allowed)
        	<alignment of RNA-Seq reads>    Alignment of RNA-Seq reads (bam)
        	<alignment of Ribo-Seq reads>   Alignment of RPFs (bam)

    	Options:

        	--start         	Start codons (spearate by comma) [default: AUG,UUG,GUG,CUG]
        	--pcov          	Minimum RPF coverage of Psites (0-1) [default: 0]
        	--nCores        	Number of multiprocessors [default: 5]
        	--transOnly     	Output all predicted ORFs (0) or only slORF (1) [default: 1]
        	--outdir        	Output directory [default: ./slORF]
        	--prefix        	Prefix of output files [default: ribosl]
        	--Rscript       	Path to Rscript
        	--bedtools      	Path to bedtools
        	--samtools      	Path to samtools
        	--slseed        	The number of soft clipped reads [default: 8]
        	--filetype      	bam/sam file type [default: bam]
	

3. USAGE

	slORFfinder requires the inputs of the reference genome sequence (in fasta format), the genome annotation (in GTF format), RNA-Seq alignments (in bam/sam format) and Ribo-Seq alignments (in bam/sam format) and a sequence or sequences of SL(s). 

	There are 11 options
	
	--start
	
	slORFfinder use cognate start codons (NUG) by default. It can also be used to predict canonical ORFs initiated by AUG when only 'AUG' is inputted.
	
	--pcov
	
	The minimum coverage of RPF at P-sites for a canidate ORF. This value should be set between 0-1. By default, slORFfinder skip if there is no RPF-supported P-site in a candiate ORF.
	
	--nCores
	
	The number of multiprocesssors. The identification of ORFs could be speed up by using more processors simultaneously. Five tasks will be processed parellely by default. If '--transOnly' is truned on (by default), it will not cost much resources and the defalt parameter is recommanded.
	
	--transOnly
	
	slORFfinder can also be used to predict ORFs other than slORFs. If the users want to predict slORFs and other ORFs such as short ORFs, this option can be turned off. Please note that if '--transOnly' is turned off, the computing tasks increased dramatically, it's better to invoke more processors by seeting '--nCores' to speed up the process.
	
	--outdir
	
	The output directory of results.
	
	--prefix
	
	The prefix of output files.
	
	--Rscript
	
	Some processes of slORFfinder require the pre-installation of some R pacakges. If the commands of Rscript cannot be aumatically invoked by the tool, the usrers can give a path of Rscript here.
	
	--bedtools
	
	slORFfinder uses bedtools to manupilate sam/bam files in some steps, if the commands of bedtools cannot be found automatically, the path to it can be given here.
	
	--samtools
	
	slORFfinder uses samtools to manupilate sam/bam files in some steps, if the commands of samtools cannot be found automatically, the path to it can be given here.
	
	--slseed
	
	slORFfinder searches for the relicts (8 bp at the tail by default) of SL sequences in reads. The size of the target relic sequence can be set here.
	
	--filetype
	
	slORFfinder takes the alignment files in bam format by default, if the files are in sam format, the users can change this option to 'sam'.
		
4. RESULTS

	slORFfinder outputs eight files and a directory "Plot" in which metagene and multitaper plots are
	included if the whole pipeline is completed.

	4.1 prefix_OffSet.txt
	
	It lists the offsets for RPFs in different sizes
	Column 1: Size of RPF
	Column 2: Offset
	Column 3: Proportion

	4.2 prefix_multitaper.Freq.txt

	The results of multitaper tests of RPFs in various sizes
	Column 1: Size of RPF
	Column 2: Minus Log10 of P-values
	Column 3: The peaked frequency
	
	4.3 DataSumForPrediction.txt
	
	It includes the RPFs and their amounts that were selected for ORF prediction, as well as
	the entropy used in this prediction.
	
	4.4 prefix_summary.txt
	In this file, the numbers of all the types of predicted ORFs are listed.

	4.5 prefix_orf.fa
	
	The DNA sequences of predicted ORFs in fasta format. The ID line is formmated as:
	
	'>'ID_of_ORF	Chromosome	Strand	Name_of_transcript		Combined_P-value	P1:P2:P3:P4
	
	The ID_of_ORF is consisted of the positional information as formmated as:
	
	Name_of_transcript:Chromsome:Strand:Coordinate_of_start:Coordinate_of_stop

	4.6 prefix_orf_aa.fa

	The amino acid sequences translated from prefix_orf.fa

	4.7 prefix_ORF.gff

	The annotation of predicted ORFs in .gff format

	4.8 prefix_orf_table.txt

	The ID, class, DNA and AA sequences of all the predicted ORFs

	4.9 Plot
	A directory including metagene and multitaper plots.

	Note:
	Classification of predicted ORFs: 
	(a) annotated ORF, the ORFs identical with the annotated ones
	(b) truncated ORF, the ORFs with the same start or stop codon but shorter than the annotated ones
	(c) extended ORF, the ORFs with the same start or stop codon but longer than the annotated ones
	(d) uORF, upstream ORF, ORFs located in 5'UTR
	(e) ouORF, overlapped uORF, ORFs located in 5'UTR and overlapped with the annotated start codon
	(f) dORF, downstream ORF, ORFs located in 3'UTR
	(g) odORF, overlapped dORF, ORFs located in 3'UTR and overlapped with the annotated stop codon
	(h) ncsORF, ORFs located in non-coding RNAs, ORFs predicted from a gene without any annotated CDS will also be classified into ncsORF
	(i) teORF, ORFs located in transponsable elements
	(j) pORF, ORFs on pseudo genes.
	(k) slORF, ORFs resulted from trans-splicing of spliced leader sequences
	
