# Course Summary
This is summary of the key learnings noted based on the course [“Command Line Tools for Genomic Data Science”](https://www.coursera.org/learn/genomic-tools?msockid=263c03b77fe161fb36cf104d7e47606b) 
## **Week 1:**
List of commands:
1.	Basic commands
Just line in every operating system we can communicate with the Unix OS using command line tools in the terminal-
-	cd : change directory, we can provide the absolute path starting from the home directory or relative path from the current working directory 
-	pwd – present working directory
-	ls – list the content of the current directory (options:
o	 -l :provides a more detailed info about each file/ directory in the current directory
o	-lt : provides the list in the order of which the files were created.
**Wild cards** are useful for finding a list of files that meets our needs. Like 
-	The ‘*’ is used to represent any combination of letters and numbers, so any string. 
-	The ‘?’ now stands for just one character.
-	The [ ] is used to represent a range of characters that are enclosed within the square brackets
-	The { } represents to check for all files that have the characters enclosed we can check for files with have 2 or more names in out directory.
2.	Content creation and remove and move content
-	mkdir : create a directory
-	cp : to copy files and place them in a specific location (creates a duplicate file)
-	mv : to move files to a different location ( replaces the original file from the initial location to the new location)
-	rm : remove files 
-	rmdir : Removing a directory ( make sure the directory is empty before deleting the directory)
3.	Accessing content
-	more : to display content page by page, only forward navigation
-	less : display content, forward and reverse navigation, search using /(forward), ?(backward) and more options
-	head : Displays the top 10 (default, if we want we can specify by -n(number)) lines of the file
-	tail : Displays the last 10 (default, if we want we can specify by -n(number)) lines of the file
-	cat : concatenates the contents of the files and gives it as an output.
4.	Redirecting content
Usually there is a standard input and standard output which is the command line input and terminal out respectively. However we can redirect this using ‘<’ (input), ’>’ (output) and | (piping)
5.	Querying content
-	sort : Sorts the content of the file in alphabetical order ( -r for reverse alphabetical order). If we want to sort by a particular column then we have to specify column number using -k followed by the “column number” and  “n” (to specify that numerically and not ascii). -u option is used to obtain unique values (Sorts the lines first, then removes all duplicates (globally)).
-	cut : to access 1 specific column in the file we use -f”column_number”. It keeps the default delimiter as a tab if we use a different one we can specify that using -d “delimiter”
-	uniq : Remove adjacent duplicates only.  Option ‘-c’ gives is the count of the number of times the word or line appears in a certain context. 
-	grep : searches lines of text for matches to a given pattern (usually a string or regular expression) and prints those lines. ‘-n’ option gives the line number in the file the pattern occurs. ‘-v’ provides the lines where the pattern is not found.
6.	Comparing content
-	diff : Show line-by-line differences between unsorted files
-	comm : Show which lines are common or unique (column wise) between sorted files. Selected columns can be displayed by ‘-“column_number”’
7.	Archiving content
-	gzip :Compresses a file using the .gz format.
-	gunzip : Decompresses a .gz file.
-	bzip2 :Compresses a file using Burrows–Wheeler algorithm, creates .bz2 files.
-	bunzip2 : Decompresses a .bz2 file.
-	tar : Archives multiple files into one file (optionally compresses with -z or -j). For archiving -cvf and dearchiving/extracting -xvf
-	zcat : Displays the contents of a .gz file without decompressing it to disk.

**Assignment1** – My original solutions to the questions and the output is shared in [bash script](assignment1.sh) and [output](assignment1.output.txt)

To execute the bash script - 

```bash
chmod +x assignment1.sh
./assignment1.sh
```
### Week 2 :  Sequences & Genomic Features. 
1.	Introduction to Molecular biology concepts and gene expression.
-	Genome  is the complete genetic material present within a cell. It is organized into 3 dimensional chromosomes (22 pairs  of autosomal chromosomes and 2 sex chromosomes) by complex folding.  Important part of this genome are genes that is the protein coding section. 
-	The constituting material, DNA is double stranded and hence genes can be on either the forward strand or reverse strand. Looking closely at the genes they are further divided into exons introns and many other elements.
-	A process called gene expression allows the conversion of the information from the DNA to be converted into proteins that run the functions within the cell. The 1st step is transcription process that converts the DNA into mRNA involving multiple steps in-between. DNA is transcribed into pre RNA that contains a 5’ cap and 3’ tail which is removed in the next step. Mainly the pre-RNA also contains introns that are to be spliced and based on the introns splicing site there are different mRNA’s made and hence different proteins. 
-	This mRNA is then translated into proteins.
-	Each gene can have multiple copies of mRNA in the nucleus(where transcription occurs) or cytoplasm and multiple copies of proteins in the cytoplasm (where translation occurs)
-	The role of a computational biologist is to find out the role of these genes and proteins their sequences and quantify them to find their expression levels.   
2.	Sequence Representation 
-	Genomic sequences are a set of strings formed by the combinations of the 4 nucleotide bases adenosine, cytosine, guanine and thymine (DNA) or uracil (RNA) represented as A, C, G and T/U respectively. Protein sequences are represented by the corresponding letters used to indicate the 20 amino acids.  
-	These sequences are generated by 1st extracting the genomic data (DNA or RNA) from the cell and fragmented into roughly a particular size, amplified and then sequenced which is the wet lab part of the process. 
-	The second part is the assembly of these fragments to get the entire sequence which is the bioinformatics part.  The reading of the sequence can be in one direction (single end sequencing) or from both sides (paired end sequencing). 
-	The assembly usually is done by finding overlaps between the different fragments (reads) and an overlap graph is generated. This graph is then transversed to get the complete sequence. 
-	Upto 2008, the sanger sequencing was used to generate the sequences of the reads. This method produced reads of longer length (500-600 bp) and was expensive, medium to high through put data generated.  It required amplification into cloning vector 
-	Later came the Next generation sequencing methods (NGS), which generated shorter reads (40-450 bps) and less expensive. It is also very fast compared to the previous method generating high throughput data. 
-	Sequences that are generated are represented in the ‘fasta’ format. Each read that is sequenced will be represented by 4 lines. The first line starts with a ‘>’ called the header line contains information (identifier) regarding the read/sequence.  Followed by multiple lines or a single line of the sequence (ACTG and occasion N (not identified base). There the 3rd line is a plus alone or followed by the sequence information(in the 1st line). The 4th line contains information relating the sequence quality called phred quality. This value gives the confidence of the base call being correct. How confidant is the sequencer that that base is present in that specific location. This is represented by  ASCII values corresponding to 33-126 and values(0 to 0.93). 
`Q =-10log<sup>10</sup>pB  
Where pB is the probability of base B occurring at that position and Q is the quality score at that position in the sequence.
-	The maximum quality score is around 40 and scores less that 20 or 25 based on the application are considered to be low quality reads. 
3.	Genomic Features:
-	Genomic annotation is to represent the specific location of a genomic feature in the genome associated with its biological information (exon, intron, non coding, protein binding sites, transcription start site, etc..)
-	When we are annotating a genefor instance, we need to add information regarding where are the exons present, how many exons are there, which strand it is present in, which chromosome, translational start site etc.
-	These annotations can be made in different formats of files like BED (browser extensible data), GTF (genomic transfer format) and GFF (genomic feature file).
-	BED file in the basic format requires only 3 columns of information, column 1 corresponding to chromosome number, column 2 and 3 represent the start and end of the feature. This file is zero based counting ( meaning the numbering starts from the space before each base (in the genomic sequence). 
-	The extended version of the file contains more columns for name of the feature, score (confidence level on the feature), strand (+ forward, - reverse) on which the feature is present, numerical values to indicate thick start and end and finally rgb_color. The thick_start an thick_end are to tell where exactly should be highlight the genome to represent the genomic feature and rgb_color to indicate the color in which the feature has to be represented. 
-	To represent features belonging to a single gene together, 2 additional columns represent the number of individual features(exons for example) and the relative positions (zero based from the start of the gene) is added. 
-	GTF format contains a few extra columns, the 1st column corresponds to chromosome number, 2nd is the source (program), 3rd is the type of feature (exon etc.), 4th and 5th is the beginning and end of the feature along the genomic access, 6th is score, 7th is strand, 8th is the frame (reading frame), 9th is a composite column with 2 different types of fields (1- gene identification, 2- transcript identification). Additional columns indicating information related to abundance, exon number can also be included. 
-	The features belonging to a single gene is grouped. The coordinate indication is based on 1 system (counts the bases from 1)
-	GFF format is evolved from GTF format, it has a header line which preceeds all the columns that specifies the version. Column 1 represents the chromosome number again, 2 represents the source(program-cufflinks, genmark, C4 etc..), 3 is the type of feature, 4 and 5 are the numerical start and end of the feature, 6 is score, 7 is strand, 8 is the frame, 9 contains information regarding the feature(ID).
4.	Alignment
-	Alignment is the mapping of a sequence read to a reference genome, identifying similarities and differences. A perfect match is not possible due to polymorphisms and sequencing errors. 
-	Two types of alignments are observed, continuous alignment where the read maps completely within the exon. Spliced alignment occurs when the read spans along the exon-intron border in case of DNA and exon-exon in case of mRNA
-	Properly paired alignment occurs when the alignment follows a specific distance and orientation between the mapping of the pairs (concordant mapping), non-concordant mapping is when these rules are violated.
-	SAM format :
o	Starts with a header file that includes information about if the file is sorted, what type of software is used to generate the file, what are the command-line options used in the generation. 
o	Each line followed by the header represents a single read. It contains fields such as QNAME- query name(header from FASTA file), FLAG, RNAME- chromosome number, Pos- start of the read in the genome(reference), MAPQ- Mapping quality score, CIGAR, RNEXT – reference name of the mate (for paired read) /next read, PNEXT- position of the mate (for paired read) /next read, [in case of unmapped mate it will be a *], TLEN- distance between the mates, SEQ – actual sequence, QUAL- ASCII encoded  base quality scores.  
-	FLAG – binary representation that conveys information about sequence alignment such as if its paired or single end sequencing, is the read mapped properly, segment is mapped or unmapped, is the mate mapped,  if the reverse complement of the read is mapped, if the mate’s reverse complement maps to the genome, if this is the 1st seqment or last segment in a pair of reads, secondary or primary alignment based on quality of alignment, if it passes the quality check,  if the read is a PCR duplicate, if it’s a supplementary alignment. 
-	The CIGAR field represents alignment in a compressed format, using a sequence of numbers and letters to denote edit operations like matches, insertions, and deletions. Different symbols used – M: Match (nucleotide matches the reference), I: Insertion (nucleotide added to the read), D: Deletion (nucleotide removed from the reference), N: Represents an intron (a gap in the reference), S: Soft clipping (bases not aligned but included in the sequence), H : Hard clipping (bases not included in the sequence), P: Padding (used for alignment purposes).
5.	Retrieving sequences and features
-	Several online databases can be used to retrieve sequence data.
-	To retrieve a specific gene sequence, such as the IL-2 gene in humans, a string-based search and filter results to find the desired mRNA sequence can be done in the NCBI GenBank database.
-	To retrieve fasta files from sra(sequence read archive) **command line tools** are handy-
-	wget : followed by the link to the project or specific file we are interested in downloads the file in the respective directory.
```bash
wget [URL]
```
-	fastq-dump:  converts the SRA files to FASTQ format that is used further. This is found in the sratoolkit. 
```bash
#for installation
sudo apt update
sudo apt install sra-toolkit
#First time users have to configure it as well
vdb-config –interactive
#Use the arrow keys to go through options. You can set the download directory or use default settings and save.
#for running
fastq-dump [SRA_file]
```
For paired end reads-
```bash
fastq-dump --split-3 [SRA_file]
```
-	nohup : is to allow process to continue running in the background even if the terminal is closed. The & at the end puts the command in the background.
```bash
nohup fastq-dump [SRA_file] &
```
-	We can also retrieve annotations files from UCSC, we can select the genomic features we are looking for and the format as well (GTF, BED, etc..)
6.	Samtools command
-	SAMtools is a software package designed for manipulating SAM and BAM formatted files in genomic data analysis. It provides various command line options for tasks such as indexing, editing headers, and removing PCR duplicates. Operations include sorting, merging, and converting BAM files to FASTQ format, as well as extracting information like depth and alignment statistics.
-	Key Commands and Options
o	**flagstat**: Provides statistics on the number of alignments in a BAM file.
```bash
samtools flagstat “file_name.bam”
```
o	**sort**: Sorts a BAM file for efficient indexing and visualization.
```bash
samtools sort “input_file_name.bam” -o “output_file_name.sorted.bam”
```
-o:  option used to specify the output file name.
o	**index**:  Creates an index for a sorted BAM file.
```bash
samtools index example.sorted.bam
```
o	**merge**: Merges multiple BAM files into one.
```bash
samtools merge “output.bam” “input1.bam” “input2.bam"
o	**view**: This command allows us to view alignments in BAM files as BAM files, and convert them to SAM format and vice versa, view alignments with a specific range.
•	Viewing Alignments
```bash
samtools view “file_name.bam”
```
Options:
-h: Include header information in the output.
-H: Include only the header information
•	Conversion Between Formats
BAM to SAM:
```bash
samtools view “input_file.bam” > “output_file.sam”
```
SAM to BAM:
```bash
samtools view -b -T hg38c “input_file.sam” > “output_file.sam.bam”
```
Options:
-b: Output in binary (BAM) format.
-T: Specify the reference genome file.
•	Extracting Alignments Within a Range
Using Coordinates:
```bash
samtools view “input_file.bam” [range]
```
Using a BED File:
```bash
samtools view -L “input_file_with coordinates.bed” “input_file.bam”
```
Options:
-L: Include only reads overlapping the specified regions in the BED file.
-	Commands can be run in the background using nohup to prevent interruptions.
-	Example for background execution: 
```bash
nohup samtools sort example.bam -o example.sorted.bam 
```
7.	Bamtools: 
It contains a set of tools for performing genome arithmetic and manipulating genomic intervals, format conversion, and sequence manipulation. We can extract sequences within specific genomic intervals and determine overlaps between gene annotations.
**Recap** - The GTF format represents each exon as a feature, while the BED format compacts multiple exons into single lines for each gene.
•	Finding Overlaps (Interval)
Use Case: Identify overlapping genomic features, such as exons and regulatory elements.The intersect command is used to find overlapping intervals between two sets of genomic features. Options for the intersect command allow users to customize the output, such as reporting overlaps and the extent of overlaps.
```bash
bedtools intersect -a “fileA” -b “fileB” [options]
```
Options – 
o	-wa: Write the original entry in A for each overlap.
o	-wb: Write the original entry in B for each overlap.
o	-wo: Write original A and B entries that overlap, concatenated along one line, including the number of bases that overlap.
o	-wao: Write all entries in A, with overlaps from B, including the number of bases that overlap.
o	-f <fraction>: Specify the minimum fraction of overlap required for features to be considered overlapping (default is 1 base pair).
o	-r: Requires that the fraction of overlap be applied equally to both A and B.
o	-split: Treats multi-block BED entries as separate annotations for intersection.
o	-c: Count the number of overlaps.
•	Extracting Sequences:
Extract sequences from a reference genome based on specific intervals.
```bash
bedtools getfasta -fi reference.fa -bed intervals.bed -fo output.fa
```
Options
o	split option ensures that the output reflects the correct lengths of exonic blocks, providing a more accurate genomic sequence representation ignoring the introns.
•	Counting Features
Count the number of features overlapping with a specific set of intervals.
```bash
bedtools intersect -a features.bed -b regions.bed -c
```
•	Format Conversion
Convert between different file formats, such as BAM to BED.
```bash
bedtools bamtobed -i input.bam > output.bed

```
Options:
o	-split : will split the alignments into a number of blocks corresponding to number of exons covered in that alignment.

Converting bed to bam file
```bash
bedtools bamtobed -i “input.bam/gff/vcf” -g “genomefile_specific_format”
```
Options:
o	-bed12 : represents the bed file along with the exon intron positions(separations-CIGAR)
•	Merging Overlapping Features
Combine overlapping intervals into a single feature.
```bash
bedtools merge -i features.bed
```
•	Subtracting Intervals
Remove specific intervals from a set of features.
```bash
bedtools subtract -a features.bed -b regions.bed
```
•	Clustering Intervals
Group nearby features into clusters.
```bash
bedtools cluster -i features.bed
```
**Assignment2** – My original solutions to the questions and the output is shared in [bash script](assignment2.sh) and [output](assignment2.output)

To execute the bash script - 

```bash
chmod +x assignment2.sh
./assignment2.sh
```
