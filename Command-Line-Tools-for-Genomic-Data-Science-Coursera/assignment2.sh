#Assignment 2 - Documented with output and explanation

echo "Q1. How many alignments does the set contain?"
#The number of lines in the alignment bam file is the no. of alignments in the bam file
echo "Answer: $(samtools view athal_wu_0_A.bam | wc -l)"

echo "Q2. How many alignments show the read’s mate unmapped?"
# The Reference position of the mate where its mapped (chromosome where the mate is mapped) (7th) in the bam file and '*' refers to an unmapped mate
echo "Answer: $(samtools view athal_wu_0_A.bam | cut -f7 |grep -c "*")"

echo "Q3. How many alignments contain a deletion (D)?"
#The 6th column in the bam file contains CIGAR value counting the number of D's representing deletions
echo "Answer: $(samtools view athal_wu_0_A.bam | cut -f6 | grep -c "D")"

echo "Q4. How many alignments show the read’s mate mapped to the same chromosome?"
#The Mate Ref pos in column 7 of bam file will have an "=" value if the mate is mapped to the same chromosome
echo "Answer :$(samtools view athal_wu_0_A.bam |cut -f7 |grep -c "=")"

echo "Q5. How many alignments are spliced?"
#The CIGAR value in the 6th column of the bam file will contain N if the alignment includes an intron in the sequence implicating splicing
echo "Answer $(samtools view athal_wu_0_A.bam |cut -f6 |grep -c "N")"

#Questions 6 to 10 require extracting a portion of the original bam file and same logic as above is used for the extracted file
samtools view -b athal_wu_0_A.bam Chr3:11777000-11794000 >athal_wu_0_A.extract.bam

echo "Q6. How many alignments does the set contain?"
echo "Answer: $(samtools view athal_wu_0_A.extract.bam | wc -l)"

echo "Q7. How many alignments show the read’s mate unmapped?"
echo "Answer $(samtools view athal_wu_0_A.extract.bam | cut -f7 | grep -c "*")"

echo "Q8. How many alignments contain a deletion (D)?"
echo "Answer $(samtools view athal_wu_0_A.extract.bam | cut -f6 | grep -c "D")"

echo "Q9. How many alignments show the read’s mate mapped to the same chromosome?"
echo "Answer $(samtools view athal_wu_0_A.extract.bam | cut -f7 |grep -c "=")"

echo "Q10. How many alignments are spliced?"
echo "Answer: $(samtools view athal_wu_0_A.extract.bam | cut -f6 |grep -c "N")"

#Question 11 to 15 use the original bam file

echo "Q11. How many sequences are in the genome file?"
#The header part of the sam file contains list of sequences and their corresponding lengths ( these lines begin with the @SQ ), hencing counting the occurances of @SQ 
echo "Answer: $(samtools view -H athal_wu_0_A.extract.bam |grep -c "@SQ")"

echo "Q12. What is the length of the first sequence in the genome file?"
#Extracting the 1st @SQ line length is given in as LN: followed by the length
echo "Answer: $(samtools view -H athal_wu_0_A.bam |grep "@SQ" | head -n 1)"

echo "Q13. What alignment tool was used?"
#The 1st line starting with PG has the informating regarding the tool used under ID:echo "Answer: $(samtools view -H athal_wu_0_A.bam |grep "PG" | head -n 1)"

echo "Q14. What is the read identifier (name) for the first alignment?"
#The 1st coloumn of the 1st line gives the read identifier  of the 1st alignment
echo "Answer: $(samtools view athal_wu_0_A.bam | head -n 1 |cut -f1)"

echo "Q15. What is the start position of this read’s mate on the genome? Give this as ‘chrom:pos’ if the read was mapped, or ‘*” if unmapped?"
#The 7th and 8th column of the alignment file contains the Ref position (chromosome and exact position) of where the read's mate mapps to in the genome, in case the 7th column has '=' it means the mate is paired to the same chromosome as the read so column 3 has this information
echo "Answer$(samtools view athal_wu_0_A.bam | head -n 1 |cut -f7-8)"

#Questions 16 to 20 use the comparision between the extracted regions and a feature gtf file 

#Covert the bam to bed file for intersect function using bamtobed
bedtools bamtobed -i athal_wu_0_A.extract.bam >athal_wu_0_A.extract.bed

echo "Q16. How many overlaps (each overlap is reported on one line) are reported?"
# count the number of lines post the bedtools intersect between the 2 files using -wo option
echo "Answer $(bedtools intersect -wo -a athal_wu_0_A.extract.bed -b athal_wu_0_A_annot.gtf | wc -l)"

echo "Q17. How many of these are 10 bases or longer?"
#the last column in the intersect file contains the no.of bases overlap , sort this numberically in the reverse order and obtaing the 1st line wher 9 ocures is the number(-1) of overlaps with length >10
echo "Answer: $(bedtools intersect -wo -a athal_wu_0_A.extract.bed -b athal_wu_0_A_annot.gtf | cut -f16 | sort -nr |grep -n "^9" |head -n 1)"

echo "Q18. How many alignments overlap the annotations?"
#The 1st 6 columns in the intersect file correspond to the alignments and getting the unique no. of rows in this is the number of alignments overlaping with the exons
echo "Answer: $(bedtools intersect -wo -a athal_wu_0_A.extract.bed -b athal_wu_0_A_annot.gtf | cut -f1-6 |sort -u |wc -l)"

echo "Q19. Conversely, how many exons have reads mapped to them?"
#The 5th to 15th column correspond to the lines in the gtf (exons) that overlap with the alignment file. so sorting this uniquely gives the number of exons having overlap
echo "Answer: $(bedtools intersect -wo -a athal_wu_0_A.extract.bed -b athal_wu_0_A_annot.gtf | cut -f7-15 |sort -u |wc -l)"

echo "Q20. If you were to convert the transcript annotations in the file “athal_wu_0_A_annot.gtf” into BED format, how many BED records would be generated?"
# We have to count the number of transcripts in the gtf file present in the 9th column of the gtf file.
echo "Answer: $(cut -f9 athal_wu_0_A_annot.gtf |cut -d " " -f1,2 |sort -u |wc -l)"




















