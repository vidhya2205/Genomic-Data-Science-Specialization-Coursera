#!/bin/bash

# Assignment 1 - Documented with output and explanation

echo "Q1. Number of chromosomes in the genome:"
#Every chromosome's start is indicated by '>' in the .genome file so counting the occurances of '>' will provide the number of chromosomes
echo "Answer: $(grep '>' apple.genome | wc -l)"

echo "Q2. Number of genes:"
#The 1st column of the .genes file has the list of genes removing the duplicates will fetch the number of genes
echo "Answer: $(cut -f1 apple.genes | sort | uniq -c | wc -l)"

echo "Q3. Number of transcript variants:"
#The 2nd column of the .genes file has the list of gene variants and hence that will give the no. of transcript variants
echo "Answer: $(cut -f2 apple.genes | sort | uniq -c | wc -l)"

echo "Q4. Number of genes with a single splice variant:"
#uniq -c give the count of occurance of every gene hence grep of '1' will provide genes with single variants
echo "Answer: $(cut -f1 apple.genes | sort | uniq -c | grep ' 1 ' | wc -l)"

echo "Q5. Number of genes with 2 or more splice variants:"
#uniq -c give the count of occurance of every gene hence grep of -v (not including)'1' will provide genes with >1 variants
echo "Answer: $(cut -f1 apple.genes | sort | uniq -c | grep -v ' 1 ' | wc -l)"

echo "Q6. Number of genes on the '+' strand:"
#the 4th column in the .genes file has the strand information and hence extracting that along with 1st column(genes) and removing duplicates and then counting the number of '+' 
echo "Answer: $(cut -f1,4 apple.genes | sort | uniq -c | grep -c '+')"

echo "Q7. Number of genes on the '-' strand:"
#the 4th column in the .genes file has the strand information and hence extracting that along with 1st column(genes) and removing duplicates and then counting the number of '-'
echo "Answer: $(cut -f1,4 apple.genes | sort | uniq -c |grep -c "-")"

echo "Q8. Number of genes on chromosome chr1:"
#the 3rd column in the .genes file has the chromosome information and hence extracting that along with 1st column(genes) and removing duplicates and then counting the number of 'chr1'
echo "Answer: $(cut -f1,3 apple.genes | sort -u | grep 'chr1' | wc -l)"

echo "Q9. Number of genes on chromosome chr2:"
#the 3rd column in the .genes file has the chromosome information and hence extracting that along with 1st column(genes) and removing duplicates and then counting the number of 'chr2'
echo "Answer: $(cut -f1,3 apple.genes | sort -u | grep 'chr2' | wc -l)"

echo "Q10. Number of genes on chromosome chr3:"
#the 3rd column in the .genes file has the chromosome information and hence extracting that along with 1st column(genes) and removing duplicates and then counting the number of 'chr3'
echo "Answer: $(cut -f1,3 apple.genes | sort -u | grep 'chr3' | wc -l)"

echo "Q11. Number of transcripts on chromosome chr1:"
#the 3rd column in the .genes file has the chromosome information and hence extracting that along with 2nd column(transcript variants) and removing duplicates and then counting the number of 'chr1'
echo "Answer: $(cut -f2,3 apple.genes | sort -u | grep 'chr1' | wc -l)"

echo "Q12. Number of transcripts on chromosome chr2:"
#the 3rd column in the .genes file has the chromosome information and hence extracting that along with 2nd column(transcript variants) and removing duplicates and then counting the number of 'chr2'
echo "Answer: $(cut -f2,3 apple.genes | sort -u | grep 'chr2' | wc -l)"

echo "Q13. Number of transcripts on chromosome chr3:"
#the 3rd column in the .genes file has the chromosome information and hence extracting that along with 2nd column(transcript variants) and removing duplicates and then counting the number of 'chr3'
echo "Answer: $(cut -f2,3 apple.genes | sort -u | grep 'chr3' | wc -l)"

echo "Q14. Number of genes in common between condition A and B"
#extract only the genes and store the deduplicated genes in each condition in a file
cut -f1 apple.conditionA |sort -u > genesconditionA
cut -f1 apple.conditionB |sort -u > genesconditionB
cut -f1 apple.conditionC |sort -u > genesconditionC
#comm compares the 2 files and present 3 columns (column1 - unique to A,column 2 -unique to B and column 3- common to all 3)
echo "Answer (method 1): $(comm -1 -2 genesconditionA genesconditionB | wc -l)"
#cat concatenates the files and the sorting and using uniq -c and get the count of each gene if its 2 then its present in both the files
echo "Answer (method 2): $(cat genesconditionA genesconditionB | sort | uniq -c | grep '2 ' | wc -l)"

echo "Q15. Genes exclusive to condition A:"
echo "Answer: $(comm -3 -2 genesconditionA genesconditionB | wc -l)"

echo "Q16. Genes exclusive to condition B:"
echo "Answer: $(comm -3 -1 genesconditionA genesconditionB | wc -l)"

echo "Q17. Genes common to all three conditions:"
echo "Answer: $(cat genesconditionA genesconditionB genesconditionC | sort | uniq -c | grep '3 ' | wc -l)"

