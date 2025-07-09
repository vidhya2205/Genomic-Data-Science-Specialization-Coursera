#Course Summary
This is summary of the key learnings noted based on the course [“Command Line Tools for Genomic Data Science”]{https://www.coursera.org/learn/genomic-tools?msockid=263c03b77fe161fb36cf104d7e47606b) 
##**Week 1:**
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
Assignment1 – My original solutions to the questions and the output is shared in [bash script](https://github.com/vidhya2205/Command-Line-Tools-for-Genomic-Data-Science-Coursera/blob/main/assignment1.output.txt) and [output](https://github.com/vidhya2205/Command-Line-Tools-for-Genomic-Data-Science-Coursera/blob/main/assignment1.sh)
