# 4DVOUS
VOUS annotation using 4D Genome data.
The app takes a text file with genomic coordinates as input (chr:start-end format) and matches them against 4D Genome database in order to find chromatin interactions outside the regions. A screening of the ASD database is also implemented in order to test associations between mapped genes and autism-associated genes.

The app is developed in R with a Python style option parsing from the command line

Usage: VOUS4D.R [options]

Options:
	
	-c COORDINATES, --coordinates=COORDINATES
		path to the text file with input coordinates in the format 'chr:start-stop'

	-m CHARACTER, --method=CHARACTER
		Choose overlapping method, options are: loose (match any region with overlaps >=1 BP), stringent (>50% overlap required). [default: loose]

	-d CHARACTER, --database=CHARACTER
		Path to the VOUS4D.db file, default assumes it is in the same folder

	-a CHARACTER, --asd=CHARACTER
		Perform match against ASD data? If TRUE an additional file ([outfile].ASD.txt) will be produced. [default: FALSE]

	-o CHARACTER, --outfile=CHARACTER
		Output file name '.txt' and '.xls' format are supported. [default: 4DG.out.txt]

	-h, --help
		Show this help message and exit

