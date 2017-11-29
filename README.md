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

# An example
```{shell}
$ Rscript ./VOUS4D.R --coordinates dummy.coordinates.txt --asd TRUE --out test.dummy.txt
MacBookNando:vous4d febo$ Rscript ./VOUS4D.R --coordinates dummy.coordinates.txt --asd TRUE --out test.dummy.txt

=======================================================================
||                                                                   ||
||  Script developed by: Ferdinando Bonfiglio (nandobonf@gmail.com)  ||
||                                                                   ||
=======================================================================

### Loading required packages... DONE!

### input coordinates:

|input.chr | input.start| input.end|
|:---------|-----------:|---------:|
|chr1      |    31584329|  31656511|
|chr16     |      100000|    110000|
|chr22     |      100000|    210000|

### Reading 4DGenome and ASD data... DONE!

### 25 interactions found in 4Genome
### Method: LOOSE overlap in H1ESC and 4 derived cells
### Method: removed 2 matches with both interactors overlapping input coordinates

### 7 unique genes mapped
### Method: removed matches without gene symbol annotation

### 7 input genes found in Genome-wide predictions of autism-associated genes (ASD) 

### Top genes (ranked by q-value):
-------------------------------------------------------------------------------
  genes         cells            methods        score    probability   q-value 
--------- ----------------- ----------------- --------- ------------- ---------
 SERINC2        H1ESC            IM-PET        -0.9271     0.3325      0.01268 

 RHBDF1       H1ESC, H1        5C, IM-PET,     -0.9238     0.3336      0.02095 
               derived       IM-PET, IM-PET,                                   
             mesenchymal     IM-PET, IM-PET,                                   
            stem cell, H1    IM-PET, IM-PET                                    
               derived                                                         
             mesendoderm                                                       
              cell, H1                                                         
               derived                                                         
            trophoblast,                                                       
              H1ESC, H1                                                        
           derived neural                                                      
           progenitors, H1                                                     
               derived                                                         
           trophoblast, H1                                                     
               derived                                                         
             mesenchymal                                                       
              stem cell                                                        

 POLR3K      H1 derived      IM-PET, IM-PET,   -0.9428     0.3297      0.1103  
             mesenchymal     IM-PET, IM-PET                                    
            stem cell, H1                                                      
               derived                                                         
           trophoblast, H1                                                     
               derived                                                         
             mesendoderm                                                       
             cell, H1ESC                                                       
-------------------------------------------------------------------------------

### ASD output file written to: test.dummy.ASD.txt 

### 4D Genome output file written to: test.dummy.txt 
### Analysis completed in 5.09094 seconds


```
