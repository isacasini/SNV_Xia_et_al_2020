ShortReadMe for “Base Editing Opportunities via PAM”
Prepared by Isabella Casini
Last Updated: 24 March 2020


Introduction:
The following file is provided as a brief ReadMe file, for someone who would like to edit and use this program. 
Everything was coded in Python 3.6 using Pycharm (2018.1) with a Conda (Anaconda) environment, “Base Editing”, 
the details for which are provided in a separate .yml file (Base_Editing_Env.yml). The final directory contains this file, 
a more detailed .docx version of this file (with a more detailed description of what program does and each function - which can also be accessed through the help() 
function in Python), the zipped PyCharm files, the input and output files that were used and produced,and the Conda environment file. 

General Notes:
1. Suffix “T” refers to Top strand
2. Suffix “B” refers to Bottom strand
3. Suffix “df” refers to the type of the variable (pandas DataFrame)
4. Suffix “dict” refers to the type of the variable (dictionary)
5. Suffix “list” refers to the type of the variable (list)
6. “PAM”, “AA”, “CDS” are always kept all capitalized

User Required Input: 
The main file has places where the user input is required, these are labeled “USER_INPUT_REQUIRED”. The instances are summarized below. 
1.	(STEP 1) Provide the path to the top strand sequence file (txt)
a.	# USER_INPUT_REQUIRED - change the path in TopSequenceFile
2.	(STEP 1) Provide a path out for the complement sequence file (filecomp)
a.	# USER_INPUT_REQUIRED - change the path in filecomp
3.	(STEP 2) To change the PAM sequence, go into the position_pam.py function, findPAM and change (both instances of) “([atc]gg)” to the letter combination of interest.
4.	(STEP 2) To change the length and position of the editable area (default is -19 to -11) in the protospacer change (both instances of) the number of base pairs (bp#) and the number in the “[start-#]”.
5.	(STEP 19) Provide the path out for the final excel sheet with all the calculations.
a.	# USER_INPUT_REQUIRED - change the path in writePathout

External Files Required:
1.	DNA top strand sequence file: 
a,	A text file with the nucleotide sequence
b.	The reader (readtxt) will skip a line that starts with “>”
2.	File with the CDS (gene) information
a.	A tab delimited text/csv file
b.	Should (bold are optional) to contain the following columns/headings “CDS, CDS length, direction, start, stop, annotation”
i.	E.g. “locus	length bp	direction	start	stop	NCBI Name	NCBI Locus	NCBI Annotation”
c.	Direction must be given as “=>” as forward and “<=” for backward.
