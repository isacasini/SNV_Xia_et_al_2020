# SNV_Xia_et_al_2020
There are two files, ReadMe_brief.txt and readme.docx that provide brief and detailed instructions, respectively.

"The algorithm [aims] to identify all editable sites at genome-scale to identify possible single-nucleotide variations and mutations at translation level." The files included here are for Clostridium ljundahlii DSM13528. This program is used in Xia, P., Casini, I., et al. 2020. Reprogramming acetogenic bacteria with CRISPR-targeted base editing via deamination. (submitted)

The directory, "Input Files" are the files that the program requires as an input. These should be modified based on the organism. 
The directory, "Output Files" are the outputs of the program, including the compliment strand, and a three excel files (DataSet.xlsx - original file, Dataset_70.xlsx - which is the updated file with information on the genes that can be prematurely stopped in the first 70% of the gene, and Dataset_70_19_to_16.xlsx - which is the same as Dataset_70.xlsx but with the editing window of -19 to -16 (rather than to -11)).
The PyCharm library with all the Python files can be found in the "Program Files" directory.
The Conda environment file is also included (Base_Editing_ENV.yml).
