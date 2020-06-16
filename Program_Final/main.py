# This is the main file that call all other files in the Base Editing program.

from read_txt import *
from complement import *
from position_PAM import *
from all_editable import *
# from collections import Counter
from check_CDS import *
from AA_dict import *
# import numpy as np
from heatmap import *
from edit_stats import *
from add_functions import *

import pandas as pd
pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)

import time

start = time.time()

# STEP 1: Read in Top strand sequence file and make the dna_complement in the 5'->3' direction (the Bottom strand)
# USER_INPUT_REQUIRED - change the path in TopSequenceFile
# Name of the file with the top sequence
TopSequenceFile = r"E:\Colleague\Pengfei\Full_genome_clostridium_ljungdahlii_DSM13528.txt"
#r"E:\Colleague\Pengfei\Full_genome_clostridium_ljungdahlii_DSM13528_short2.txt"#r"E:\Colleague\Pengfei\Full_genome_clostridium_ljungdahlii_DSM13528.txt"
#returns file as a string
TopSequence = readtxt(TopSequenceFile)

# makes dna_complement sequence
BottomSequence = dna_complement(TopSequence) #this is the dna_complement and in the reverse
# USER_INPUT_REQUIRED - change the path in filecomp
# Writes out a file with the dna_complement sequence
filecomp = open(r"E:\Colleague\Pengfei\Comp_Sequence.txt",'w')
filecomp.write(BottomSequence)
filecomp.close()

print("The original (top) sequence runs 5' -> 3'. The first 10 bp and the last 10 bp are as follows:\n5' %s ... %s 3'"
      "\n\n" % (str(TopSequence[:10]), str(TopSequence[-10:])))
print("The bottom sequence runs 5' -> 3'. The first 10 bp and the last 10 bp are as follows:\n5' %s ... %s 3'"
      "" % (str(BottomSequence[:10]),str(BottomSequence[-10:])))

#------------------------------------------------------------------------------------------------------------
# STEP 2: Find all the PAM Sequence information (considering "ngg" 5'->3')
# this step will return a list of lists with [PAM ID #, PAM seq, start, stop, 20 bp upstream, bp 2-10*] *(-19 to -11)
# each list is strand dependent, therefore it needs to be run twice

#TOP SEQUENCE
listPAMinfo_T = findPAM(TopSequence, "+")  # this returns a list, in module "position_PAM.py"
# print("LIST INFO TOP\n", listPAMinfo_T)

#BOTTOM SEQUENCE
listPAMinfo_B = findPAM(BottomSequence, "-")  # this returns a list, in module "position_PAM.py"
# print("LIST INFO BOTTOM\n", listPAMinfo_B)

#------------------------------------------------------------------------------------------------------------
# STEP 3: Find all the Cs that are in the 2-10 (-19 to -11) bp range of the PAM sequence, done on each strand

# TOP SEQUENCE
numC_T, locC_T, locCdict_T, listPAM_T = countCedit(listPAMinfo_T, "+")
# BOTTOM SEQUENCE
numC_B, locC_B, locCdict_B, listPAM_B = countCedit(listPAMinfo_B, "-")

#------------------------------------------------------------------------------------------------------------
# STEP 4: This step reads in a text sheet with the CDS information from Snapgene and NCBI (CDS, CDS length, direction, start, stop, annotation)
pathtocds = r"E:\Colleague\Pengfei\Combined_annotation_page.txt"
dataCDS = readtxt2df(pathtocds) # this returns a data frame with info from an annotation sheet

#------------------------------------------------------------------------------------------------------------
# STEP 5: Makes a series of dictionaries and lists of the CDS and the range (start,stop), ["CDS",[start,stop]] or {"CDS":[start,stop]}
dictranges, listlistrange, dictCDSinfo = makeCDSdict(dataCDS) #dictranges, listranges, listdictrange, listlistrange

#------------------------------------------------------------------------------------------------------------
# STEP 6: This step makes a list that shows the [C position, C strand, Locus, Locus strand]

# Top Strand; is comparing the CDS to the "C"s in the top strand
listCCDS_T = checkinCDS(listlistrange, locC_T, "+")  # returns a list of [C position, C position strand ,locus,strand])
# Bottom Strand; is comparing the CDS to the "C"s in the bottom strand
listCCDS_B = checkinCDS(listlistrange, locC_B, "-")  # returns a list of [C position, C position strand ,locus,strand])

#------------------------------------------------------------------------------------------------------------
# STEP 7: Remove duplicates from the listCCDS_T/B, there are duplicates because the same "C" can be edited sometimes by multiple PAM sequences

# Top Strand
noduplist_T = remove_dup(listCCDS_T) # returns a list of [[location of C, CDS/Locus]]
# Bottom Strand
noduplist_B = remove_dup(listCCDS_B) # returns a list of [[location of C, CDS/Locus]]

#------------------------------------------------------------------------------------------------------------
# STEP 8: This will make a dataframe for each strand with the following headers
# (8 columns): 'Position on Top Strand', '"C" Strand', 'Locus', 'Locus Strand', 'AA Number', 'Codon',
# 'C Position in Codon', 'C Position in Codon Print', 'Range in Top Strand'

# Top Strand; this will check the C's in the top stand
df_T = make_codon(noduplist_T, dictranges, TopSequence, BottomSequence)
# Bottom Strand; this will check the C's in the bottom stand
df_B = make_codon(noduplist_B, dictranges, TopSequence, BottomSequence)

#------------------------------------------------------------------------------------------------------------
# STEP 9: Determine the Amino Acid (AA) (and the single code) that was orignially in the coding strand and then edited one
# AA_determinator_df pulls out the following information: 'C Position in Codon'

# Top Strand
new_codon_T, oldAA_T, oldAAcode_T, newAA_T, newAAcode_T, new_term_T = AA_determinator_df(df_T)
# Bottom Strand
new_codon_B, oldAA_B, oldAAcode_B, newAA_B, newAAcode_B, new_term_B = AA_determinator_df(df_B)

#------------------------------------------------------------------------------------------------------------
# STEP 10: Export the PAM data for the given C Positions (want the repeats etc)

# Top Strand
CPAMdict_T, listseq_T, listrange_T, list20bp_T = addPAM(listPAM_T, df_T)
# Bottom Strand
CPAMdict_B, listseq_B, listrange_B, list20bp_B = addPAM(listPAM_B, df_B)

#------------------------------------------------------------------------------------------------------------
# STEP 11: Add the new AA data to the dataframe
# this figures out which column to start adding the info to
# Combine the Term lists to be able to countAAchange them later

# In theory the number of coloumns in the top and bottom strand DF should be the same
numcol_T = len(df_T.columns)
# numcol_B = len(df_T.columns)
numcol_B = len(df_B.columns)

# Top Strand
df_T.insert(numcol_T + 0, "Old AA", oldAA_T, True)
df_T.insert(numcol_T + 1, "Old AA Code", oldAAcode_T, True)
df_T.insert(numcol_T + 2, "New AA", newAA_T, True)
df_T.insert(numcol_T + 3, "New AA Code", newAAcode_T, True)
df_T.insert(numcol_T + 4, "New Codon", new_codon_T, True)
df_T.insert(numcol_T + 5, "Term", new_term_T, True)
df_T.insert(numcol_T + 6, "PAM Sequence", listseq_T, True)
df_T.insert(numcol_T + 7, "PAM Range", listrange_T, True)
df_T.insert(numcol_T + 8, "20bp up PAM", list20bp_T, True)
# Bottom Strand
df_B.insert(numcol_B + 0, "Old AA", oldAA_B, True)
df_B.insert(numcol_B + 1, "Old AA Code", oldAAcode_B, True)
df_B.insert(numcol_B + 2, "New AA", newAA_B, True)
df_B.insert(numcol_B + 3, "New AA Code", newAAcode_B, True)
df_B.insert(numcol_B + 4, "New Codon", new_codon_B, True)
df_B.insert(numcol_B + 5, "Term", new_term_B, True)
df_B.insert(numcol_B + 6, "PAM Sequence", listseq_B, True)
df_B.insert(numcol_B + 7, "PAM Range", listrange_B, True)
df_B.insert(numcol_B + 8, "20bp up PAM", list20bp_B, True)

#------------------------------------------------------------------------------------------------------------
# STEP 12: Clean up of the DF, remove the 'C Position in Codon' column (the user doesn't care about it)
# Combine the two into 1 DF, this will be used for the CDS information page
# Sort the DF by the the C position (decending)

# Top Strand
df_T_print = df_T.drop(columns=['C Position in Codon'])
# Bottom Strand
df_B_print = df_B.drop(columns=['C Position in Codon'])

df_all = pd.concat([df_T_print, df_B_print], sort=False)

df_all = df_all.sort_values('Position on Top Strand')

# resets the index from 0-end value
df_all.reset_index(inplace = True, drop = True)

#------------------------------------------------------------------------------------------------------------
# STEP 13: Put together the summary information
# First it extracts all the locuses that can be edited from the DF
# Because each locus can be edited in multiple spots, we remove the duplicates this is the unique list
# But we are interested in how many sites are editable in each locus, so these are counted

# Top Strand
locus_list_T = list(df_T["Locus"])
unique_locus_list_T = remove_dup(locus_list_T)
countCDS_T = dict(Counter(locus_list_T))  # dictionary with locus:number of instances
# Bottom Strand
locus_list_B = list(df_B["Locus"])
unique_locus_list_B = remove_dup(locus_list_B)
countCDS_B = dict(Counter(locus_list_B))  # dictionary with locus:number of instances

#------------------------------------------------------------------------------------------------------------
# STEP 14: Turn information into a DF
# Top Strand
countCDSdf_T = pd.DataFrame.from_dict(countCDS_T, orient='index', columns=['Total Number Editing Sites'])
# print(countCDSdf_T)
# Bottom Strand
countCDSdf_B = pd.DataFrame.from_dict(countCDS_B, orient='index', columns=['Total Number Editing Sites'])
# print(countCDSdf_B)

#------------------------------------------------------------------------------------------------------------
# STEP 15: combine the CDS/number of sites between top and bottom strand into a DF
# This combines the data frames by adding together the values with the same index

countCDSdf_all = countCDSdf_T.combine(countCDSdf_B, np.add, fill_value=0)

#------------------------------------------------------------------------------------------------------------
# STEP 16: Get stats and dataframe with the following information

# countsCDSedits returns a CDS with their strand, the total number of editing sites, the unique nonsense, missense and stop edits,
# strand of the CDS, the num/count of prematurely stopable CDSs in top strand,
# the num/count of prematurely stopable CDSs in bottom strand, the num/count of non-editable CDSs on the top strand,
# the num/count of non-editable CDSs on the bottom strand
countCDS_df, strand_list, count_stop_T, count_stop_B, count_nonedit_T, count_nonedit_B = countsCDSedits(countCDSdf_all, df_all, dictCDSinfo)

#####
dictCDSinfo_new = cds_aa(dictCDSinfo)

# finds the number (unique) of CDS that are stopable in the first 70%
unique_cds_stop_t, unique_cds_stop_b = first_70(dictCDSinfo_new,df_all)
count_cds_stop_t = len(unique_cds_stop_t)  # top strand (from CDS perspective)
count_cds_stop_b = len(unique_cds_stop_b)  # bottom strand (from CDS perspective)

#------------------------------------------------------------------------------------------------------------
# STEP 17: Create information for the summary sheet of the workbook.
# Pull together stats information, the index (rows), top (columns)
# Then created the DF

# sum_index = ["Total Edits in Genome", "Total Number of Editable CDSs", "Total Number of Noneditable CDSs",
#              "Total Number of Stopable CDSs", "Total Number of Stopable CDSs in 1st 70%", "",
#              "Silent Mutation", "Nonsense Mutation", "Missense Mutation"]
# sum_T = [numC_T, strand_list.count("+"), count_nonedit_T, count_stop_T, count_cds_stop_t, "",
#          new_term_T.count("Silent Mutation"), new_term_T.count("Nonsense Mutation"),
#          new_term_T.count("Missense Mutation")]
# sum_B = [numC_B, strand_list.count("-"), count_nonedit_B, count_stop_B, count_cds_stop_b, "",
#          new_term_B.count("Silent Mutation"), new_term_B.count("Nonsense Mutation"),
#          new_term_B.count("Missense Mutation")]
# sum_all = [numC_T + numC_B, strand_list.count("+") + strand_list.count("-"), count_nonedit_T + count_nonedit_B,
#            count_stop_T + count_stop_B, count_cds_stop_t+count_cds_stop_b, "",
#            new_term_T.count("Silent Mutation") + new_term_B.count("Silent Mutation"),
#            new_term_T.count("Nonsense Mutation") + new_term_B.count("Nonsense Mutation"),
#            new_term_T.count("Missense Mutation") + new_term_B.count("Missense Mutation")]

sum_index = ["Total Edits in Genome", "Total Number of Editable CDSs", "Total Number of Noneditable CDSs",
             "Total Number of Stopable CDSs","Total Number of Stopable CDSs in 1st 70%", "", "Nonsense", "Nonsense (stop)", "Missense"]
sum_T = [numC_T, strand_list.count("+"), count_nonedit_T, count_stop_T, count_cds_stop_t,"", new_term_T.count("Silent"),
         new_term_T.count("Nonsense"), new_term_T.count("Missense")]
sum_B = [numC_B, strand_list.count("-"), count_nonedit_B, count_stop_B, count_cds_stop_b, "", new_term_B.count("Silent"),
         new_term_B.count("Nonsense"), new_term_B.count("Missense")]
sum_all = [numC_T + numC_B, strand_list.count("+") + strand_list.count("-"), count_nonedit_T + count_nonedit_B,
           count_stop_T + count_stop_B,count_cds_stop_t+count_cds_stop_b, "", new_term_T.count("Silent") + new_term_B.count("Silent"),
           new_term_T.count("Nonsense") + new_term_B.count("Nonsense"), new_term_T.count("Missense")
           + new_term_B.count("Missense")]
#





# Functioning one

# sum_index = ["Total Edits in Genome", "Total Number of Editable CDSs", "Total Number of Noneditable CDSs",
#              "Total Number of Stopable CDSs","Total Number of Stopable CDSs in 1st 70%", "", "Nonsense", "Nonsense (stop)", "Missense"]
# sum_T = [numC_T, strand_list.count("+"), count_nonedit_T, count_stop_T, count_cds_stop_t,"", new_term_T.count("Nonsense"),
#          new_term_T.count("Nonsense (stop)"), new_term_T.count("Missense")]
# sum_B = [numC_B, strand_list.count("-"), count_nonedit_B, count_stop_B, count_cds_stop_b, "", new_term_B.count("Nonsense"),
#          new_term_B.count("Nonsense (stop)"), new_term_B.count("Missense")]
# sum_all = [numC_T + numC_B, strand_list.count("+") + strand_list.count("-"), count_nonedit_T + count_nonedit_B,
#            count_stop_T + count_stop_B,count_cds_stop_t+count_cds_stop_b, "", new_term_T.count("Nonsense") + new_term_B.count("Nonsense"),
#            new_term_T.count("Nonsense (stop)") + new_term_B.count("Nonsense (stop)"), new_term_T.count("Missense")
#            + new_term_B.count("Missense")]
# #








# sum_index = ["Total Edits in Genome", "Total Number of Editable CDSs", "Total Number of Noneditable CDSs",
#              "Total Number of Stopable CDSs", "", "Nonsense", "Nonsense (stop)", "Missense"]
# sum_T = [numC_T, strand_list.count("+"), count_nonedit_T, count_stop_T, "", new_term_T.count("Nonsense"),
#          new_term_T.count("Nonsense (stop)"), new_term_T.count("Missense")]
# sum_B = [numC_B, strand_list.count("-"), count_nonedit_B, count_stop_B, "", new_term_B.count("Nonsense"),
#          new_term_B.count("Nonsense (stop)"), new_term_B.count("Missense")]
# sum_all = [numC_T + numC_B, strand_list.count("+") + strand_list.count("-"), count_nonedit_T + count_nonedit_B,
#            count_stop_T + count_stop_B, "", new_term_T.count("Nonsense") + new_term_B.count("Nonsense"),
#            new_term_T.count("Nonsense (stop)") + new_term_B.count("Nonsense (stop)"), new_term_T.count("Missense")
#            + new_term_B.count("Missense")]



# These give the column headings
sum_all_dict = {'Top Strand': sum_T, 'Bottom Strand': sum_B, "Total": sum_all}
sum_all_df = pd.DataFrame(sum_all_dict, index=sum_index)

#------------------------------------------------------------------------------------------------------------
# STEP 18: Find the number of instances of each type of change (old amino acid to new amino acid)

countAAchange = numAAchanges(df_all)
countAAchange_df = pd.DataFrame.from_dict(countAAchange, orient='index')
countAAchange_df = countAAchange_df.rename(columns={0: 'countAAchange'})
countAAchange_df.index.name = 'Original AA, New AA'
count_matrix_df = make_matrix(countAAchange_df)

#------------------------------------------------------------------------------------------------------------
# STEP 19: Edit the column names of the dataframes for the final excel sheet to make more user-friendly
countCDS_df.rename(columns={'Unique Nonsense': 'Silent Mutations', 'Unique Missense': 'Missense Mutations',
                            'Unique Nonsense (stop)': 'Nonsense Mutations'}, inplace=True)
df_all.rename(columns={'"C" Strand': 'Editing Strand ("+" indicates top strand in the DNA file)',
                       'Locus': 'CDS/Gene','Locus Strand': 'CDS Strand','AA Number': 'AA Position',
                       'Old AA': 'Original AA','Old AA Code': 'Original AA (One letter)','New AA': 'Replaced AA',
                       'New AA Code': 'Replaced AA (One letter)','New Codon': 'Replaced Codon','Term': 'Mutation Type',
                       'PAM Range': 'PAM Position','20bp up PAM': 'Protospacer'}, inplace=True)

sum_all_df.rename(index={'Nonsense': 'Silent Mutations', 'Missense': 'Missense Mutations',
                            "Nonsense (stop)": 'Nonsense Mutations'}, inplace=True)

#------------------------------------------------------------------------------------------------------------
# STEP 20: Write out all the information to a workbook in 5 different sheets
# USER_INPUT_REQUIRED - change the path in writePathout
writePathout = r'E:\Colleague\Pengfei\Program_Final\Dataset_70.xlsx'
writer = pd.ExcelWriter(writePathout,engine='xlsxwriter')

sum_all_df.to_excel(writer, sheet_name='Summary', header=True, index=True, index_label=None)
countCDS_df.to_excel(writer, sheet_name='Base Editing in CDS', header=True, index=True, index_label=None)
count_matrix_df.to_excel(writer, sheet_name='Amino Acid Replacement', startcol=1,startrow =1, header=True, index=True, index_label=None)
df_all.to_excel(writer, sheet_name='List of Editable Sites', header=True, index=True, index_label=None)
countAAchange_df.to_excel(writer, sheet_name='Heatmap', header=True, index=True, index_label=None)

# Want to add labels on each side of the "matrix" in the sheet Amino Acid Replacement
# Get the xlsxwriter workbook and worksheet objects.
workbook = writer.book
worksheet = writer.sheets['Amino Acid Replacement']
cell_format = workbook.add_format({'bold': 1,'border': 1,'align': 'center','valign': 'vcenter'})
cell_format.set_rotation(90)
cell_format2 = workbook.add_format({'bold': 1,'border': 1,'align': 'center','valign': 'vcenter'})

# Merge cells
worksheet.merge_range('A3:A23', 'Original Amino Acids', cell_format)
worksheet.merge_range('C1:W1', 'Replaced Amino Acids', cell_format2)

writer.save()

#------------------------------------------------------------------------------------------------------------
# STEP 21: Check that the program exited correctly
# Check and print how long the program takes
end = time.time()
print("\n\nprocess exited successfully")
print("It took: ", (end-start)/60, " minutes")