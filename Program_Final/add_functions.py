# This file contains the following functions:
# They are used to calculate which premature stop codons occur in the first 70% of the gene

from check_CDS import *

def cds_aa(dictCDSinfo):
    ''' This function takes in a dictionary with the CDS/locus (key) and a list which contains the range of the locus
    (also in a list) and whether the locus coding strand is the top (+) or bottom (-). The CDS is kept as the key and
    the number of amino acids in the locus is appended the list of the value (value = [[start position, stop position],
     “strand (either “+” or “-“)])'''

    # get the list of keys (all the locuses)
    locus_list = list(dictCDSinfo.keys())

    # create new dictionary where the 3 item in the list of values is the number of amino acids
    dictCDSinfo_new = {}

    for l in locus_list:
        old_value = dictCDSinfo[l]  # a list with the range (list), and the strand
        # print("old_value",old_value)
        range_locus = dictCDSinfo[l][0]
        # print("range_locus",range_locus)

        # strand_locus = dictCDSinfo[l][1]

        start_locus = range_locus[0]  # start position (numbering based on top strand)
        stop_locus = range_locus[1]  # stop position (numbering based on top strand)

        # takes the absolute value because if the CDS is coded in the bottom strand then it will be a negative number
        num_aa = abs(stop_locus - start_locus)/ 3  # calculates how many AAs in the locus
        # print("num_aa",num_aa)
        old_value.append(num_aa)
        # print("new_value",old_value)
        dictCDSinfo_new[l] = old_value

    return dictCDSinfo_new

def first_70(dictCDSinfo_new,df_all):
    '''This program will go through all the edits and if the edit creates a premature stop codon then it will check
    if the codon is in the first 70% of the gene.
    It will count how many unique instances there are sorted by strand - based on the CDS.'''

    # make a list of the index to loop through all the genes

    index_edit = list(df_all.index)
    cds_stop_t = []  # initializes the list of the CDSs in the TOP strand that can be stopped in the 1st 70% of the gene
    cds_stop_b = []  # initializes the list of the CDSs in the TOP strand that can be stopped in the 1st 70% of the gene

    cut_off = 0.7  # we want all edits that are in the first 70% of the CDS
    for edit in index_edit:
        # checks the column with the title: "Replaced AA"
        # edit_type =df_all["Replaced AA"][edit]
        edit_type =df_all["New AA"][edit]  # new name is "Replaced AA"

        # checks if the edit is stop
        if edit_type == "STOP":
            aa_num = df_all["AA Number"][edit]  # this gets the amino acid number in the locus, new name "AA Position"
            locus = df_all["Locus"][edit]  # this gets the locus where there is the edit, new name is "CDS/Gene"

            locus_aa = dictCDSinfo_new[locus][2]  # this gets the number of amino acids in the locus

            # defines the strand, as the strand of the CDS
            # rather than the strand on which the edit occurs?
            strand_edit = df_all["Locus Strand"][edit]  # new name is "CDS Strand"

            # calculates the what part of the first part of the CDS the edit is in
            percent = aa_num/locus_aa

            if percent <= cut_off:
                if strand_edit == "+":
                    cds_stop_t.append(locus)
                elif strand_edit == "-":
                    cds_stop_b.append(locus)
    unique_cds_stop_t = remove_dup(cds_stop_t)
    unique_cds_stop_b = remove_dup(cds_stop_b)
    return unique_cds_stop_t, unique_cds_stop_b




