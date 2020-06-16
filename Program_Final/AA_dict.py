# This file contains the following function(s): make_codon, AA_determinator_df, addPAM, inversebp

from builtins import print

import math
import pandas as pd
pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)

AAdict = {"att": ["ile","I"], "atc": ["ile","I"], "ata": ["ile","I"], "ctt": ["leu","L"], "ctc": ["leu","L"],
           "cta": ["leu","L"], "ctg": ["leu","L"], "tta": ["leu","L"], "ttg": ["leu","L"], "gtt": ["val","V"],
           "gtc": ["val","V"], "gta": ["val","V"], "gtg": ["val","V"], "ttt": ["phe","F"], "ttc": ["phe","F"],
           "atg": ["met_START","M"], "tgt": ["cys","C"], "tgc": ["cys","C"], "gct": ["ala","A"], "gcc": ["ala","A"],
           "gca": ["ala","A"], "gcg": ["ala","A"], "ggt": ["gly","G"], "ggc": ["gly","G"], "gga": ["gly","G"],
           "ggg": ["gly","G"], "cct": ["pro","P"], "ccc": ["pro","P"], "cca": ["pro","P"], "ccg": ["pro","P"],
           "act": ["thr","T"], "acc": ["thr","T"], "aca": ["thr","T"], "acg": ["thr","T"], "tct": ["ser","S"],
            "tcc": ["ser","S"], "tca": ["ser","S"], "tcg": ["ser","S"], "agt": ["ser","S"], "agc": ["ser","S"],
           "tat": ["tyr","Y"], "tac": ["tyr","Y"], "tgg": ["trp","W"], "caa": ["gln","Q"], "cag": ["gln","Q"],
           "aat": ["asn","N"], "aac": ["asn","N"], "cat": ["his","H"], "cac": ["his","H"], "gaa": ["glu","E"],
           "gag": ["glu","E"], "gat": ["asp","D"], "gac": ["asp","D"], "aaa": ["lys","K"], "aag": ["lys","K"],
           "cgt": ["arg","R"], "cgc": ["arg","R"], "cga": ["arg","R"], "cgg": ["arg","R"], "aga": ["arg","R"],
           "agg": ["arg","R"], "taa": ["STOP","*"], "tag": ["STOP","*"], "tga": ["STOP","*"]}

def make_codon(listCCDS, dictranges, SequenceTop, SequenceBottom):
    """returns a dataframe with the following columns: 'Position on Top Strand', '"C" Strand', 'locus_list',
    'locus_list Strand', 'AA Number', 'codon_list',  'C Position in codon_list', 'C Position in codon_list Print',
    'Range in Top Strand'"""

    print("Start make_codon")
    # listCCDS = a list made from all the editable Cs (near a PAM) on either the top or bottom sequence --> locC_T list
    # listCCDS has the structure = [[position c, c strand, locus, locus strand],[position c2, c2 strand, locus, locus strand]]
    # dictranges = {"locus":[start, stop],"locus2":[start, stop]} this is for all genes in the genome, numbers in terms of top strand

    #Strand can equal either "Top" or "Bottom"
    # ListCCDS and dictranges are in terms of the Top strand (numbering)
    # Sequence is in terms of the respective strand

    # Initalize all the lists that will be used to make the dataframe
    position_top_strand_list = [] # will be used as the first column
    which_strand_C_list = [] # will either have a + or - for top or bottom repectively
    locus_list = [] # will say which gene/CDS it is in
    which_strand_CDS_list = [] # will either have a + or - for top or bottom repectively
    AA_num_list = [] # will say which AA in the respective gene/CDS it is in
    codon_list = [] # Will say what the codon is on the coding strand (5'->3')
    range_top_list = [] # will give the range of the AA in terms of the top strand (5'->3')
    #RangeBottom = [] # will give the range of the AA in terms of the bottom strand (5'->3') only relavent for bottom strand
    position_codon_list = [] # will give you the position 0, 1, 2 of the c being investigated in the codon
    position_codon_print_list = [] # will give you the position 1, 2, 3 of the c being investigated in the codon; for user

    #UniqueLocus = [] # This list is to gather the CDS that can be edited (no duplicates) and how many sites they have

    for c in range(len(listCCDS)):
        # print(listCCDS_T[c])
        Cstrand = listCCDS[c][1]
        locus_strand = listCCDS[c][3]
        locus = listCCDS[c][2]
        positionC = listCCDS[c][0]

        if Cstrand == "+" and locus_strand == "+":
            # print("C position and locus_list are both on the TOP strand")
            listrange = dictranges[locus]  # should be the range as a list [start,stop]
            start = listrange[0]
            #stop = listrange_T[1]

            codoncal = (positionC - start) / 3 + 1
            #print("\n\nListCCDS[{}]:\t {}".format(c, ListCCDS[c]))
            # print("\n\nLocus: {}\t".format(ListCCDS[c][1]))
            # print("Position (Top Strand): {}\t".format(ListCCDS[c][0]))
            # #print("codoncal\t", codoncal)

            trucnum = math.trunc(codoncal)
            AA_number = trucnum  # this gives you the AA number of the CDS

            decimal = abs(round(codoncal - trucnum, 2))  # this gives you the decimal of the AA will be used to find BP of the AA
            positionCpy = positionC - 1  # gives the position to find in the Sequence file in python (numbering starts at 0)

            if decimal == 0.67:
                #print("equal than 0.67 take -2")
                bp_1 = SequenceTop[positionCpy - 1]
                bp_2 = SequenceTop[positionCpy - 2]
                codon = bp_2 + bp_1 + "c"
                codonposition = [positionC - 2, positionC]
                positioncodon = 2

            elif decimal == 0.33:
                #print("equal to 0.33, take +1 and -1")
                bp_1 = SequenceTop[positionCpy - 1]
                bp_2 = SequenceTop[positionCpy + 1]
                codon = bp_1 + "c" + bp_2
                codonposition = [positionC - 1, positionC + 1]
                positioncodon = 1

            else:
                #print("equal to 0 take +2")
                bp_1 = SequenceTop[positionCpy + 1]
                bp_2 = SequenceTop[positionCpy + 2]
                codon = "c" + bp_1 + bp_2
                codonposition = [positionC, positionC + 2]
                positioncodon = 0

            # print("Amino Acid number in CDS:\t", AA_number)
            # print("codon_list 5'->3':\t", codon)
            # print("Position 5'->3'(Top Strand):\t", codonposition)

            position_top_strand_list.append(positionC)
            which_strand_C_list.append(Cstrand)
            locus_list.append(locus)
            which_strand_CDS_list.append(locus_strand)
            AA_num_list.append(AA_number)
            codon_list.append(codon)
            range_top_list.append(codonposition)
            position_codon_list.append(positioncodon)
            position_codon_print_list.append(positioncodon+1)

        elif Cstrand == "-" and locus_strand == "-":
            # print("C position and locus_list are both on the BOTTOM strand")
            numbp = len(SequenceBottom) # the length of the sequence, will be used to adjust positioning in case of Bottom strand

            listrange = dictranges[locus]  # should be the range as a list [stop,start] in top strand (but because bottom they go the other direction)
            # start and stop are the opposite of that in the top strand!!!
            start = listrange[1]
            #stop = listrange_T[0]

            codoncal = (start - positionC) / 3 + 1
            #print("\n\nListCCDS[{}]:\t {}".format(c, ListCCDS[c]))
            # print("\n\nLocus: {}\t".format(ListCCDS[c][1]))
            # print("Position (Top Strand): {}\t".format(ListCCDS[c][0]))
            # #print("codoncal\t", codoncal)

            trucnum = math.trunc(codoncal)
            AA_number = trucnum  # this gives you the AA number of the CDS

            decimal = abs(round(codoncal - trucnum, 2))  # this gives you the decimal of the AA will be used to find BP of the AA

            if decimal == 0.67:
                #print("equal than 0.67 take +2")
                bp_1 = SequenceBottom[numbp-positionC-1]
                bp_2 = SequenceBottom[numbp-positionC-2]
                codon = bp_2 + bp_1+"c"
                # The position numbers are fine with out a correction
                codonposition = [positionC +2, positionC]
                codonpositionB = [numbp-positionC,numbp-positionC+2]
                positioncodon = 2

            elif decimal == 0.33:
                #print("equal to 0.33, take +1 and -1")
                bp_1 = SequenceBottom[numbp-positionC - 1]
                bp_2 = SequenceBottom[numbp-positionC + 1]
                codon = bp_1 + "c" + bp_2
                # The position numbers are fine with out a correction
                codonposition = [positionC + 1, positionC - 1]
                codonpositionB = [numbp-positionC - 1, numbp-positionC + 1]
                positioncodon = 1

            else:
                # print("equal to 0 take +2")
                bp_1 = SequenceBottom[numbp-positionC + 1]
                bp_2 = SequenceBottom[numbp-positionC + 2]
                # The position numbers are fine with out a correction
                codon = "c" + bp_1 + bp_2 # 5'-> 3'
                codonposition = [positionC, positionC-2]
                codonpositionB = [numbp-positionC-2, numbp-positionC]
                positioncodon = 0

            # print("Amino Acid number in CDS:\t", AA_number)
            # print("codon_list 5'->3':\t", codon)
            # print("Position 5'->3'(Top Strand):\t", codonposition)
            # print("Position 5'->3'(Bottom Strand):\t",codonpositionB)

            position_top_strand_list.append(positionC)
            which_strand_C_list.append(Cstrand)
            locus_list.append(locus)
            which_strand_CDS_list.append(locus_strand)
            AA_num_list.append(AA_number)
            codon_list.append(codon)
            range_top_list.append(codonposition)
            position_codon_list.append(positioncodon)
            position_codon_print_list.append(positioncodon + 1)


        elif Cstrand == "+" and locus_strand == "-":
            # print("C position is on the TOP strand and the locus_list on the BOTTOM strand")

            listrange = dictranges[locus]  # should be the range as a list [stop,start] in top strand (but because bottom they go the other direction)
            # start and stop are the opposite of that in the top strand!!!
            stop = listrange[1] # actually the stop value

            codoncalB = (stop - positionC) / 3 + 1
            trucnumB = math.trunc(codoncalB)
            AAnumberB = trucnumB  # this gives you the AA number of the CDS
            # print("AANUMBERB", AAnumberB) # This is the AA number because the CDS is on the BOTTOM Strand

            # Now the calculation for the codon in the TOP STRAND
            start = listrange[0]
            codoncal = (positionC - start) / 3 + 1
            trucnum = math.trunc(codoncal)

            decimal = abs(round(codoncal - trucnum, 2))  # this gives you the decimal of the AA will be used to find BP of the AA
            positionCpy = positionC - 1  # gives the position to find in the Sequence file in python (numbering starts at 0)

            if decimal == 0.67:
                # print("equal than 0.67 take -2")
                bp_1 = SequenceTop[positionCpy - 1]
                bp_2 = SequenceTop[positionCpy - 2]
                codon = bp_2 + bp_1 + "c"
                codonposition = [positionC - 2, positionC]
                positioncodon = 0
                bp_1b = inversebp(bp_1)
                bp_2b = inversebp(bp_2)
                codonbottom = "g" + bp_1b + bp_2b

            elif decimal == 0.33:
                # print("equal to 0.33, take +1 and -1")
                bp_1 = SequenceTop[positionCpy - 1]
                bp_2 = SequenceTop[positionCpy + 1]
                codon = bp_1 + "c" + bp_2
                codonposition = [positionC - 1, positionC + 1]
                positioncodon = 1
                bp_1b = inversebp(bp_1)
                bp_2b = inversebp(bp_2)
                codonbottom = bp_2b + "g" + bp_1b

            else:
                # print("equal to 0 take +2")
                bp_1 = SequenceTop[positionCpy + 1]
                bp_2 = SequenceTop[positionCpy + 2]
                codon = "c" + bp_1 + bp_2
                codonposition = [positionC, positionC + 2]
                positioncodon = 2
                bp_1b = inversebp(bp_1)
                bp_2b = inversebp(bp_2)
                codonbottom = bp_2b + bp_1b + "g"


            #print("Amino Acid number in CDS:\t", AA_number)
            # print("codon_list 5'->3':\t", codon)
            # print("Position 5'->3'(Top Strand):\t", codonposition)
            # print("codon_list Bottom 5'->3':\t", codonbottom)

            position_top_strand_list.append(positionC)
            which_strand_C_list.append(Cstrand)
            locus_list.append(locus)
            which_strand_CDS_list.append(locus_strand)
            AA_num_list.append(AAnumberB)
            codon_list.append(codonbottom)
            range_top_list.append(codonposition)
            position_codon_list.append(positioncodon)
            position_codon_print_list.append(positioncodon + 1)


        elif Cstrand == "-" and locus_strand == "+":
            # print("C position is on the BOTTOM strand and the locus_list on the TOP strand")

            listrange = dictranges[locus]  # should be the range as a list [stop,start] in top strand (but because bottom they go the other direction)

            # Calculated the AA number of the CDS in the TOP strand
            start = listrange[0]
            codoncal = (positionC - start) / 3 + 1
            trucnum = math.trunc(codoncal)
            AA_number = trucnum  # this gives you the AA number of the CDS

            # Now the CODON Calculation
            numbp = len(SequenceBottom)  # the length of the sequence, will be used to adjust positioning in case of Bottom strand

            # start and stop are the opposite of that in the top strand!!!
            start = listrange[1]
            # stop = listrange_T[0]

            codoncal = (start - positionC) / 3 + 1
            trucnum = math.trunc(codoncal)

            decimal = abs(round(codoncal - trucnum, 2))  # this gives you the decimal of the AA will be used to find BP of the AA

            if decimal == 0.67:
                # print("equal than 0.67 take +2")
                bp_1 = SequenceBottom[numbp - positionC - 1]
                bp_2 = SequenceBottom[numbp - positionC - 2]
                codon = bp_2 + bp_1 + "c"
                # The position numbers are fine with out a correction
                codonposition = [positionC + 2, positionC]
                codonpositionB = [numbp - positionC, numbp - positionC + 2]
                positioncodon = 0
                bp_1b = inversebp(bp_1)
                bp_2b = inversebp(bp_2)
                codontop = "g" + bp_1b + bp_2b

            elif decimal == 0.33:
                # print("equal to 0.33, take +1 and -1")
                bp_1 = SequenceBottom[numbp - positionC - 1]
                bp_2 = SequenceBottom[numbp - positionC + 1]
                codon = bp_1 + "c" + bp_2
                # The position numbers are fine with out a correction
                codonposition = [positionC + 1, positionC - 1]
                codonpositionB = [numbp - positionC - 1, numbp - positionC + 1]
                positioncodon = 1
                bp_1b = inversebp(bp_1)
                bp_2b = inversebp(bp_2)
                codontop = bp_2b + "g" + bp_1b

            else:
                # print("equal to 0 take +2")
                bp_1 = SequenceBottom[numbp - positionC + 1]
                bp_2 = SequenceBottom[numbp - positionC + 2]
                # The position numbers are fine with out a correction
                codon = "c" + bp_1 + bp_2  # 5'-> 3'
                codonposition = [positionC, positionC - 2]
                codonpositionB = [numbp - positionC - 2, numbp - positionC]
                positioncodon = 2
                bp_1b = inversebp(bp_1)
                bp_2b = inversebp(bp_2)
                codontop = bp_2b + bp_1b + "g"


            # print("Amino Acid number in CDS:\t", AA_number)
            # # print("codon_list 5'->3':\t", codon)
            # print("Position 5'->3'(Top Strand):\t", codonposition)
            # # print("Position 5'->3'(Bottom Strand):\t", codonpositionB)
            # print("codon_list Top 5'->3':\t", codontop)

            position_top_strand_list.append(positionC)
            which_strand_C_list.append(Cstrand)
            locus_list.append(locus)
            which_strand_CDS_list.append(locus_strand)
            AA_num_list.append(AA_number)
            codon_list.append(codontop)
            range_top_list.append(codonposition)
            position_codon_list.append(positioncodon)
            position_codon_print_list.append(positioncodon + 1)

    # intialise data of lists.
    data = {'Position on Top Strand': position_top_strand_list, '"C" Strand': which_strand_C_list, 'Locus': locus_list, 'Locus Strand': which_strand_CDS_list,
            'AA Number': AA_num_list, 'Codon': codon_list, 'C Position in Codon': position_codon_list,
            'C Position in Codon Print': position_codon_print_list, 'Range in Top Strand': range_top_list}

    # Create DataFrame
    df = pd.DataFrame(data)

    # Print the output.
    df

    return df

def AA_determinator_df(dataframeStrand):
    """Takes the dataframe from make_codon() and finds the corresponding amino acids (AA), both the original and the
    edited ones (3 letter code, 1 letter code, and they type of edit)"""

    print("Start AA_determinator_df")
    # DF header (columns)
    # 'Position on Top Strand', '"C" Strand', 'Locus', 'Locus Strand', 'AA Number', 'Codon', 'C Position in Codon', 'C Position in Codon Print', 'Range in Top Strand'

    oldAAlist = [] # list of the orginial AA
    newAAlist = [] # list of the new AA
    term_list = [] # list of what changes
    oldAAcodelist = [] # list of the orginal AA single letter codes
    newAAcodelist = [] # list of the new AA single letter codes
    new_codon_list = [] # list of the new codons 5'->3'

    # counters for the number of same amino acids, stop codons that are made and different amino acids
    count_same = 0
    count_stop = 0
    count_rest = 0

    # skip the first line because it has the headers
    for r in range(len(dataframeStrand)):
        positionincodon = dataframeStrand.iloc[r][6] # gives the position of the c of interest in the codon
        key = dataframeStrand.iloc[r][5] # pulls out the codon
        locus_strand = dataframeStrand.iloc[r][3]
        cstrand = dataframeStrand.iloc[r][1]
        aa = AAdict[key][0] # value of the dictionary # the amino acid 3 letter code
        aacode = AAdict[key][1] # gives the AA acid single letter code

        keylist = list(key) # make the codon into a list so that it can be mutable

        if cstrand == "+" and locus_strand == "+":
            keylist[positionincodon] = 't'
        elif cstrand == "-" and locus_strand == "-":
            keylist[positionincodon] = 't'
        elif cstrand == "+" and locus_strand == "-":
            keylist[positionincodon] = 'a'
        elif cstrand == "-" and locus_strand == "+":
            keylist[positionincodon] = 'a'

        # make new codon into a string
        newcodon = "".join(keylist)

        newAA = AAdict[newcodon][0]
        newAAcode = AAdict[newcodon][1]
        newAAlist.append(newAA)
        newAAcodelist.append(newAAcode)
        oldAAlist.append(aa)
        oldAAcodelist.append(aacode)
        new_codon_list.append(newcodon)

        # print("\n\nOld codon: ", key, "\nOld AA: ", aa)
        # print("New codon:\t",newcodon, "\nNew AA: ", newAA)

        if newAA == aa:
            # print("The amino acids are the same")
            count_same += 1
            # term = "Nonsense"
            term = "Silent"

        elif newAA == "STOP":
            # print("The new amino acid is a stop codon")
            count_stop += 1
            # term = "Nonsense (stop)"
            term = "Nonsense"

        else:
            # print("print amino acids are different")
            count_rest += 1
            # term = "Missense"
            term = "Missense"

        term_list.append(term)

    # print("\n\nCountSame:\t", count_same)
    # print("CountStop:\t", count_stop)
    # print("CountRest:\t", count_rest,"\n\n")

    return new_codon_list, oldAAlist, oldAAcodelist, newAAlist, newAAcodelist, term_list

def addPAM(listPAM, df):
    """Takes in the list with the C position and the PAM information, and combines it to a dictionary. This is done
    while preserving/listing the multiple PAM sequences for each C position if there are any – throughout the whole
    genome. It returns the following dictionary and three lists: CPAMdict, (key is the C position (top strand),
    the value has the list of the PAM number, the PAM range, the PAM sequence and the 20 bp before the PAM sequence)
    listseq (a list of the all the PAM sequences, a list (of lists) if there are multiple that can edit the same C)
    listrange (a list of the range of the PAM sequences, a list (of lists) if there are more than one for a given C)
    list20bp (the 20 bp before the PAM sequence, a list (of lists) if there are more than one for a given C)."""

    print("Start addPAM")
    #this program will first take in the list with the C position and the PAM information, and combine it to a dictionary
    #this is done while preserving/listing the multiple PAM sequences for each C position if there are any
    #this is throughout the whole genome
    CPAMdict = {}
    for i in range(len(listPAM)):
        c = listPAM[i][0]  # this is the position of the C in reference to the top strand
        # print("c ",c)
        if c in CPAMdict:  # this means there are multiple PAM sequences tha can edit the same C
            #templist = []
            # print(c, " is in the dictionary already")
            templist = CPAMdict[c] #this should be a list of lists will have values in it
            #print("templist and type", templist, type(templist))
            #print("listPAM_T[i][1]: ",listPAM_T[i][1])
            #CPAMdict_T[c] = [oldvalue, listPAM_T[i][1]]
            templist.append(listPAM[i][1])  # this is the additional PAM sequence
            CPAMdict[c] = templist
            #print("CPAMdict_T[c] ",CPAMdict_T[c])

        else:
            templist=[]
            # print("new key added. Key: ", c)
            #print("templist",templist)
            templist.append(listPAM[i][1])  # add the appropriate PAM sequence
            CPAMdict[c] = templist
            #print("CPAMdict_T[c] ",CPAMdict_T[c])

    # print("CPAMdictFinal",CPAMdict_T)

    #up until this point, it is throughout the whole geonome

    # print("START THE NEW LISTS")
    # now we go specifically
    listCposition = df.iloc[:, 0]  # order of the C positions from the dataframes
    listseq = []  # the PAM sequence
    listrange = []  # the range of the PAM sequence
    list20bp = []  # the 20 bp upstream of the sequence
    for c in listCposition:
        # print(c)
        #print(len(CPAMdict_T[c]))
        # print(CPAMdict_T[c]) # this should be a list of lists

        if len(CPAMdict[c]) == 1:  # there is only one PAM sequence that can be used to edit this C
            # print("There is only 1 PAM sequence")
            listseq.append(CPAMdict[c][0][2])
            listrange.append(CPAMdict[c][0][1])
            list20bp.append(CPAMdict[c][0][3])

        elif len(CPAMdict[c]) > 1:  # There are more than one PAM sequences that can be used to edit this C
            # print("there are {} pam sequences".format(len(CPAMdict_T[c])))
            tempseqlist = []
            temprangelist = []
            temp20bplist = []
            for i in range(len(CPAMdict[c])):  # the number of C that can be edited
                # print("CPAMDICT for multiple pam",CPAMdict_T[c][i])
                # print("length of the thing",len(CPAMdict_T[c][i]))

                tempseqlist.append(CPAMdict[c][i][2])
                temprangelist.append(CPAMdict[c][i][1])
                temp20bplist.append(CPAMdict[c][i][3])

            listseq.append(tempseqlist)
            listrange.append(temprangelist)
            list20bp.append(temp20bplist)

    return CPAMdict, listseq, listrange, list20bp

def inversebp(bp):
    """This program makes the dna complement of bps (when passed a base pair)"""
    if bp == "a":
        bp = 't'
    elif bp == "t":
        bp = 'a'
    elif bp == "g":
        bp = 'c'
    elif bp == "c":
        bp = 'g'
    return bp