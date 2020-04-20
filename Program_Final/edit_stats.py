# This file contains the following function(s): find_duplicates, getnoneditable, countsCDSedits

import pandas as pd
pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)
from check_CDS import remove_dup

def find_duplicates(df_all):
    """finds the count of the “unique” number of nonsense, nonsense (stop), and missense edits for each CDS.
    “Not unique” is defined as the same type of edit at the same amino acid number within the CDS. This function
    returns a dataframe (df_CDS) with the CDSs as the row indexes and the unique edit types as the columns. """

    Locusarray = df_all.Locus.unique() # can index into

    dictCDS = {}  # The key is the CDS, and the values is the list of the row index
    dictCDScount = {}  # The key is the CDS, [the number of unique Nonsense (same), Missense edits, Nonsense stop]

    for CDS in Locusarray:
        # get index (row id) of the rows with this CDS
        listCDS = df_all.index[df_all['Locus'] == CDS].tolist()

        dictCDS[CDS] = listCDS

    for cds in dictCDS.keys():
        Templist = [] # [[AA number,[Change],[AA number, [Change]]
        countnonsense = 0
        countstop = 0
        countmissense = 0

        for row in dictCDS[cds]:  # for every row in the list of indexes
            #print(row)
            AAnum = df_all.loc[row, "AA Number"]
            OAA = df_all.loc[row, "Old AA Code"]
            NAA = df_all.loc[row, "New AA Code"]
            Term = df_all.loc[row, "Term"]
            tempcode = [AAnum,[OAA,NAA]]

            if tempcode not in Templist:
                Templist.append(tempcode)
                if Term == "Missense":
                    countmissense+=1
                elif Term == "Nonsense":
                    countnonsense+=1
                else:
                    countstop+=1
        #print(Templist)
        dictCDScount[cds] = [countnonsense,countmissense,countstop]

    countCDSunique_df = pd.DataFrame.from_dict(dictCDScount, orient='index')
    countCDSunique_df.rename(columns={0: 'Unique Nonsense', 1: 'Unique Missense', 2: 'Unique Nonsense (stop)'}, inplace=True)

    return countCDSunique_df

def getnoneditable(df_all, dictCDSinfo):
    """makes a dictionary of all the CDS that are not editable and the count based on top and bottom strand."""

    TotalCDSlist = dictCDSinfo.keys()  # list of all the genes
    EditedCDSlist = df_all.loc[:, "Locus"].to_list()
    EditedCDSlistnondup = remove_dup(EditedCDSlist)
    nonedit_dict = {}  # key is the CDS, [number of editing, strand]

    count_T = 0  # counter for the non-editable on top strand
    count_B = 0  # counter for the non-editable on bottom strand
    for gene in TotalCDSlist:
        if gene not in EditedCDSlistnondup:
            nonedit_strand = dictCDSinfo[gene][1]
            nonedit_dict[gene] = [nonedit_strand, 0]

            if nonedit_strand == "+":
                count_T +=1
            else:
                count_B +=1

    return nonedit_dict, count_T, count_B

def countsCDSedits(countCDSdf_all,df_all,dictCDSinfo):
    """Returns a dataframe (with all CDSs, their strand, the total number of editing sites, the unique nonsense,
    missense and nonsense (stop) edits), strand of the CDS, the number/count of prematurely stoppable CDSs in
    top strand, the number/count of prematurely stoppable CDSs in bottom strand, the number/count of non-editable CDSs
    on the top strand, the number/count of non-editable CDSs on the bottom strand. The function also prints the CDSs
    that are not stoppable and the total count. This function uses both find_duplicates and getnoneditable."""

    # Find which strand the CDSs are on
    # Makes a new list with the order of the CDSs in the combined CDS dataframe (countCDSdf_all)
    # Then you add the "strand_list" to the dataframe

    strand_list = []
    for i in countCDSdf_all.index:
        # i = the CDS/locus
        strand = dictCDSinfo[i][1]
        strand_list.append(strand)

    countCDSdf_all.insert(0, "Strand", strand_list, True)
    # gets a dictionary of non-editable genes and the number based on top/bottom strand
    noneditCDSdict, count_nonedit_T, count_nonedit_B = getnoneditable(df_all, dictCDSinfo)
    noneditCDS_df = pd.DataFrame.from_dict(noneditCDSdict, orient='index', columns=['Strand', "Total Number Editing Sites"])

    allCDS_df =pd.concat([countCDSdf_all, noneditCDS_df], sort=False)

    #find the unique counts
    countCDSunique_df  = find_duplicates(df_all)
    countCDS_df = pd.concat([allCDS_df, countCDSunique_df], axis=1, sort=False)

    # find total number of CDSs that are prematurely stopable
    count_stop_T = 0
    count_stop_B = 0

    CDScount_df_index = countCDS_df.index

    count_nonstop= 0
    for i in CDScount_df_index:
        value = countCDS_df.at[i, "Unique Nonsense (stop)"]
        if value > 0:
            # print("it can be stopped")
            if countCDS_df.at[i, "Strand"] == "+":
                count_stop_T += 1
            elif countCDS_df.at[i, "Strand"] == "-":
                count_stop_B += 1
            else:
                print(i, " can be stopped but strand unknown")
        else:
            print(i," can't be stopped")
            count_nonstop += 1

    countCDS_df.index.name = 'CDS'

    print("Number of non-stoppable CDS ", count_nonstop)

    return countCDS_df, strand_list, count_stop_T, count_stop_B, count_nonedit_T, count_nonedit_B