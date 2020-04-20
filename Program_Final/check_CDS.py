# This file contains the following function(s): makeCDSdict, checkinCDS, remove_dup

import itertools

def makeCDSdict(dataCDS):
    """Takes the dataframe of the CDS information and puts that information on the ranges of the locus in lists
    and dictionaries to be used later."""

    print("Start makeCDSdict")
    # Top strand lists:
    dictrangesT = {} #initialize the dicitonary of key = locus, values = [start,stop]
    listlistrangeT = [] #[[key,[start,stop]],[key,[start,stop]]]

    # Bottom strand lists:
    dictrangesB = {} #initialize the dicitonary of key = locus, values = [start,stop]
    listlistrangeB = [] #[[key,[start,stop]],[key,[start,stop]]]

    # Combined strand lists:
    dictranges = {}  # initialize the dicitonary of key = locus, values = [start,stop]
    listlistrange = []  # [[key,[start,stop],direction],[key,[start,stop],direction]]

    # Combined dictionary with CDS info: the CDS will be the key and the other info will be a list [[start,stop],strand]
    dictCDSinfo = {}

    for r in range(len(dataCDS)): #for every row in dataCDS
        #print(dataCDS.loc[[r],["locus", "start ","stop"]])
        key= dataCDS["locus"][r]
        start = dataCDS["start"][r]
        stop = dataCDS["stop"][r]
        value = [start,stop] # list of the start and stop
        #this step will sort between CDS in top vs bottom strand
        direction = dataCDS["direction"][r]

        if direction == "=>":

            # add a new entry to the dictionary
            dictrangesT[key] = value
            dictranges[key] = value

            # makes a list of lists
            templist = [key, value]
            listlistrangeT.append(templist)

            templistT = [key, value,"+"]
            listlistrange.append(templistT)

            dictCDSinfo[key] = [value,"+"]

        elif direction == "<=":
            # add a new entry to the dictionary
            dictrangesB[key] = value
            dictranges[key] = value

            # makes a list of lists
            templist = [key, value]
            listlistrangeB.append(templist)

            # makes a list of lists
            templistT = [key, value,"-"]
            listlistrange.append(templistT)

            dictCDSinfo[key] = [value,"-"]

        else:
            print("Direction not given for: ", key)

    # return dictranges, listranges, listdictrange, listlistrange
    # return dictranges, listlistrange, dictrangesB, listlistrangeB, dictranges, listlistrange
    return dictranges, listlistrange, dictCDSinfo

def checkinCDS(listlistrange, locC, strandC):
    """checks to see which of the Cs in editable regions are also in a CDS (on either strand). It returns a list
    (listCCDS) with the following information: location of the C, strand of the C, locus, strand of the locus."""

    print("Start checkinCDS")
    # listlistranges is a list with this structure [[key,[start,stop]],[key,[start,stop]]], where key = locus
    # locC_T is the list of lists of the PAM sequence id# and then the position of the "C" that is nearby
    # the locC_T and the strand should match, ie, locC_T is strand dependent

    listCCDST = []
    listCCDStop = []
    listCCDSbottom = []
    counter = 0
    # if strandC == "TOP":
    #     stranduse = "+"
    # elif strandC == "BOTTOM":
    #     stranduse = "-"
    # else:
    #     print("strandC doesn't give proper value")
    # print("stranduse", stranduse)
    for c in range(len(locC)):
        testvalue = locC[c][1] # this is the position of the C value
        for r in range(len(listlistrange)):
            locus = listlistrange[r][0]
            start = listlistrange[r][1][0]
            stop = listlistrange[r][1][1]
            strand = listlistrange[r][2]  # this is the strand of the CDS

            if start <= testvalue <= stop: #stop position must be further down than start position
                # print("C in position: \t", testvalue, "\t is in the CDS:\t",locus)
                # print("locus has positions:\t",listlistrange[r][1])
                counter += 1
                stringtoprint = '{}\t{}\n'.format(testvalue, locus)
                listCCDST.append([testvalue, strandC, locus, strand])

                if strand == "+":
                    # print('this gene is on the top strand')
                    listCCDStop.append([testvalue,locus])
                elif strand == "-":
                    # print("this gene is on the bottom strand")
                    listCCDSbottom.append([testvalue,locus]) # the position of the c, the gene/CDS/locus

            # else:
            #     print("C in position: \t", testvalue, "\t is not in a CDS")
    print("checkinCDS finished, counter = ", counter)

    # return listCCDST, listCCDStop, listCCDSbottom
    return listCCDST

def remove_dup(listCCDS):
    """removes any duplicate values e.g., any C’s that can be edited by the same PAM sequence (thus show up >1x)."""
    # https://www.w3resource.com/python-exercises/list/python-data-type-list-exercise-69.php

    listCCDS.sort()
    noduplist = list(listoflists for listoflists,_ in itertools.groupby(listCCDS))

    return noduplist