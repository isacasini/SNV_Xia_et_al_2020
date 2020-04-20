# This file contains the following function(s): countCedit

def countCedit(listPAMinfo, strand):
    """This program finds all instances of the base C in the 2nd --> 10th base pairs before the PAM sequence(counting
     backward from the PAM they are identified as base pairs -19 to -11)."""

    print("Start countCedit")
    if strand == "+":

        locC = []  # location of the C in terms of the top strand numbering (list)
        numC = 0  # number of Cs in the -19 to -11 (count)
        locCdict = {}  # this info in dictionary form, key = locC, value = [numC, startpam, stoppam, pamseq, bp20] (dict)
        listPAM = []  # this info in list form, [locC, [numC, [startpam, stoppam], pamseq, bp20]] (list)
        # skip the first line because it has the headers
        for pam in listPAMinfo[1:]:
            num = pam[0]  # pam id #
            # print("NUM EQUALS: ",num)
            # print(type(num)) num is a string
            startpam = pam[2]
            stoppam = pam[3]
            pamseq = pam[1]
            bp20 = pam[4]
            key = pam[5]  # the bp 2-10 before the pam sequence

            for l in range(len(key)):
                # l = the number within the 9 bp that are listed in the 2-10 (out of 20) sequence before the pam
                bp = key[l]
                if bp == "c":
                    # print("bp is c")
                    # print(l)
                    numC += 1
                    locofc = startpam - 19 + l  # this is the corrected location of the bp "c" in the 2-10 bp (out of 20) pre pam
                    locC.append([num, locofc])
                    locCdict[locofc] = [num, startpam, stoppam, pamseq, bp20]
                    listPAM.append([locofc, [num, [startpam, stoppam], pamseq, bp20]])

    elif strand == "-":

        locC = []
        numC = 0
        locCdict = {}
        listPAM = []
        # skip the first line because it has the headers
        for pam in listPAMinfo[1:]:
            num = pam[0]  # pam id #
            # print("NUM EQUALS: ",num)
            # print(type(num)) num is a string
            startpam = pam[2]
            stoppam = pam[3]
            pamseq = pam[1]
            bp20 = pam[4]
            key = pam[5]  # the bp 2-10 before the pam sequence

            for l in range(len(key)):
                # l = the number within the 9 bp that are listed in the 2-10 (out of 20) sequence before the pam
                bp = key[l]
                if bp == "c":
                    # print("bp is c")
                    # print(l)
                    numC += 1
                    locofc = startpam + 19 - l  # this is the corrected location of the bp "c" in the 2-10 bp (out of 20) pre pam
                    locC.append([num, locofc])
                    locCdict[locofc] = [num, startpam, stoppam, pamseq, bp20]
                    listPAM.append([locofc, [num, [startpam, stoppam], pamseq, bp20]])

    # numC = number of Cs in the -19 to -11 (countAAchange)
    # locC = location of the C in terms of the top strand numbering (list)
    # locCdict = this info in dictionary form (dict)
    # listPAM = this info in list form (list)

    return numC, locC, locCdict, listPAM
