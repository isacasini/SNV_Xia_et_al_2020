# This file contains the following function(s): findPAM

import re

def findPAM(sequence, strand):
    """finds all the PAM sequences with the sequence “ngg” in the 5’->3’ direction. The function returns a list of
    lists with the following information: “Sequence #”, "PAM Sequence 5'->3'", "Start Position 5'->3'",
    "Stop Position 5'->3'", "20bp before PAM 5'->3'", "bp 2-10 5'->3'". The function determines if the file is the top
    or bottom strand in order to adjust all numbering to be in terms of the top strand (as is convention)."""

    if strand == "+":
        # initialize row
        row = 0

        # initialize the list (of lists)
        PAM_info = [['Sequence #', "PAM Sequence 5'->3'", "Start Position 5'->3'", "Stop Position 5'->3'",
                    "20bp before PAM 5'->3'", "bp 2-10 5'->3'"]]
        for i in re.finditer("([atc]gg)", sequence):  # The letters in "([atc]gg)" can be changed to alter the PAM
            # print('This is the start{}, this is the end {}, this is the group {}'.format(i.start(), i.end(), i.group()))
            start = int(i.start())
            bp = sequence[start - 20:start]
            # these extra the BP in the 2nd, 3rd, and 4th position from the 5' direction in the 20 BP before the PAM Sequence
            bp2 = sequence[start - 19]
            bp3 = sequence[start - 18]
            bp4 = sequence[start - 17]
            bp5 = sequence[start - 16]
            bp6 = sequence[start - 15]
            bp7 = sequence[start - 14]
            bp8 = sequence[start - 13]
            bp9 = sequence[start - 12]
            bp10 = sequence[start - 11]

            templist = [str(row + 1), i.group(), i.start() + 1, i.end(), bp,
                        str(bp2 + bp3 + bp4 + bp5 + bp6 + bp7 + bp8 + bp9 + bp10)]
            PAM_info.append(templist)
            row += 1

    elif strand == "-":
        # initialize column and row
        row = 0

        # adjustment for being bottom sequence
        numbp = len(sequence)

        # initialize the list (of lists)
        PAM_info = [['Sequence #', "PAM Sequence 5'->3'", "Start Position 5'->3'", "Stop Position 5'->3'",
                    "20bp before PAM 5'->3'", "bp 2-10 5'->3'"]]

        # print(type(PAM_info))
        for i in re.finditer("([atc]gg)", sequence):  # The letters in "([atc]gg)" can be changed to alter the PAM
            # print('This is the start{}, this is the end {}, this is the group {}'.format(i.start(), i.end(), i.group()))

            start = int(i.start())
            bp = sequence[start - 20:start]
            # these extra the BP in the 2nd, 3rd, and 4th position from the 5' direction in the 20 BP before the PAM Sequence
            bp2 = sequence[start - 19]
            bp3 = sequence[start - 18]
            bp4 = sequence[start - 17]
            bp5 = sequence[start - 16]
            bp6 = sequence[start - 15]
            bp7 = sequence[start - 14]
            bp8 = sequence[start - 13]
            bp9 = sequence[start - 12]
            bp10 = sequence[start - 11]

            # this is to give the actual position of the PAM sequence on the Bottom Strand based on the numbering from the top strand
            adj_start = numbp - i.start()
            adj_end = numbp - i.end() + 1

            templist = [str(row + 1), i.group(), adj_start, adj_end, bp,
                        str(bp2 + bp3 + bp4 + bp5 + bp6 + bp7 + bp8 + bp9 + bp10)]

            PAM_info.append(templist)

            row += 1

    return PAM_info