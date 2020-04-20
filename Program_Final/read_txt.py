# This file contains the following function(s): readtxt, readtxt2df

import pandas as pd

def readtxt(filepathname):
    """Skips a line that starts with “>” and then reads in the rest of the file while removing any new lines (“\n”)"""

    sequence = ""
    for line in open(filepathname, 'r'):
        if line[0] != ">":
            sequence += line.replace("\n","")

    return sequence


def readtxt2df(filepathname):
    """Makes a pandas dataframe from a tab delimited csv file, skipping the 1st row and using the 2nd as a header"""

    #  Needs to be a string with the path and file
    file = filepathname

    # skips the first row, sets the second row as header, uses tab as a delimiter
    data = pd.read_csv(filepathname, sep="\t", header=0 ,encoding="ISO-8859-1")#,dtype={"PubChem":int,"PubChem":np.nan})

    return data