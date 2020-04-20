# This file contains the following function(s): numAAchanges, make_matrix

import pandas as pd
pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)
import numpy as np
from collections import Counter
from check_CDS import remove_dup



def numAAchanges(df_all):
    """Counts the number of types a certain edit occurs from x AA to y AA."""

    change_list = []  # ["Old","New"]

    for index, row in df_all.iterrows():
        # access data using column names
        #print(index, row['Old AA Code'], row['New AA Code'], row['Term'])
        templist = [row['Old AA Code'], row['New AA Code']]
        change_list.append(templist)

    # convert to tuple
    change_tup = map(tuple, change_list)
    count = Counter(change_tup)
    return count

def make_matrix(countAAchange_df):
    """Creates a matrix with the old amino acids as rows and the new/edited amino acids as columns, and fills in the
    number of edits of each type there are. It returns a dataframe with a structure that can be used as a heatmap. """

    originalAAlist=[]
    editAAlist = []

    countAAchange_df_index = countAAchange_df.index
    for i in countAAchange_df_index:
        originalAAlist.append(i[0])
        editAAlist.append(i[1])

    #Sorts them alphabetically
    originalAAlist.sort()
    editAAlist.sort()

    #Removes duplicates
    listoriginalnodup = remove_dup(originalAAlist)

    # Makes dataframe with the original amino acids as the rows and the new as the columns, starts by filling all w/ zero
    zero_df = pd.DataFrame(np.zeros((len(listoriginalnodup), len(listoriginalnodup))), columns=listoriginalnodup,
                          index=listoriginalnodup)

    for i in countAAchange_df_index:
        value = countAAchange_df.at[i, "countAAchange"]
        r = i[0] #(original)
        c = i[1] # (edited)
        zero_df.loc[r, c] = value #replaces the zero with the value from the countAAchange

    return zero_df