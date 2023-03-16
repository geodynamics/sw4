#!/usr/bin/env python3
import sys
import hatchet as ht

def main():
    # add /usr/workspace/wsb/ramesh/Project6/MAMMOTH/MERGE_TESTS/hatchet-2022.2.2 to PYTHONPATH on Tioga
    print("Hatchet version ::",ht.__version__)
    filename=sys.argv[1]
    print("Opening ",filename)
    gf = ht.GraphFrame.from_caliperreader(filename)
    print(gf.tree())

    with open("test.dot", "w") as dot_file:
        dot_file.write(gf.to_dot())

    # Drop all index levels in the DataFrame except ``node``.
    gf.drop_index_levels()
        
    # Group DataFrame by ``name`` column, compute sum of all rows in each
    # group. This shows the aggregated time spent in each function.
    #gf.dataframe.rename(columns = {'sum#sum#time.duration':'time'}, inplace = True)
    grouped = gf.dataframe.groupby('name').sum()
    counts=gf.dataframe['name'].value_counts()
    
    # Sort DataFrame by ``time`` column in descending order.
    sorted_df = grouped.sort_values(by=['sum#sum#time.duration'],
                                   ascending=False)
    
    # Display resulting DataFrame.
    print(sorted_df)
    print(counts)
    cols=[]
    cols.append('sum#sum#time.duration')
    with open("profile.dat","w") as p_file:
        df2S = sorted_df.to_string(header=True, index=True,columns=cols)
        p_file.write(df2S)


if __name__ == "__main__":
    main()
 
