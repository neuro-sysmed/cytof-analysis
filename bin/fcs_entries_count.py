#!/usr/bin/env python3
#
# Kim Brugger (17.08.2022) kbr(at)brugger.dk

import sys
import os
import pprint as pp
import fcsparser
import flowkit as fk
from bokeh.plotting import show
import matplotlib.pyplot as plt
_ = plt.ioff()

import argparse


def rev_dict(input:dict) -> dict:
    res = {}
    for k in input.keys():
        res[ input[k]] = k

    return res


def pns_to_pnn(pns:list, pnn:list) -> dict:
    res = {}
    for i in range(0, len(pns)):
        res[pns[i]] = pnn[i]

    return res

def main():

    parser = argparse.ArgumentParser(description='fcs_entries_counter tool')
    parser.add_argument('-t', '--time-points', default=False, action="store_true", help="print times for each entry point")
    parser.add_argument('filenames', nargs='*', help="List of files to find entries in")

    args = parser.parse_args()

    total_entries = 0

    for filename in args.filenames:
        if os.path.isdir( filename ):
            continue
        
        meta, data = fcsparser.parse(filename, meta_data_only=False)
        if args.time_points:
#            print(data)
            for t in data['Time']:
                print(t)
            
        entries = len(data)
        total_entries += entries
        print( f"{filename}: {entries}" )

    print(f"Total: {total_entries}")

    sys.exit()


#get_channel_number_by_label
if __name__ == "__main__":
    main()
