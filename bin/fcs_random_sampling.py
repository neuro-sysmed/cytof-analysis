#!/usr/bin/env python3
#
# Kim Brugger (17.08.2022) kbr(at)brugger.dk

import sys
import os
import pprint as pp
import flowkit as fk

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
    parser.add_argument('-n', '--number', default=55, help="Number of points to randomly extract")
    parser.add_argument('-o', '--outfile-pattern', default="{sample}_rand{number}.fcs", help="Filename pattern")
    parser.add_argument('filenames', nargs='*', help="List of files to find entries in")

    args = parser.parse_args()

    total_entries = 0

    random_seed = 11

    for filename in args.filenames:
        if os.path.isdir( filename ):
            continue

        sce = fk.Sample(filename, cache_original_events=True)
#       print(sce)
        sce.subsample_events(args.number, random_seed)
        print(sce)

        outfile_name = args.outfile_pattern.format(sample='t/test', number=args.number )
        sce.export( outfile_name, source='orig', subsample=True)

#        print(f'Number of events: {list(sce.get_events())}')

        print( sce.get_events(source='raw') )

    sys.exit()


#get_channel_number_by_label
if __name__ == "__main__":
    main()
