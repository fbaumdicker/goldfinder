#!/usr/bin/env python3
"""Author: Anne Kupczok

Description: Script to combine pangenomes and interaction matrix
"""

import argparse

def parsing_cmd_line():

    parser = argparse.ArgumentParser(description="Combine pangenomes and interaction matrix")

    parser.add_argument("-i", "--interactions", type=str, required=True,
                        help="Interaction matrix (empty cell or 0 means negative interactions, all other entries indicate positive interaction)", dest="inter")
    parser.add_argument("-s", "--separator", type=str, required=False,
                        help="Separator in interactions matrix, default ','", dest="sep", default=",")
    parser.add_argument("-j", "--phage_column", help="Interaction matrix has phages as columns and bacteria as lines instead of phages as lines and bacteria as columns (default)",default=False,action='store_true')
    parser.add_argument("-b","--bact_pan", type=str, required=True, help="Bacterial pangenome in roary format",dest="bact")
    parser.add_argument("-p","--phage_pan", type=str, required=False, help="Phage pangenome in roary format (optional)",dest="phag", default="")
    parser.add_argument("-m","--bact_map", type=str, required=False, help="Map (tab-separated file) with bacteria names in interaction matrix (first column) and pangenome (second column)",dest="bact_map",default="")
    parser.add_argument("-n","--phage_map", type=str, required=False, help="Map (tab-separated file) with phage names in interaction matrix (first column) and pangenome (second column)",dest="bact_map",default="")
    parser.add_argument("-o", "--out", type=str, required=True, help="output file", dest="outf")
    parser.add_argument("-r", "--out_rev", type=str, required=False, help="output file for reverse encoding of phages", dest="outf_rev",default="")

    return parser

def main():
    parser = parsing_cmd_line()
    args = parser.parse_args()

    #todo: read phage map
    #todo: phage pangenome

    #if maps are present read them
    map={} #key: name in pangenome, value: name in matrix
    if args.bact_map:
        for line in open(args.bact_map):
            spl=line.split()
            map[spl[1]]=spl[0]
        print("Read map of {} bacteria".format(len(map)))

    #read interaction matrix
    pids=[]
    phages={} #key: phage, value: list of bact
    if args.phage_column:
        b=0
        for line in open(args.inter):
            spl=line.rstrip().split(args.sep)
            if not pids: pids=spl
            else:
                if spl:
                    b+=1
                    ids=[i for i in range(1,len(spl)) if spl[i] not in ["","0"]]
                    aphages=[pids[i] for i in ids] #list of phages
                    for p in aphages:
                        if p in phages: phages[p].append(spl[0])
                        else:
                            #print(p)
                            phages[p]=[spl[0]]
        porg=len(pids)-1
    else:
        porg=0
        for line in open(args.inter):
            spl = line.rstrip().split(args.sep)
            if not pids: pids=spl
            else:
                if spl:
                    porg+=1
                    ids = [i for i in range(1, len(spl)) if spl[i] not in ["", "0"]]
                    phages[spl[0]]=[pids[i] for i in ids]
        b=len(pids)-1
    print ("Read interaction matrix of {} phages and {} bacteria, continue with {} phages with positive interactions".format(porg,b,len(phages)))

    #if no phage pangenome - modify pangenome and exit
    gnames = []
    outf = open(args.outf, 'w')
    outf2 = None
    if args.outf_rev:
        outf2=open(args.outf_rev, 'w')
    for line in open(args.bact):
        if not gnames:
            gnames = line.rstrip().split(",")
            if map: gnames[14:]=[map[g] for g in gnames[14:]]
            #print(gnames)
        outf.write(line)
        if outf2: outf2.write(line)
    for p in phages:
        aline = [""] * len(gnames)
        aline[0] = p + "_Phag"
        found = False
        for g in phages[p]:
            if g in gnames:
                aline[gnames.index(g)] = "x"
                found = True
        if found:
            outf.write(",".join(aline) + "\n")
            if outf2:
                aline2 = ["y" if x == "" else x for x in aline]
                aline2 = ["" if x == "x" else x for x in aline2]
                outf2.write(",".join(aline2) + "\n")

    #read phage pangenome

    #write pangenome



if __name__ == '__main__':
    main()