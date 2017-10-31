#!/usr/bin/env python
import os, sys, re
from Bio import SeqIO

rex = re.compile('(\S+\/\d+)\/\S+')
def get_seqs_from_list(fastafile, listfile, by_zmw):
    if by_zmw:
        seqs = []
        for line in open(listfile):
            m = rex.match(line.strip())
            if m is None: 
                raise Exception, "expected ZMW format but saw {0} instead!".format(line.strip())
            seqs.append(m.group(1))
    else:
        seqs = [line.strip() for line in open(listfile)]
    for r in SeqIO.parse(open(fastafile), 'fasta'):
        if r.id in seqs or r.id.split('|')[0] in seqs or (by_zmw and r.id[:r.id.rfind('/')] in seqs):
            print ">" + r.id
            print r.seq

if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser("Get sequences from a fasta file from a list")
    parser.add_argument("fasta_filename", help="Input fasta filename to extract sequences from")
    parser.add_argument("list_filename", help="List of sequence IDs to extract")
    parser.add_argument("--by_zmw", action="store_true", default=False, help="Use ZMW name instead of full ID names")

    args = parser.parse_args()

    get_seqs_from_list(args.fasta_filename, args.list_filename, args.by_zmw)
