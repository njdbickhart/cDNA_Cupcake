#!/usr/bin/env python

__version__ = '2.0'

from Bio import SeqIO
import os, sys

def main():
    from argparse import ArgumentParser
    parser = ArgumentParser("Convert fasta to fastq")
    parser.add_argument("fasta_filename", help="input fasta (must end with .fasta or .fa)")
    parser.add_argument("--fq", default=None, help="fastq file that contains QV information for input fasta. If given, seq IDs must be in format XXX/<start>_<end>_CCS.")
    args = parser.parse_args()

    input = args.fasta_filename
    fq_src = args.fq

    fa2fq(input, fq_src)

def fa2fq(input, fq_src):
    try:
        assert input.lower().endswith('.fasta') or input.lower().endswith('.fa')
    except AssertionError:
        print >> sys.stderr, "Input {0} does not end with .fasta or .fa! Abort".format(input)
        sys.exit(-1)

    output = input[:input.rfind('.')] + '.fastq'
    if os.path.exists(output):
        print >> sys.stderr, "Output file {0} already exists. Abort!".format(output)
        sys.exit(-1)

    if fq_src is None:
        f = open(output, 'w')
        for r in SeqIO.parse(open(input),'fasta'):
            r.letter_annotations['phred_quality'] = [60]*len(r.seq)
            SeqIO.write(r, f, 'fastq')
        f.close()
    else:
        f = open(output, 'w')
        ids_seen = dict((r.id[:r.id.rfind('/')],r.description) for r in SeqIO.parse(open(input),'fasta'))
        for r in SeqIO.parse(open(fq_src),'fastq'):
            p = r.id[:r.id.rfind('/')]
            if p in ids_seen:
                r.description = ids_seen[p]
                r.id = r.description.split()[0]
                s, e, junk = r.id.split('/')[-1].split('_')
                SeqIO.write(makerec(r, int(s), int(e)), f, 'fastq')
                del ids_seen[p]
        f.close()

    print >> sys.stderr, "Output written to", f.name
    return f.name


# ex: m150803_002149_42161_c100745121910000001823165807071563_s1_p0/14/31_1114_CCS
def makerec(r, s, e):
    """
    :param r: Bio.Seq record
    :param s: 0-based start
    :param e: 1-based end
    :return: new Bio.Seq record that is subset from s-e (if +) or (e-s) if - strand

    *note:* if s > e, is - strand and will reverse complement.
    """
    if s < e:
        return r[s:e]
    else:
        r2 = r[e:s].reverse_complement()
        r2.id = r.id
        r2.description = r.description
        return r2



if __name__ == "__main__":
    main()
