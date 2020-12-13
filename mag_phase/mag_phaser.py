import os
import sys

try:
    import vcf
except ImportError:
    print("Cannot import vcf! Please install pyvcf!", file=sys.stderr)
    sys.exit(-1)
from Bio import SeqIO
import io.SAMMPileUpReader as sam
import io.MPileUpVariantCaller as VC
from io import VariantPhaser
from io import VariantPhaseCleaner
import argparse
import subprocess as sp
import re

MIN_COVERAGE = 10     # minimum number of FL reads for a gene to do SNP calling and phasing
ERR_SUB = 0.005
MAX_DIFF_ALLOWED = 3  # maximum difference in bases allowed for two haplotype strings
MIN_PERC_ALLOWED = .25  # minimum percent of total count for an allele
PVAL_CUTOFF = 0.01
MIN_AF_AT_ENDS = 0.10 # minimum minor allele frequency for SNPs at ends, which tend to have unreliable alignments

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("fastx_filename", help="Input FLNC fasta or fastq")
parser.add_argument("sam_filename")
parser.add_argument("mpileup_filename")
parser.add_argument("read_stat")
parser.add_argument("mapping_filename")
parser.add_argument("-o", "--output_prefix", required=True)
parser.add_argument("--strand", choices=['+', '-'], required=True)
parser.add_argument("--partial_ok", default=False, action="store_true")
parser.add_argument("-p", "--pval_cutoff", default=PVAL_CUTOFF, type=float)
parser.add_argument("-n", "--ploidy", type=int, default=2)
#parser.add_argument("-e", "--err_sub", default=ERR_SUB, type=float, help="Estimated substitution error rate (default: 0.005)")
def parse_user_input():
    parser = argparse.ArgumentParser(
            description = "A pipeline for aligning sequence data on a slurm cluster"
            )
    parser.add_argument('-a', '--assembly',
                        help="The mag assembly file in fasta format",
                        required=True, type=str
                        )
    parser.add_argument('-b', '--bamfile',
                        help="Aligned reads in bam file format [full path needed!]",
                        required=True, type=str
                        )
    parser.add_argument('-o', '--output',
                        help="output prefix",
                        required=True, type=str
                        )
    parser.add_argument('-g', '--genes',
                        help='SCG gene bed file',
                        required =True, type=str
                        )
    parser.add_argument('-p', '--pval_cutoff',
                        help="P value cutoff for variant calls",
                        default=PVAL_CUTOFF, type=float
                        )

    return parser.parse_args(), parser

def main(args, parser):
    args = parser.parse_args()

    # remove potential past run output
    past_files = [args.output_prefix+'.NO_SNPS_FOUND',
             args.output_prefix+'.NO_HAPS_FOUND',
             args.output_prefix+'.log',
             args.output_prefix+'.human_readable.txt',
             args.output_prefix+'.vcf',
             args.output_prefix+'.cleaned.human_readable.txt',
             args.output_prefix+'.cleaned.vcf']

    for file in past_files:
        if os.path.exists(file):
            os.remove(file)

    # (0) generate pileups
    mpileupFile, coords = elitePileups(args.bamfile, args.genes, args.assembly, args.output)

    # (1) read the mpileup and vall variants
    reader = sam.MPileUpReader(mpileupFile)
    recs = [r for r in reader]
    vc = VC.MPileUPVariant(recs, min_cov=MIN_COVERAGE, err_sub=ERR_SUB, expected_strand='+', pval_cutoff=args.pval_cutoff)
    vc.call_variant()
    print(vc.variant)

    if len(vc.variant) == 0:
        os.system("touch {out}.NO_SNPS_FOUND".format(out=args.output))
        print("No SNPs found. END.", file=sys.stderr)
        sys.exit(0)

    # (2) for each CCS read, assign a haplotype (or discard if outlier)
    pp = VariantPhaser.VariantPhaser(vc)
    pp.phase_variant(args.sam_filename, coords, args.output, partial_ok=args.partial_ok)
    pp.haplotypes
    pp.haplotypes.get_haplotype_vcf_assignment()

    # (3) Write phase to output files
    #seqids = set([r.id for r in SeqIO.parse(open(args.fastx_filename), VariantPhaser.type_fa_or_fq(args.fastx_filename))])
    #isoform_tally = VariantPhaser.phase_isoforms(args.read_stat, seqids, pp)
    #if len(isoform_tally) == 0:
    #    os.system("touch {out}.NO_HAPS_FOUND".format(out=args.output_prefix))
    #    print("No good haps found. END.", file=sys.stderr)
    #    sys.exit(0)
    pp.haplotypes.write_haplotype_to_vcf("tempfakefilename", dict(), args.output)

    # temporary exit
    sys.exit()
    # (4) clean isoforms
    hap_count = VariantPhaseCleaner.make_haplotype_counts(isoform_tally)

    # (5) error correct haplotypes
    #  if diploid, use exhaustive search
    #  otherwise, use hap counts (ToDo: make this work with exhaustive search later)
    variants = [ [base.upper() for base,count in vc.variant[pos]] for pos in pp.accepted_pos]

    if args.ploidy == 2 and all(len(vars)==2 for vars in variants):
        diff_arr, hap_count_ordered = VariantPhaseCleaner.infer_haplotypes_via_exhaustive_diploid_only(pp.haplotypes, variants)
    else:
        diff_arr, hap_count_ordered = VariantPhaseCleaner.infer_haplotypes_via_min_diff(pp.haplotypes.haplotypes, hap_count, args.ploidy, MAX_DIFF_ALLOWED, MIN_PERC_ALLOWED)

    if diff_arr is None:
        os.system("touch {out}.cleaned.NO_HAPS_FOUND".format(out=args.output_prefix))
        print("No good haps found. END.", file=sys.stderr)
        sys.exit(0)

    m, new_hap, new_isoform_tally = VariantPhaseCleaner.error_correct_haplotypes(pp.haplotypes, isoform_tally, diff_arr, hap_count_ordered)
    # write out the mapping relationship between: FL CCS --> (pre-corrected) hap --> error-corrected hap
    with open(args.output_prefix+'.cleaned.hap_info.txt', 'w') as f:
        f.write("id,hap_preclean,hap_postclean\n")
        for seqid, old_i in pp.seq_hap_info.items():
            f.write("{0},{1},{2}\n".format(seqid, pp.haplotypes.haplotypes[old_i], new_hap.haplotypes[m[old_i]]))


    new_hap.get_haplotype_vcf_assignment()
    new_hap.write_haplotype_to_vcf(args.mapping_filename, new_isoform_tally, args.output_prefix+'.cleaned')

def elitePileups(bam : str, elites : str, assembly : str, outprefix : str) -> str:
    # Get the sample name from the basename of the bed file
    bfile = os.path.basename(bam)
    bsegs = re.split("\.", bfile)
    change = re.compile('[:-]')
    files = []
    coords = []
    with open(elites, 'r') as bedfh:
        for region in getFormatBedLine(bedfh):
            cv = re.split(r'(.+):(\d+)-(\d+)', region)
            coords.append(list(cv[0], int(cv[1]), int(cv[2])))
            safe = re.sub(change, '_', region)
            cmd = ["samtools", "mpileup", "-r",
                   region, "-f", assembly, bam ]
            print(f'Cmd: {" ".join(cmd)}')
            sp.run(cmd, check = True, stdout=open(bsegs[0] + "." + safe + ".pileup", "w"))
            files.append(bsegs[0] + "." + safe + ".pileup")

    retpileup = f'{outprefix}.{bsegs[0]}.pileup'
    with open(retpileup, 'w') as out:
        for f in files:
            with open(f, 'r') as fh:
                for l in fh:
                    l = l.rstrip()
                    out.write(l + '\n')
            os.remove(f)
    return retpileup, coords

def getFormatBedLine(bedfh : TextIO) -> str:
    for l in bedfh:
        segs = l.split()
        yield segs[0] + ":" + segs[1] + "-" + segs[2]

if __name__ == "__main__":
    args, parser = parse_user_input()
    main(args, parser)
