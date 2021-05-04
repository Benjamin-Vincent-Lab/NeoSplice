import pysam
import argparse
import collections
import logging


class CigarOp:
    """Enumeration of CIGAR codes."""
    MATCH = 0
    INS = 1
    DEL = 2
    REF_SKIP = 3
    SOFT_CLIP = 4


def retrieve_splices(read, splices, chromosome):
    ref_pos = read.reference_start
    read_pos = 0
    for op, count in read.cigartuples:
        if op == CigarOp.MATCH:
            ref_pos += count
            read_pos += count
        elif op == CigarOp.INS:
            read_pos += count
        elif op == CigarOp.DEL:
            ref_pos += count
        elif op == CigarOp.REF_SKIP:
            start = ref_pos
            end = ref_pos + count
            splices[chromosome, start, end] += 1
            ref_pos += count
        elif op == CigarOp.SOFT_CLIP:
            read_pos += count
        else:
            logging.warn('Unexpected cigar op {0}'.format(op))


def main():
    logging.basicConfig(level=logging.DEBUG, format="%(asctime)-15s [%(processName)s.%(levelname)s] %(message)s")
    parser = argparse.ArgumentParser(description="Utility for finding splice junctions from bam file.")
    parser.add_argument("input_bam", type=str, nargs='?', help="The input bam file")
    parser.add_argument("out_file", type=str, nargs='?', help="The output file")
    args = parser.parse_args()

    outf = open(args.out_file, 'w')
    samfile = pysam.AlignmentFile(args.input_bam, "rb")

    splices = collections.defaultdict(int)
    for chromosome in range(1, 23):
        chromosome = "chr" + str(chromosome)
        for read in samfile.fetch(chromosome):
            if not read.is_unmapped and not read.is_secondary and read.is_proper_pair and read.is_paired and\
                    not read.is_duplicate and not read.is_supplementary:
                retrieve_splices(read, splices, chromosome)

    for splice in splices:
        outf.write(splice[0] + '\t' + str(splice[1]) + '\t' + str(splice[2]) + '\t' + str(splices[splice]) + '\n')

    outf.close()
    samfile.close()


if __name__ == '__main__':
    main()