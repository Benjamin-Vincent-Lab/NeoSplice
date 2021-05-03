import pysam
import math
import argparse


def get_len(bam_file_tumor, bam_file_normal):
    bam_file = pysam.AlignmentFile(bam_file_tumor, 'rb')
    read_lens_tumor = []
    for read in bam_file.fetch("chr22"):
        read_lens_tumor.append(read.query_length)
    bam_file.close()
    bam_file = pysam.AlignmentFile(bam_file_normal, 'rb')
    read_lens_normal = []
    for read in bam_file.fetch("chr22"):
        read_lens_normal.append(read.query_length)
    bam_file.close()
    return int(math.floor(min(read_lens_tumor + read_lens_normal) * 0.9))


def main():
    parser = argparse.ArgumentParser(description="Utility for determining max k-mer search length.")
    parser.add_argument('bam_file_tumor', type=str, nargs='?', help='provide tumor bam file path here')
    parser.add_argument('bam_file_normal', type=str, nargs='?', help='provide normal bam file path here')
    args = parser.parse_args()
    print get_len(args.bam_file_tumor, args.bam_file_normal)


if __name__ == '__main__':
    main()