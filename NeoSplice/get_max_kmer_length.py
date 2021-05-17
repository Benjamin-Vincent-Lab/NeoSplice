import pysam
import math
import argparse


def get_len(tumor_bam, normal_bam):
    bam_file = pysam.AlignmentFile(tumor_bam, 'rb')
    read_lens_tumor = []
    for read in bam_file.fetch("chr22"):
        read_lens_tumor.append(read.query_length)
    bam_file.close()
    bam_file = pysam.AlignmentFile(normal_bam, 'rb')
    read_lens_normal = []
    for read in bam_file.fetch("chr22"):
        read_lens_normal.append(read.query_length)
    bam_file.close()
    return int(math.floor(min(read_lens_tumor + read_lens_normal) * 0.9))


def main():
    parser = argparse.ArgumentParser(description="Utility for determining max k-mer search length.")
    parser.add_argument('--tumor_bam', required=True, type=str, nargs='?', help='provide tumor bam file path here')
    parser.add_argument('--normal_bam', required=True, type=str, nargs='?', help='provide normal bam file path here')
    args = parser.parse_args()
    print get_len(args.tumor_bam, args.normal_bam)


if __name__ == '__main__':
    main()