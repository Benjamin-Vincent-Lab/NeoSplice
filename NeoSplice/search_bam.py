import pysam
import ahocorasick
import argparse
import logging
import timeit


def main():
    start = timeit.default_timer()
    logging.basicConfig(level=logging.DEBUG, format="%(asctime)-15s [%(processName)s.%(levelname)s] %(message)s")
    parser = argparse.ArgumentParser(description="Utility for retrieving tumor specific kmers in RNA-seq reads.")
    parser.add_argument("Kmer_file", type=str, nargs='?', help="provide Kmer file here")
    parser.add_argument("input_bam_file", type=str, nargs='?', help="provide input bam file path here")
    parser.add_argument("out_bam_file", type=str, nargs='?', help="provide output bam file path here")
    args = parser.parse_args()

    cigar_map = {0:'M', 1:'I', 2:'D', 3:'N', 4:'S', 5:'H', 6:'P', 7:'=', 8:'X', 9:'B'}

    samfile = pysam.AlignmentFile(args.input_bam_file, "rb")
    kmer_reads = pysam.AlignmentFile(args.out_bam_file, "wb", template=samfile)

    trie = ahocorasick.Automaton()

    with open(args.Kmer_file) as f:
        for line in f:
            line_split = line.strip().split('\t')
            trie.add_word(line_split[0], line_split[0])

    trie.make_automaton()
    logging.info("finished making automaton")

    for read in samfile.fetch():
        if not read.is_unmapped and not read.is_secondary and read.is_proper_pair and read.is_paired and\
                not read.is_duplicate and not read.is_supplementary and 'N' in read.cigarstring:

            for end, kmer in trie.iter(read.query_sequence.upper()):
                start_index = end - len(kmer) + 1
                end_index = end + 1

                if all(ref_pos is None for ref_pos in
                       read.get_reference_positions(full_length=True)[start_index:end_index]):
                    continue

                quality_string = read.to_string().split('\t')[10][start_index:end_index]

                kmer_read = pysam.AlignedSegment()
                kmer_read.query_name = read.query_name
                kmer_read.query_sequence = read.query_sequence[start_index:end_index].upper()
                kmer_read.flag = read.flag
                kmer_read.reference_id = read.reference_id

                for ref_pos in read.get_reference_positions(full_length=True)[start_index:end_index]:
                    if ref_pos is not None:
                        kmer_read.reference_start = ref_pos
                        break

                kmer_read.mapping_quality = read.mapping_quality

                current_ind = 0
                in_kmer = False
                cigarString_temp = ""

                for operation, count in read.cigartuples:
                    if current_ind >= end_index:
                        break

                    if cigar_map[operation] == 'N' or cigar_map[operation] == 'D':
                        if in_kmer:
                            cigarString_temp += str(count) + cigar_map[operation]

                    elif cigar_map[operation] == 'M' or cigar_map[operation] == 'I' or cigar_map[operation] == 'S':
                        if current_ind + count > start_index:
                            cigarString_temp += str(
                                min(end_index, current_ind + count) - max(current_ind, start_index)) + cigar_map[
                                                    operation]
                            in_kmer = True
                        current_ind += count
                    else:
                        logging.warn('Unexpected cigar op {}'.format(cigar_map[operation]))

                if "S" in cigarString_temp or "N" not in cigarString_temp:
                    continue

                kmer_read.cigarstring = cigarString_temp
                kmer_read.query_qualities = pysam.qualitystring_to_array(quality_string)
                kmer_reads.write(kmer_read)

    kmer_reads.close()
    samfile.close()
    logging.info("Done!")
    stop = timeit.default_timer()
    logging.info("total search time: {}".format(stop - start))


if __name__ == '__main__':
    main()
