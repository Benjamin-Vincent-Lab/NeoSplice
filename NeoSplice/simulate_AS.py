import numpy as np
import os
import argparse
import collections
import itertools
from BCBio import GFF
import logging


class Gene(object):
    def __init__(self, chromosome, startposition, endposition, strand, ID):
        self.chromosome = chromosome
        self.startposition= startposition
        self.endposition = endposition
        self.strand = strand
        self.frame = '.'
        self.transcripts = []
        self.score = '.'
        self.ID = ID
        self.left_splices = set()
        self.right_splices = set()
        self.all_splices = set()


class Transcript(object):
    def __init__(self, chromosome, startposition, endposition, strand, ID, genename, cds_list, start_codon, stop_codon,
                 splice_list, cds_set):
        self.chromosome = chromosome
        self.startposition = startposition
        self.endposition = endposition
        self.strand = strand
        self.ID = ID
        self.genename = genename
        self.cds_list = cds_list
        self.start_codon = start_codon
        self.stop_codon = stop_codon
        self.splice_list = splice_list
        self.score = '.'
        self.frame = '.'
        self.cds_set = cds_set


strand_symbol = {+1: '+', -1: '-'}


def _find_annotated_transcripts(gff_in_file, chromosome):
    in_file = gff_in_file
    in_handle = open(in_file)
    limit_info = dict(gff_id=[chromosome])
    transcripts_dic = {}
    for chromosome_inf in GFF.parse(in_handle, limit_info=limit_info):
        for gene in chromosome_inf.features:
            if "protein_coding" in gene.qualifiers["gene_type"]:
                gene_name = gene.qualifiers["gene_name"][0]
                gene_data = Gene(chromosome, gene.location.start.position, gene.location.end.position, gene.strand,
                                 gene_name)
                transcripts_dic[gene_name] = gene_data
                qualified = False
                for transcript in gene.sub_features:
                    if "protein_coding" in transcript.qualifiers["transcript_type"]:
                        start_codon = []
                        stop_codon = []
                        cds_list = []
                        splice_list = []
                        strand = transcript.strand

                        for feature in transcript.sub_features:
                            if feature.type == "start_codon":
                                start_codon.append(feature)
                            elif feature.type == "stop_codon":
                                stop_codon.append(feature)
                            elif feature.type == "CDS":
                                cds_list.append(feature)

                        cds_list.sort(key=lambda x: x.location.start.position)
                        start_codon.sort(key=lambda x: x.location.start.position)
                        stop_codon.sort(key=lambda x: x.location.start.position)
                        cds_set = set([(cds.location.start.position, cds.location.end.position) for cds in cds_list])
                        i = 0

                        while i < len(cds_list) - 1:
                            splice_list.append(
                                (cds_list[i].location.end.position, cds_list[i + 1].location.start.position))
                            i += 1

                        if len(start_codon) != 0 and len(cds_list) > 3:
                            gene_data.transcripts.append(Transcript(chromosome,
                                                                    transcript.location.start.position,
                                                                    transcript.location.end.position, strand,
                                                                    transcript.qualifiers["transcript_id"][0],
                                                                    gene.qualifiers["gene_name"][0], cds_list,
                                                                    start_codon, stop_codon, splice_list, cds_set))

                            qualified = True
                            gene_data.all_splices |= set(splice_list)
                            gene_data.left_splices |= set([splice[0] for splice in splice_list])
                            gene_data.right_splices |= set([splice[1] for splice in splice_list])

                if not qualified:
                    transcripts_dic.pop(gene_name)

    return transcripts_dic


def get_transcript_occurrence(transcripts, novel_splice):
    count = 0
    for transcript in transcripts:
        if novel_splice[0] in [splice[0] for splice in transcript.splice_list] and novel_splice[1] in [splice[1] for
                                                                                    splice in transcript.splice_list]:
            count += 1
    return count


def get_possible_transcript_occurrence_donor_acceptor(transcripts, splice_pos, is_doner, pos_type, annotated=False):
    count = 0
    for transcript in transcripts:
        annotated_ends = [splice[0] for splice in transcript.splice_list] if is_doner else [splice[1] for splice in
                                                                                            transcript.splice_list]
        if pos_type == "intronic":
            if is_from_intron(transcript, splice_pos, is_doner) or annotated and splice_pos in annotated_ends:
                count += 1
        elif pos_type == "exonic":
            if is_from_exon(transcript, splice_pos, is_doner) or annotated and splice_pos in annotated_ends:
                count += 1
    return count


def is_from_intron(transcript, splice_pos, is_doner):
    for intron in transcript.splice_list:
        if is_doner and intron[0] < splice_pos <= intron[1] - 1:
            return True
        elif not is_doner and intron[0] <= splice_pos - 1 < intron[1] - 1:
            return True
    else:
        return False


def is_from_exon(transcript, splice_pos, is_doner):
    for cds in transcript.cds_list:
        if is_doner and cds.location.start.position < splice_pos <= cds.location.end.position - 1:
            return True
        elif not is_doner and cds.location.start.position <= splice_pos - 1 < cds.location.end.position - 1:
            return True
    else:
        return False


def write_gtf(gene_data, tumor_outf, normal_outf, splice_type=None):
    transcript_left_splices = gene_data.left_splices
    transcript_right_splices = gene_data.right_splices
    if splice_type == "exon_skipping":
        possible_novel_splices = sorted([(a, b, get_transcript_occurrence(gene_data.transcripts, (a, b))) for a, b in
                                  itertools.product(list(transcript_left_splices), list(transcript_right_splices)) if
                                  a + 100 < b and (a, b) not in gene_data.all_splices], key=lambda x: x[2])

        possible_novel_splices = possible_novel_splices[-3:]
        if not possible_novel_splices:
            return
        novel_splice_idx = np.random.choice(len(possible_novel_splices))
        novel_splice = possible_novel_splices[novel_splice_idx][:-1]

    elif splice_type == "exon_loss":
        donor_valid_count = collections.Counter()
        acceptor_valid_count = collections.Counter()
        donor_possible_count = collections.Counter()
        acceptor_possible_count = collections.Counter()
        for splice in gene_data.all_splices:
            donor_valid_count[splice[0]] = get_possible_transcript_occurrence_donor_acceptor(gene_data.transcripts,
                                                                                             splice[0], True, "exonic")
            acceptor_valid_count[splice[1]] = get_possible_transcript_occurrence_donor_acceptor(gene_data.transcripts,
                                                                                        splice[1], False, "exonic")
            donor_possible_count[splice[0]] = get_possible_transcript_occurrence_donor_acceptor(gene_data.transcripts,
                                                                                        splice[0], True, "exonic", True)
            acceptor_possible_count[splice[1]] = get_possible_transcript_occurrence_donor_acceptor(
                                                                gene_data.transcripts, splice[1], False, "exonic", True)

        donor_valid_top3 = [ele[0] for ele in donor_valid_count.most_common(3) if ele[1] != 0]
        acceptor_valid_top3 = [ele[0] for ele in acceptor_valid_count.most_common(3) if ele[1] != 0]
        donor_possible_top3 = [ele[0] for ele in donor_possible_count.most_common(3) if ele[1] != 0]
        acceptor_possible_top3 = [ele[0] for ele in acceptor_possible_count.most_common(3) if ele[1] != 0]

        left_comb = [(a, b) for a, b in itertools.product(donor_valid_top3, acceptor_possible_top3) if
                                                                    a + 100 < b and (a, b) not in gene_data.all_splices]
        right_comb = [(a, b) for a, b in itertools.product(donor_possible_top3, acceptor_valid_top3) if
                                                                    a + 100 < b and (a, b) not in gene_data.all_splices]

        possible_novel_splices = left_comb + right_comb
        if not possible_novel_splices:
            return
        novel_splice_idx = np.random.choice(len(possible_novel_splices))
        novel_splice = possible_novel_splices[novel_splice_idx]

    elif splice_type == "intron_gain":
        donor_valid_count = collections.Counter()
        acceptor_valid_count = collections.Counter()
        donor_possible_count = collections.Counter()
        acceptor_possible_count = collections.Counter()
        for splice in gene_data.all_splices:
            donor_valid_count[splice[0]] = get_possible_transcript_occurrence_donor_acceptor(gene_data.transcripts,
                                                                                            splice[0], True, "intronic")
            acceptor_valid_count[splice[1]] = get_possible_transcript_occurrence_donor_acceptor(gene_data.transcripts,
                                                                                        splice[1], False, "intronic")
            donor_possible_count[splice[0]] = get_possible_transcript_occurrence_donor_acceptor(gene_data.transcripts,
                                                                                    splice[0], True, "intronic", True)
            acceptor_possible_count[splice[1]] = get_possible_transcript_occurrence_donor_acceptor(
                                                            gene_data.transcripts, splice[1], False, "intronic", True)

        donor_valid_top3 = [ele[0] for ele in donor_valid_count.most_common(3) if ele[1] != 0]
        acceptor_valid_top3 = [ele[0] for ele in acceptor_valid_count.most_common(3) if ele[1] != 0]
        donor_possible_top3 = [ele[0] for ele in donor_possible_count.most_common(3) if ele[1] != 0]
        acceptor_possible_top3 = [ele[0] for ele in acceptor_possible_count.most_common(3) if ele[1] != 0]

        left_comb = [(a, b) for a, b in itertools.product(donor_valid_top3, acceptor_possible_top3) if
                     a + 100 < b and (a, b) not in gene_data.all_splices]
        right_comb = [(a, b) for a, b in itertools.product(donor_possible_top3, acceptor_valid_top3) if
                      a + 100 < b and (a, b) not in gene_data.all_splices]

        possible_novel_splices = left_comb + right_comb
        if not possible_novel_splices:
            return
        novel_splice_idx = np.random.choice(len(possible_novel_splices))
        novel_splice = possible_novel_splices[novel_splice_idx]

    else:
        logging.info("type not supported")
        return

    tumor_outf.write('\t'.join(
        [gene_data.chromosome, "test", "gene", str(gene_data.startposition + 1), str(gene_data.endposition),
         gene_data.score, strand_symbol[gene_data.strand],
         gene_data.frame,
         """gene_id "{}";""".format(gene_data.ID) + """ transcript_id "{}";""".format(gene_data.ID)]) + '\n')

    normal_outf.write('\t'.join(
        [gene_data.chromosome, "test", "gene", str(gene_data.startposition + 1), str(gene_data.endposition),
         gene_data.score, strand_symbol[gene_data.strand],
         gene_data.frame,
         """gene_id "{}";""".format(gene_data.ID) + """ transcript_id "{}";""".format(gene_data.ID)]) + '\n')

    logging.info("selected novel splice {}".format(novel_splice))

    for transcript in gene_data.transcripts:
        if novel_splice[0] <= transcript.cds_list[0].location.start.position or novel_splice[1] >= transcript.cds_list[
                -1].location.end.position or novel_splice in transcript.splice_list:
            continue

        if splice_type == "exon_skipping" and not (
                            novel_splice[0] in [splice[0] for splice in transcript.splice_list] and novel_splice[1] in [
                            splice[1] for splice in transcript.splice_list]):
            continue

        if splice_type == "exon_loss":
            donor_annotated = [splice[0] for splice in transcript.splice_list]
            acceptor_annotated = [splice[1] for splice in transcript.splice_list]
            if is_from_exon(transcript, novel_splice[0], True) and (
                    novel_splice[1] in acceptor_annotated or is_from_exon(transcript, novel_splice[1], False)):
                pass
            elif is_from_exon(transcript, novel_splice[1], False) and (
                    novel_splice[0] in donor_annotated or is_from_exon(transcript, novel_splice[0], True)):
                pass
            else:
                continue

        if splice_type == "intron_gain":
            donor_annotated = [splice[0] for splice in transcript.splice_list]
            acceptor_annotated = [splice[1] for splice in transcript.splice_list]
            if is_from_intron(transcript, novel_splice[0], True) and (
                    novel_splice[1] in acceptor_annotated or is_from_intron(transcript, novel_splice[1], False)):
                pass
            elif is_from_intron(transcript, novel_splice[1], False) and (
                    novel_splice[0] in donor_annotated or is_from_intron(transcript, novel_splice[0], True)):
                pass
            else:
                continue

        left_index = None
        right_index = None
        left_type = None
        right_type = None

        tumor_outf.write('\t'.join(
            [transcript.chromosome, "test", "transcript", str(transcript.startposition + 1),
             str(transcript.endposition),
             transcript.score, strand_symbol[transcript.strand], str(transcript.frame),
             """gene_id "{}";""".format(gene_data.ID) + """ transcript_id "{}";""".format(
                 transcript.ID)]) + '\n')

        normal_outf.write('\t'.join(
            [transcript.chromosome, "test", "transcript", str(transcript.startposition + 1),
             str(transcript.endposition),
             transcript.score, strand_symbol[transcript.strand], str(transcript.frame),
             """gene_id "{}";""".format(gene_data.ID) + """ transcript_id "{}";""".format(
                 transcript.ID)]) + '\n')

        for i, exon in enumerate(transcript.cds_list):
            # splice not originate from end of an exist intron
            if exon.location.start.position < novel_splice[0] <= exon.location.end.position - 1:
                left_index = i
                left_type = "exon"
            # splice not end right before start of an exist intron
            if exon.location.start.position <= novel_splice[1] - 1 < exon.location.end.position - 1:
                right_index = i
                right_type = "exon"

        for i, intron in enumerate(transcript.splice_list):
            if intron[0] <= novel_splice[0] <= intron[1] - 1:
                left_index = i
                left_type = "intron"
            if intron[0] <= novel_splice[1] - 1 <= intron[1] - 1:
                right_index = i
                right_type = "intron"

        if right_index is None or left_index is None:
            continue

        skipped_index = range(left_index + 1, right_index)

        i = 0
        temp_list = []
        ref_list = []

        while i < len(transcript.cds_list):

            if i in skipped_index:
                ref_list.append('\t'.join(
                    [transcript.chromosome, "test", "exon",
                     str(transcript.cds_list[i].location.start.position + 1),
                     str(transcript.cds_list[i].location.end.position), '.',
                     strand_symbol[transcript.cds_list[i].strand], '.',
                     """gene_id "{}";""".format(gene_data.ID) + """ transcript_id "{}";""".format(
                         transcript.ID)]) + '\n')

            elif i != left_index and i != right_index:
                temp_list.append('\t'.join(
                    [transcript.chromosome, "test", "exon",
                     str(transcript.cds_list[i].location.start.position + 1),
                     str(transcript.cds_list[i].location.end.position), '.',
                     strand_symbol[transcript.cds_list[i].strand], '.',
                     """gene_id "{}";""".format(gene_data.ID) + """ transcript_id "{}";""".format(
                         transcript.ID)]) + '\n')
                ref_list.append('\t'.join(
                    [transcript.chromosome, "test", "exon",
                     str(transcript.cds_list[i].location.start.position + 1),
                     str(transcript.cds_list[i].location.end.position), '.',
                     strand_symbol[transcript.cds_list[i].strand], '.',
                     """gene_id "{}";""".format(gene_data.ID) + """ transcript_id "{}";""".format(
                         transcript.ID)]) + '\n')
            # left splice right splice in same intron/exon
            elif i == left_index == right_index:
                temp_list.append('\t'.join(
                    [transcript.chromosome, "test", "exon",
                     str(transcript.cds_list[i].location.start.position + 1),
                     str(novel_splice[0]), '.',
                     strand_symbol[transcript.cds_list[i].strand], '.',
                     """gene_id "{}";""".format(gene_data.ID) + """ transcript_id "{}";""".format(
                         transcript.ID)]) + '\n')
                ref_list.append('\t'.join(
                    [transcript.chromosome, "test", "exon",
                     str(transcript.cds_list[i].location.start.position + 1),
                     str(transcript.cds_list[i].location.end.position), '.',
                     strand_symbol[transcript.cds_list[i].strand], '.',
                     """gene_id "{}";""".format(gene_data.ID) + """ transcript_id "{}";""".format(
                         transcript.ID)]) + '\n')
                #connect to next exon
                if right_type == "intron":
                    i += 1
                    ref_list.append('\t'.join(
                        [transcript.chromosome, "test", "exon",
                         str(transcript.cds_list[i].location.start.position + 1),
                         str(transcript.cds_list[i].location.end.position), '.',
                         strand_symbol[transcript.cds_list[i].strand], '.',
                         """gene_id "{}";""".format(gene_data.ID) + """ transcript_id "{}";""".format(
                             transcript.ID)]) + '\n')

                temp_list.append('\t'.join(
                    [transcript.chromosome, "test", "exon",
                     str(novel_splice[1] + 1),
                     str(transcript.cds_list[i].location.end.position), '.',
                     strand_symbol[transcript.cds_list[i].strand], '.',
                     """gene_id "{}";""".format(gene_data.ID) + """ transcript_id "{}";""".format(
                         transcript.ID)]) + '\n')

            elif i == left_index:
                temp_list.append('\t'.join(
                    [transcript.chromosome, "test", "exon",
                     str(transcript.cds_list[i].location.start.position + 1),
                     str(novel_splice[0]), '.',
                     strand_symbol[transcript.cds_list[i].strand], '.',
                     """gene_id "{}";""".format(gene_data.ID) + """ transcript_id "{}";""".format(
                         transcript.ID)]) + '\n')
                ref_list.append('\t'.join(
                    [transcript.chromosome, "test", "exon",
                     str(transcript.cds_list[i].location.start.position + 1),
                     str(transcript.cds_list[i].location.end.position), '.',
                     strand_symbol[transcript.cds_list[i].strand], '.',
                     """gene_id "{}";""".format(gene_data.ID) + """ transcript_id "{}";""".format(
                         transcript.ID)]) + '\n')

            elif i == right_index:
                if right_type == "intron":
                    ref_list.append('\t'.join(
                        [transcript.chromosome, "test", "exon",
                         str(transcript.cds_list[i].location.start.position + 1),
                         str(transcript.cds_list[i].location.end.position), '.',
                         strand_symbol[transcript.cds_list[i].strand], '.',
                         """gene_id "{}";""".format(gene_data.ID) + """ transcript_id "{}";""".format(
                             transcript.ID)]) + '\n')
                    i += 1

                temp_list.append('\t'.join(
                        [transcript.chromosome, "test", "exon",
                         str(novel_splice[1] + 1),
                         str(transcript.cds_list[i].location.end.position), '.',
                         strand_symbol[transcript.cds_list[i].strand], '.',
                         """gene_id "{}";""".format(gene_data.ID) + """ transcript_id "{}";""".format(
                             transcript.ID)]) + '\n')
                ref_list.append('\t'.join(
                        [transcript.chromosome, "test", "exon",
                         str(transcript.cds_list[i].location.start.position + 1),
                         str(transcript.cds_list[i].location.end.position), '.',
                         strand_symbol[transcript.cds_list[i].strand], '.',
                         """gene_id "{}";""".format(gene_data.ID) + """ transcript_id "{}";""".format(
                             transcript.ID)]) + '\n')

            i += 1

        if transcript.strand == -1:
            for rec in temp_list[::-1]:
                tumor_outf.write(rec)
            for rec in ref_list[::-1]:
                normal_outf.write(rec)
        else:
            for rec in temp_list:
                tumor_outf.write(rec)
            for rec in ref_list:
                normal_outf.write(rec)


def main():
    logging.basicConfig(level=logging.DEBUG, format='%(asctime)-15s [%(processName)s.%(levelname)s] %(message)s')
    parser = argparse.ArgumentParser(description='Utilities for simulating reads')
    parser.add_argument('sample_id', type=str, nargs='?', help='provide sample id here')
    parser.add_argument('gff_file', type=str, nargs='?', help='provide gff3 file path here')
    parser.add_argument('outdir', type=str, nargs='?', help='provide outdir here')
    parser.add_argument("num_alterations", type=int, nargs='?', help='provide number of altered genes per chromosome')
    parser.add_argument("simulation_type", type=str, nargs='?', help='provide simulation type')
    args = parser.parse_args()

    gff_in_file = args.gff_file
    alteration_number_lim = args.num_alterations
    tumor_GTF_output = args.outdir + "tumor_GTF_output.gtf"
    normal_GTF_output = args.outdir + "normal_GTF_output.gtf"

    tumor_outf = open(tumor_GTF_output, 'w')
    normal_outf = open(normal_GTF_output, 'w')

    annotatedTransctipts = {}

    for i in range(1, 23):
        if args.simulation_type == "random":
            splice_type = np.random.choice(["exon_skipping", "exon_loss", "intron_gain"])
        else:
            splice_type = args.simulation_type
        chromosome = "chr" + str(i)
        logging.info("processing {}".format(chromosome))
        annotatedTransctipts[chromosome] = _find_annotated_transcripts(gff_in_file, chromosome)
        keys_to_choose = [key for key in annotatedTransctipts[chromosome] if
                          len(annotatedTransctipts[chromosome][key].transcripts) >= 8]
        genes_to_alter = np.random.choice(keys_to_choose, alteration_number_lim, replace=False)
        logging.info("selected genes: {}".format(','.join(genes_to_alter)))

        for gene in genes_to_alter:
            write_gtf(annotatedTransctipts[chromosome][gene], tumor_outf, normal_outf, splice_type)

    tumor_outf.close()
    normal_outf.close()


if __name__ == '__main__':
    main()

