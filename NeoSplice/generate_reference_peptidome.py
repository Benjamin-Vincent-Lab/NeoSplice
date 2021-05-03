import os
from pyfaidx import Fasta
from BCBio import GFF
import argparse
import logging


class Transcripts(object):
    def __init__(self, chrnum, startposition, endposition, strand, ID, genename, cds_list, start_codon, stop_codon):
        self.chrnum = chrnum
        self.startposition = startposition
        self.endposition = endposition
        self.strand = strand
        self.ID = ID
        self.genename = genename
        self.cds_list = cds_list
        self.start_codon = start_codon
        self.stop_codon = stop_codon


class Polypeptide(object):
    def __init__(self, peptide_sequence, DNA_sequence, transcript):
        self.peptide_sequence = peptide_sequence
        self.DNA_sequence = DNA_sequence
        self.transcript = transcript


def reverse_complement(seq):
    return u''.join(reversed([complement(base) for base in seq]))


def complement(base):
    if base.upper() == 'A':
        return 'T'
    elif base.upper() == 'T':
        return 'A'
    elif base.upper() == 'G':
        return 'C'
    elif base.upper() == 'C':
        return 'G'
    elif base.upper() == 'N':
        return 'N'
    else:
        raise AssertionError("Unexpected value of nucleotide {}".format(base.upper()))


def find_sequence(chrnum, start, end, genome):
    return genome[chrnum][start:end].seq.upper()


#translate a sequence to pepetides
def translate(seq, transcript, num):
    peptide_list = list()
    gencode = {
      'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
      'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
      'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
      'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
      'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
      'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
      'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
      'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
      'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
      'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
      'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
      'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
      'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
      'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
      'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
      'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W'}

    for i in range(0, len(seq) - 3 * num + 1, 3):
        frameseq = seq[i:i+3*num]
        if 'N' not in frameseq:
            peptide = ''.join([gencode.get(frameseq[3*j:3*j+3], 'X') for j in range(len(frameseq)//3)])
            if '_' in peptide:
                logging.info("wrong peptide: {}".format(peptide))
                break
            peptide_list.append(Polypeptide(peptide, frameseq, transcript))

    return peptide_list


def find_annotated_transcripts(gff_in_file):
    in_file = gff_in_file
    in_handle = open(in_file)
    transcripts_dic = {}
    for chromosome_inf in GFF.parse(in_handle):
        for gene in chromosome_inf.features:
            if "protein_coding" in gene.qualifiers["gene_type"]:
                for transcript in gene.sub_features:
                    if "protein_coding" in transcript.qualifiers["transcript_type"]:
                        start_codon = []
                        stop_codon = []
                        cds_list = []
                        strand = "+" if transcript.strand == +1 else "-"
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

                        if len(start_codon) != 0:
                            transcripts_dic[transcript.id] = Transcripts(chromosome_inf.id,
                                                                         transcript.location.start.position,
                                                                         transcript.location.end.position, strand,
                                                                         transcript.qualifiers['transcript_id'],
                                                                         gene.qualifiers['gene_name'], cds_list,
                                                                         start_codon, stop_codon)
    return transcripts_dic


def output_peptide(path, num, transcript_dic, genome):
    seen = set()

    logging.info("total number of transcripts {}".format(len(transcript_dic)))

    outf = open(path + "reference_peptidome_{}.txt".format(num), 'w')

    for transcript in transcript_dic:
        final_sequence_list = []
        if transcript_dic[transcript].strand == '+':
            for exon in transcript_dic[transcript].cds_list:
                final_sequence_list.append(
                    find_sequence(transcript_dic[transcript].chrnum, exon.location.start.position,
                                  exon.location.end.position, genome))
        elif transcript_dic[transcript].strand == '-':
            for exon in reversed(transcript_dic[transcript].cds_list):
                final_sequence_list.append(reverse_complement(
                    find_sequence(transcript_dic[transcript].chrnum, exon.location.start.position,
                                  exon.location.end.position, genome)))

        final_seq = ''.join(final_sequence_list)

        if len(final_seq) >= 3 * num:
            peptides = translate(final_seq, transcript_dic[transcript], num)
            for i in range(0, len(peptides)):
                if peptides[i].peptide_sequence not in seen:
                    outf.write(peptides[i].peptide_sequence + '\t' + peptides[i].DNA_sequence + '\t' + transcript + '\n')
                    seen.add(peptides[i].peptide_sequence)
    outf.close()


def main():
    logging.basicConfig(level=logging.DEBUG, format='%(asctime)-15s [%(processName)s.%(levelname)s] %(message)s')
    parser = argparse.ArgumentParser(description="Utility for getting reference peptidome")
    parser.add_argument("gff_file", type=str, nargs='?', help="provide the GFF file path")
    parser.add_argument("genome", type=str, nargs='?', help="provide the reference FASTA file path")
    parser.add_argument("out_dir", type=str, nargs='?', help="provide output path here")
    args = parser.parse_args()

    genome = Fasta(args.genome)
    gff_file = os.path.abspath(args.gff_file)
    transcript_dic = find_annotated_transcripts(gff_file)
    data_store_directory = os.path.join(args.out_dir)
    path = os.path.join(data_store_directory + "peptidome_result/")
    if not os.path.isdir(path) and not os.path.exists(path):
        os.makedirs(path, 0777)

    for num in [8, 9, 10, 11, 15]:
        output_peptide(path, num, transcript_dic, genome)


if __name__ == '__main__':
    main()
