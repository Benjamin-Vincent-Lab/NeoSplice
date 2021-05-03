import sys
egg_path='/nas/longleaf/home/shengjie/MSBWT/msbwt-0.3.0-py2.7-linux-x86_64.egg'
sys.path.append(egg_path)
from MUSCython import MultiStringBWTCython as MSBWT
import multiprocessing as mp
import argparse
import timeit
import os
import itertools
import math
import logging


def find_Kmer(Kmer):
    outf = open(outdir + 'Tumor_kmers_{}.txt'.format(Kmer), 'w')

    msbwt_tumor = MSBWT.loadBWT(args.tumor_bwt_1)
    msbwt_normal = MSBWT.loadBWT(args.normal_bwt_1)

    logging.info("finished loading BWTs")

    tLow, tHigh = msbwt_tumor.findIndicesOfStr(Kmer)
    nLow, nHigh = msbwt_normal.findIndicesOfStr(Kmer)

    def Kmer_count(tLow, tHigh, nLow, nHigh, Kmer, outf):
        tLow, tHigh = msbwt_tumor.findIndicesOfStr(Kmer[0],(tLow, tHigh))
        nLow, nHigh = msbwt_normal.findIndicesOfStr(Kmer[0],(nLow, nHigh))

        tumor_count = tHigh - tLow 
        normal_count = nHigh - nLow 

        if tumor_count > tumor_threshold and normal_count < normal_threshold:
            outf.write(Kmer + '\t' + str(tumor_count) + '\t' + str(normal_count) + '\t' + Kmer + '\n')
            return

        elif tumor_count <= tumor_threshold or len(Kmer) == read_length:
            return

        for nucleotide in nucleotide_list:
            Kmer = nucleotide + Kmer
            Kmer_count(tLow, tHigh, nLow, nHigh, Kmer, outf)
            Kmer = Kmer[1:]

    for nucleotide in nucleotide_list:
        Kmer = nucleotide + Kmer
        Kmer_count(tLow, tHigh, nLow, nHigh, Kmer, outf)
        Kmer = Kmer[1:]

    outf.close()


logging.basicConfig(level=logging.DEBUG, format="%(asctime)-15s [%(processName)s.%(levelname)s] %(message)s")
parser = argparse.ArgumentParser(description='Utility for finding tumor specific Kmers with msbwt.')
parser.add_argument('tumor_bwt_1', type=str, nargs='?', help='provide tumor bwt path here')
parser.add_argument('normal_bwt_1', type=str, nargs='?', help='provide normal bwt path here')
parser.add_argument('processors', type=str, nargs='?', help='provide number of processors to be used')
parser.add_argument('max_length', type=str, nargs='?', help='provide maximum k-mer length')
parser.add_argument('tumor_threshold', type=str, nargs='?', help='provide tumor threshold')
parser.add_argument('normal_threshold', type=str, nargs='?', help='provide normal threshold')
parser.add_argument('outdir', type=str, nargs='?', help='provide output directory')
args = parser.parse_args()

start = timeit.default_timer()

outdir = args.outdir

if not os.path.isdir(outdir) and not os.path.exists(outdir):
    os.makedirs(outdir)

read_length = int(args.max_length)
tumor_threshold = int(args.tumor_threshold)
normal_threshold = int(args.normal_threshold)
nucleotide_list = ['T', 'C', 'G', 'A']
processes = []

pattern_suffixes = list(itertools.product(nucleotide_list, repeat=int(math.log(int(args.processors), 4))))

logging.info("start suffixes: {}".format(pattern_suffixes))

if len(list(pattern_suffixes)) == 0:
    find_Kmer('')
else:
    for suffix in pattern_suffixes:
        Kmer = ''.join(suffix)
        p = mp.Process(target=find_Kmer,
                       args=(Kmer,))
        processes.append(p)

    for p in processes:
        p.start()
    for p in processes:
        p.join()

stop = timeit.default_timer()
logging.info("total search time: {}".format(stop - start))
