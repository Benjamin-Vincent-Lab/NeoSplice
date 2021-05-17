import os
import sys
import subprocess
import itertools
import collections
import re
import argparse

#reference_directory = os.path.join('/datastore/nextgenout5/share/labs/Vincent_Lab/ref/ref_peptidome/peptidome_results_hg38/')
#data_store_directory = os.path.join('/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/AML_cell_lines/splice/neosplice_hg38_abra_5_0_2_0.000005_ensemble_avgSJ_final/')

class Sample(object):
    def __init__(self,sample_barcode):
        self.sample_barcode = sample_barcode
        self.Mutation_list = list()
        self.strong_binders = 0
        self.moderate_binders = 0
        self.weak_binders = 0
        self.total_binders = 0
        self.peptides = collections.defaultdict(list)
#        self.peptides_class_II = collections.defaultdict(list)
        self.group = ''
        self.header = list()


def peptide_prediction(num, chromosome, sample, path, ref_peptidome):
    #global header
    with open(path +'%s_peptide_%s_%s.xls'%(sample.sample_barcode, chromosome, num),"r") as f1:
        dat_MHC = [row for row in f1]

    with open(path +'%s_outcome_peptide_%s_%s.txt'%(sample.sample_barcode, chromosome, num), "r") as f:
        dat = [row for row in f]

        n = 1
        if (len(sample.header) == 0):
            sample.header = dat[0]
        elif (len(dat[0].strip().split('\t')) > len(sample.header.strip().split('\t'))):
            sample.header = dat[0]
        for line in dat[1:]:
            line_split = line.strip().split('\t')
            #print line_split[0]
            #print dat_MHC[n+1].strip().split('\t')[1]
            if line_split[0] != dat_MHC[n+1].strip().split('\t')[1] or round(float(line_split[9]),2) != round(float(dat_MHC[n+1].strip().split('\t')[6]),2):
                print "error"
                print chromosome
                print n
                break
            if ('strong_binder' in line_split or 'moderate_binder' in line_split or 'weak_binder' in line_split) and line_split[0] not in ref_peptidome:
                sample.peptides[line_split[0]].append((line, min([int(a[3]) for a in eval(line_split[3])])))
            n += 1



def main():
    #logging.basicConfig(level=logging.DEBUG, format='%(asctime)-15s [%(processName)s.%(levelname)s] %(message)s')
    parser = argparse.ArgumentParser(description="Utility for filtering peptides against the reference peptidome and summarizing output to a single text file.")
    parser.add_argument("--ref_dir", required=True, type=str, nargs='?', help="provide the file path to the reference peptidome")
    parser.add_argument("--data_dir", required=True, type=str, nargs='?', help="provide the path to the outdir from Step 8 ('kmer_graph_inference.py')")
    args = parser.parse_args()

    reference_directory = args.ref_dir
    data_store_directory = args.data_dir

    Sample_dic = collections.OrderedDict()
    samples_to_run = [item for item in os.listdir(data_store_directory)]

    for sample_id in samples_to_run:
        Sample_dic[sample_id] = Sample(sample_id)

    ref_peptidome = set()
    for num in [8,9,10,11]:
        with open(os.path.join(reference_directory + "reference_peptidome_{}.txt".format(num)), 'r') as f:
            for line in f:
                line_split = line.strip().split('\t')
                ref_peptidome.add(line_split[0])


    for sample in Sample_dic.values():
        print sample.sample_barcode
        try:
            for chromosome in range(1, 23):
                chromosome = "chr" + str(chromosome)
                for num in [8,9,10,11]:
                    path = os.path.join(data_store_directory + sample.sample_barcode + '/')
                    #print path
                    peptide_prediction(num, chromosome, sample, path, ref_peptidome)
            #print sample.sample_barcode
            #print 'the order of sample is '+ str(Sample_dic.values().index(sample))

        except IndexError:
            print 'wrong for %s'%(sample.sample_barcode)
            continue
        except IOError:
            print 'incomplete files removed for %s'%(sample.sample_barcode)
            continue

    #for sample in Sample_dic.values():
        outf  = open(os.path.join(data_store_directory + 'SVAgI_{}.txt'.format(sample.sample_barcode)),'w')
        outf.write(sample.header.strip() + '\t' + 'min_expression' + '\n')
        for pep in sample.peptides:
            if pep != 0:
            #if len(sample.peptides[pep]) > 1:
                #print sample.peptides[pep]
                #print pep
                sample.peptides[pep].sort(key=lambda x: x[1])
                outf.write(sample.peptides[pep][-1][0].strip() + '\t' + str(sample.peptides[pep][-1][1]) + '\n')

        outf.close()

    #outf  = open('neoantigen_peptides_class_II_new.txt','w')
    #for sample in Sample_dic.values():
    #    for pep in sample.peptides_class_II:
    #        outf.write(sample.peptides_class_II[pep][-1])

    #outf.close()

if __name__ == '__main__':
    main()
