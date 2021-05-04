NeoSplice 
===
**Introduction**
---
NeoSplice is a bioinformatics algorithm for prediction of splice variant neoantigens using tumor 
and matched normal RNA-seq data.

**System Requirements**
---
1. multi-string BWT (MSBWT) https://github.com/holtjma/msbwt
2. MSBWT-IS https://github.com/holtjma/msbwt-is
3. NetMHCpan 4.0
4. NetMHCIIpan 3.2

python packages:

* networkx 1.11
* pyahocorasick 1.4.0
* bcbio-gff 0.6.4
* pyfaidx 0.5.3.1
* pysam 0.14.1
* biopython 1.70
* scipy 1.2.0


**Usage**
---

Step1: Build splice graph using augmented_splice_graph.py 
from input RNA-seq bam file

Step2: 
* Convert RNA-seq bam file to FASTA file using
convert_bam_to_fasta.py
* Build MSBWT for tumor and normal RNA-seq data 
using MSBWT-IS
* Search for tumor specific K-mer with Kmer_search_bwt.py
* Obtain output bam file with K-mer containing reads using search_bam.py

Step3: Run augmented_splice_graph.py using splice graphs generated in step 1 and 
aligned K-mer bam file generated in step 2 to predict splice variant neoantigens