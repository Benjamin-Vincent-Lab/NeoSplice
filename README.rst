""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
NeoSplice: A bioinformatics method for prediction of splice variant neoantigens
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

NeoSplice is a novel computational method for splice variant neoantigen prediction based on 1) prediction of tumor-specific k-mers from RNA-seq data, 2) alignment of differentially expressed k-mers to the splice graph, and 3) inference of the variant transcript with MHC binding prediction.  NeoSplice demonstrates high sensitivity and precision (>80% on average across all splice variant classes) through in silico simulated RNA-seq data.  Additionally, predicted splice variant neoantigens in the K562.A2 cell line have been validated using mass spectrometry immunopeptidome analysis.  NeoSplice provides a well-validated platform for prediction of TSA vaccine targets for future cancer antigen vaccine studies to evaluate the clinical efficacy of splice variant neoantigens.

The NeoSplice workflow is currently set up for prediction of MHC-I antigens, with future plans to additionally expand to MHC-II prediction.

NeoSplice is free for academic and non-profit use.

============
Installation
============

-------------------
System requirements
-------------------
- multi-string BWT (MSBWT): https://github.com/holtjma/msbwt
- MSBWT-IS: https://github.com/holtjma/msbwt-is
- NetMHCpan 4.0

------------
Dependencies
------------
- Python 2.7
- networkx 1.11
- pyahocorasick 1.4.0
- bcbio-gff 0.6.4
- pyfaidx 0.5.3.1
- pysam 0.14.1
- biopython 1.70
- scipy 1.2.0

------------------
Installation steps
------------------

Using a Python 2.7 VirtualEnv:
::
    git clone https://github.com/Benjamin-Vincent-Lab/NeoSplice.git
    cd NeoSplice
    virtualenv .
    source bin/activate
    pip install -r requirements.txt

Using Singularity with Docker image:
::
    singularity pull docker://benjaminvincentlab/neosplice:0.0.2
    singularity shell -B </path/with/input_bams_gffs_fa_etc> neosplice_0.0.2.sif

========
Workflow
========

-----------------
Summary of steps
-----------------
Multiple steps are needed to identify a novel splice that occurs specifically in tumor cell transcripts whose translation will result in a neopeptide that can be targeted by T cells.  Broadly, these steps are: 

- **Step 1:** Tumor-specific k-mer generation. Using one RNA-seq dataset T of tumor cells and one RNA-seq dataset N of normal cells, tumor-specific k-mer sequences present abundantly in the transcriptome of the tumor cell, but not/rarely expressed in the normal cells are identified.

- **Step 2:** Prediction of splice variant transcripts.  The splice graph G from the tumor cell RNA-seq data is built, and tumor specific k-mers from above are mapped to novel splice variant transcripts.  Gencode annotations are used to determine whether the novel splice lies within a protein coding region and infer the reading frame of the transcript. 

- **Step 3:** Prediction of splice variant neoantigens. Novel splice junctions contained within each splice variant transcript are translated in the inferred open reading frame.  MHC binding affinity prediction is performed on translated peptide sequences to determine which novel regions may yield a neopeptide.

.. image:: https://github.com/Benjamin-Vincent-Lab/NeoSplice/blob/master/images/Neosplice_fig1.PNG


-------------------
Summary of workflow
-------------------

Functionally, these above steps are accomplished by individual Python2 scripts, alongside the prior listed dependencies.  This workflow is summarized in the below figure:

.. image:: https://github.com/Benjamin-Vincent-Lab/NeoSplice/blob/master/images/Neosplice_workflow.jpg

This workflow is summarized step-by-step below. Additionally, an example **Nextflow** script is provided in the ``./Nextflow_example`` directory of the GitHub repo, which provides the entire workflow as an .nf script.

0. Input files
----------------------------
The following input files will be referenced in the below workflow steps:

- **tumor.bam**: Aligned RNA-seq file for tumor sample of interest using a splice-aware aligner.  Currently, we recommend the STAR aligner: https://github.com/alexdobin/STAR
- **normal.bam**: Aligned RNA-seq file for matched-normal sample of interest using a splice-aware aligner.  Currently, we recommend the STAR aligner: https://github.com/alexdobin/STAR
- **reference.fa**: Reference genome fasta file, preferably the same file used to generate the genome index for STAR.
- **reference.gff**: Reference gff3 file, preferably the same build as the .gtf file that was used to generate the genome index for STAR.
- **./Reference_peptidome**: A GRCh38 reference peptidome of 8-11mer peptides is contained within the repository at ``./Reference_peptidome``. If necessary, you can use the script  ``./NeoSplice/generate_reference_peptidome.py`` to generate a different reference build using matching fasta and gff files. Prior to use (step 9), unzip the files contained in this directory using the following command:

.. code-block:: python
    gunzip ./Reference_peptidome/*

1. augmented_splice_graph.py
----------------------------
This step builds the splice graph for the tumor, with ``augmented_splice_graph.py`` run for each individual chromosome of interest.  The output for each instance (i.e. chromosome) is a ``.json`` file.  There are several arguments included for this step: **p-error**, **cutoff**, **min-coverage**, and **min-variants**.  While we cannot provide optimal argument recommendations for every sample, below are the values used for simulated read data benchmarking and mass spectrometry validated K562.A2 cell line splice variant neoantigens.  Below is an example for chromosome 1:

.. code-block:: python

    mkdir ./tumor1_splice_graph
    python /NeoSplice/augmented_splice_graph.py build \
        --bam ./path/to/tumor.bam \
        --seq chr1 \
        --genome ./path/to/reference.fa \
        --min-variants 10 \
        --cutoff 0.000005 \
        --gff  ./path/to/reference.gff \
        --out ./tumor_splice_graph

2. convert_bam_to_fasta.py
----------------------------
This step is a simple script to back-convert the STAR-aligned **tumor.bam** and **normal.bam** files back into fasta format:

.. code-block:: python

    python /NeoSplice/convert_bam_to_fasta.py \
        --bam_file tumor.bam \
        --R1_out tumor_R1.fasta \
        --R2_out tumor_R2.fasta
    python /NeoSplice/convert_bam_to_fasta.py \
        -bam_file normal.bam \
        --R1_out normal_R1.fasta \
        --R2_out normal_R2.fasta

3. Run multi-string BWT
----------------------------
This step uses the MSBWT-IS tool developed by Holt and colleagues (https://github.com/holtjma/msbwt-is), followed by a bash script to convert the output format for downstream compatibility:

.. code-block:: python
     
    mkdir ./tumor_bwt/
    mkdir ./normal_bwt/
    mkdir ./tumor_bwt_temp/
    mkdir ./normal_bwt_temp/
    ./msbwt-is/msbwtis tumor_bwt_temp/ tumor_R1.fasta tumor_R2.fasta
    ./msbwt-is/msbwtis normal_bwt_temp/ normal_R1.fasta normal_R2.fasta
    bash ./NeoSplice/convert_BWT_format.bash ./tumor_bwt_temp ./tumor_bwt 
    bash ./NeoSplice/convert_BWT_format.bash ./normal_bwt_temp ./normal_bwt

4. get_max_kmer_length.py
----------------------------
This step searches for the maximum read length contained within either the tumor or matched-normal files, returning an output value for use in step 5.  If you know this value already, this step can be skipped:

.. code-block:: python

     python /NeoSplice/get_max_kmer_length.py \
         --tumor_bam tumor.bam \
         --normal_bam normal.bam

5. Kmer_search_bwt.py
----------------------------
This step uses the MSBWTs generated in step 3 and searches for differentially expressed Kmers between tumor and matched-normal samples.  There are two argument variables that can be adjusted here -- **Tmin** (minimum expression of a given Kmer in the tumor) and **Nmax** (maximum expression of a given Kmer in the normal).  For a Kmer to be considered differentially expressed, it must be > **Tmin** AND < **Nmax**.  Typically, you may consider setting **Tmin** to 20-35 and **Nmax** to 1-4.  The **max_length** argument should be set to the value obtained from **step 4**, or the maximum read length of the input files.

 .. code-block:: python

    mkdir .tumor_kmers
    python ./NeoSplice/Kmer_search_bwt.py \
        --tumor_bwt = ./tumor_bwt/ \
        --normal_bwt ./normal_bwt/ \
        --processors 1 \
        --max_length $read_length \
        --tumor_threshold 20 \
        --normal_threshold 4  \
        --outdir ./tumor_kmers/
    cat ./tumor_kmers/Tumor_kmers_* >  ./tumor_kmers/merged_Tumor_kmers.txt

6. search_bam.py and Samtools sort/index
----------------------------------------
This step uses an Aho–Corasick algorithm (pyahocorasick 1.4.0) to search for the reads that contain tumor specific Kmers in the tumor RNA-seq BAM file.  This method runs in time linear in the size of the BAM file.  For each occurrence, the Kmer-containing portion of the read along with corresponding quality scores and Cigar strings is written to a new BAM file.  This output BAM is then sorted and indexed using Samtools.

 .. code-block:: python

    python ./NeoSplice/search_bam.py \
        --Kmer_file ./tumor_kmers/merged_Tumor_kmers.txt \
        --input_bam_file tumor.bam \
        --out_bam_file tumor_Kmer.bam 
    samtools sort -m 15G -o tumor_Kmer_sorted.bam tumor_Kmer.bam
    samtools index tumor_Kmer_sorted.bam

7. get_splice_junctions.py
----------------------------------------
This step collects a list of all splice junctions from the tumor and normal BAM files, storing these in a text file for downstream use.

 .. code-block:: python

    python /NeoSplice/get_splice_junctions.py \
        --input_bam tumor.bam \
        --out_file tumor_junctions.txt
    python /NeoSplice/get_splice_junctions.py \
        --input_bam normal.bam \
        --out_file normal_junctions.txt

8. kmer_graph_inference.py
----------------------------------------
In this step, each splice variant transcript sequence is identified by depth-first search.  This is then concatenated with the tumor specific Kmer sequence and translated into 8-11mer peptides for MHC-I neoantigen prediction.  Binding affinity to MHC molecules expressed by the tumor for in-silico generated peptides is predicted using NetMHCpan-4.0.  Arguments to consider in this step include **HLA_I** (provide list of NetMHCpan-compatible alleles for antigen prediction), as well as **transcript_min_coverage** (the minimum Kmer coverage necessary for a transcript to be considered).  This command is run for each chromosome of interest, with an example for chromsome 1 shown below:

 .. code-block:: python

    python /NeoSplice/kmer_graph_inference.py \
        --sample tumor \
        --chromosome chr1 \
        --bam_file tumor.bam \
        --gff_file reference.gff \
        --genome_fasta reference.fasta \
        --kmer_bam tumor_Kmer_sorted.bam \
        --splice_graph ./tumor_splice_graph \
        --tumor_junction_file tumor_junctions.txt \
        --normal_junction_file normal_junctions.txt \
        --transcript_min_coverage 15 \
        --HLA_I ${HLA_i} \
        --netMHCpan_path ./netMHCpan-4.0-docker/netMHCpan \
        --outdir ./tumor_output_dir
        
9. SV_summarization.py
----------------------------------------
In this final step, predicted splice variant peptides from above are filtered against the reference peptidome, filtered to peptides with predicted binding affinity >500nM by NetMHCpan-4.0, and summarized into a single output file.  The **data_dir** argument should point to the working directory, one level above the ``outdir`` argument from step 8 (``kmer_graph_inference.py``).  The output from this step provides a summarized text file containing all predicted splice variant neoantigens.

 .. code-block:: python

    python /NeoSplice/SV_summarization.py \
        --ref_dir ./Reference_peptidome \
        --data_dir . \
