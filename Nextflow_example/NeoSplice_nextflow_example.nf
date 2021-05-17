// Skeleton nextflow workflow for NeoSplice

params.inventory = '/NeoSplice/Nextflow/pairs_info.txt'
params.chromosomes = '/NeoSplice/Nextflow/chroms.txt'
params.reference_fa = '/path/to/ref/GRCh38_v35/GRCh38.p13.genome.fa'
params.gff = '/path/to/ref/GRCh38_v35/gencode.v35.annotation.gff3'
params.outputdir = '/NeoSplice/Nextflow/run_example/final'
params.intermediate_dir_step1 = '/NeoSplice/Nextflow/run_example/step1'
params.intermediate_dir_step2 = '/NeoSplice/Nextflow/run_example/step2'
params.intermediate_dir_step3 = '/NeoSplice/Nextflow/run_example/step3'
inventory = file(params.inventory)
chromosomes = file(params.chromosomes)
reference_fa = file(params.reference_fa)
gff = file(params.gff)

// Read in BAM file paths
Channel
    .fromPath(inventory)
    .splitCsv(header:true, sep:'\t')
    .map{ row-> tuple(row.sample, file(row.normal), file(row.normal_idx), file(row.tumor), file(row.tumor_idx), row.HLA) }
    .set { bam_channel }

bam_channel.into { bam_channel1; bam_channel2 }

// Get list of chromosomes from 'chroms.txt' file
Channel
    .fromPath(chromosomes)
    .splitCsv(header:true, sep:'\t')
    .map{ row-> tuple(row.chromosome) }
    .set { chrom_channel }

chrom_channel.into { chrom_channel1; chrom_channel2 }

// Cartesian product of bam inventory and chromosomes
bam_chromosomes = bam_channel1.combine(chrom_channel1)

// Process per sample / chromosome
// Registers splice_graph directory as output
process build_splice_graph {
    cpus 1
    memory '12G'
    tag "build_splice_graph_${sample}_${chromosome}"
    publishDir "${params.intermediate_dir_step1}/${sample}", mode: "copy"
//    machineType 'n1-highmem-4'

    input:
    set val(sample), file(normal), file(normal_idx), file(tumor), file(tumor_idx), val(HLA), val(chromosome) from bam_chromosomes
    file reference_fa
    file gff

    output:
    set val(sample), file("${chromosome}_graph.json"), val(chromosome) into splice_graph_files
    """
    mkdir ${sample}_${chromosome}_splice_graph
    python /NeoSplice/augmented_splice_graph.py build --bam ${tumor} --seq ${chromosome} --genome ${reference_fa} --min-variants 10 --cutoff 0.000005 --gff  ${gff} --out ${sample}_${chromosome}_splice_graph

    cp ${sample}_${chromosome}_splice_graph/${chromosome}_graph.json ${chromosome}_graph.json
    """
}

// Process per sample
process tumor_specific_kmers {
    cpus 1
    memory '40G'
    tag "tumor_specific_kmers_${sample}"
    publishDir "${params.intermediate_dir_step2}/${sample}", mode: "copy", pattern: "*_{Kmer_sorted.bam*,bwt,junctions.txt}"
//    machineType 'n1-highmem-8'

    input:
    set val(sample), file(normal), file(normal_idx), file(tumor), file(tumor_idx), val(HLA) from bam_channel2
    file reference_fa

    output:
    set val(sample), file(tumor), file(tumor_idx) , file("${sample}_Kmer_sorted.bam"), file("${sample}_Kmer_sorted.bam.bai"), file("${sample}_tumor_junctions.txt"), file("${sample}_normal_junctions.txt"), file("${sample}_tumor_bwt"), file("${sample}_normal_bwt"), val(HLA) into tumor_specific_kmers_files

    """
    read_length=\$(python /NeoSplice/get_max_kmer_length.py --tumor bam ${tumor} --normal_bam ${normal})

    mkdir ${sample}_tumor_bwt/
    mkdir ${sample}_normal_bwt/
    mkdir ${sample}_tumor_bwt_temp/
    mkdir ${sample}_normal_bwt_temp/
    mkdir ${sample}_kmers/
    python /NeoSplice/convert_bam_to_fasta.py --bam_file ${tumor} --R1_out ${sample}_tumor_R1.fasta --R2_out ${sample}_tumor_R2.fasta
    python /NeoSplice/convert_bam_to_fasta.py --bam_file ${normal} --R1_out ${sample}_normal_R1.fasta --R2_out ${sample}_normal_R2.fasta
    /msbwt-is/msbwtis ${sample}_tumor_bwt_temp/ ${sample}_tumor_R1.fasta ${sample}_tumor_R2.fasta
    /msbwt-is/msbwtis ${sample}_normal_bwt_temp/ ${sample}_normal_R1.fasta ${sample}_normal_R2.fasta
    bash /NeoSplice/convert_BWT_format.bash ${sample}_tumor_bwt_temp ${sample}_tumor_bwt/ 
    bash /NeoSplice/convert_BWT_format.bash ${sample}_normal_bwt_temp ${sample}_normal_bwt/
    python /NeoSplice/Kmer_search_bwt.py --tumor_bwt ${sample}_tumor_bwt/ --normal_bwt ${sample}_normal_bwt/ --processors 1 --max_length \${read_length} --tumor_threshold 20 --normal_threshold 4  --outdir ${sample}_kmers/
    cat ${sample}_kmers/Tumor_kmers_* >  ${sample}_kmers/merged_Tumor_kmers.txt
    python /NeoSplice/search_bam.py --Kmer_file ${sample}_kmers/merged_Tumor_kmers.txt --input_bam_file ${tumor} --out_bam_file ${sample}_Kmer.bam
	samtools sort -m 15G -o ${sample}_Kmer_sorted.bam ${sample}_Kmer.bam
    samtools index ${sample}_Kmer_sorted.bam
    python /NeoSplice/get_splice_junctions.py --input_bam ${tumor} --out_file ${sample}_tumor_junctions.txt
    python /NeoSplice/get_splice_junctions.py --input_bam ${normal} --out_file ${sample}_normal_junctions.txt
    """
}

// Join output of tumor_specific_kmers with output of build_splice_graph using the sample id (cross)
// Flatten the 2 lists into the individual elements (flatten)
// Group together the 13 input variables  (collate)
// Runs once per sample / chromosome combination
process splice_variants {
    cpus 1
    memory '5G'
    tag "splice_variants_${sample}_${chromosome}"
    publishDir "${params.outputdir}/${sample}", mode: "copy"
//    machineType 'n1-highmem-4'

    input:
    set val(sample), file(tumor), file(tumor_idx), file(tumor_specific_kmers), file(bam_index), file(tumor_junctions), file(normal_junctions), file(tumor_bwt), file(normal_bwt), val(HLA), val(sample2), file(splice_graph), val(chromosome) from tumor_specific_kmers_files.cross(splice_graph_files).flatten().collate(13)
    file reference_fa
    file gff
    
    output:
    file("${sample}_outcome_peptide_${chromosome}_8.txt")
    file("${sample}_outcome_peptide_${chromosome}_9.txt")
    file("${sample}_outcome_peptide_${chromosome}_10.txt")
    file("${sample}_outcome_peptide_${chromosome}_11.txt")
    file("${sample}_peptide_${chromosome}_8.xls")
    file("${sample}_peptide_${chromosome}_9.xls")
    file("${sample}_peptide_${chromosome}_10.xls")	
    file("${sample}_peptide_${chromosome}_11.xls")
    
    """
    python /NeoSplice/kmer_graph_inference.py --sample ${sample} --chromosome ${chromosome} --bam_file ${tumor} --gff_file ${gff} --genome_fasta ${reference_fa} --kmer_bam ${tumor_specific_kmers} --splice_graph ${splice_graph} --tumor_junction_file ${tumor_junctions} --normal_junction_file ${normal_junctions} --transcript_min_coverage 15 --HLA_I ${HLA} --netMHCpan_path /netMHCpan-4.0-docker/netMHCpan outdir .

    cp neoantigen_result/${sample}/${sample}_outcome_peptide_${chromosome}_8.txt  ${sample}_outcome_peptide_${chromosome}_8.txt
    cp neoantigen_result/${sample}/${sample}_outcome_peptide_${chromosome}_9.txt ${sample}_outcome_peptide_${chromosome}_9.txt
    cp neoantigen_result/${sample}/${sample}_outcome_peptide_${chromosome}_10.txt ${sample}_outcome_peptide_${chromosome}_10.txt
    cp neoantigen_result/${sample}/${sample}_outcome_peptide_${chromosome}_11.txt ${sample}_outcome_peptide_${chromosome}_11.txt
    cp neoantigen_result/${sample}/${sample}_peptide_${chromosome}_8.xls ${sample}_peptide_${chromosome}_8.xls
    cp neoantigen_result/${sample}/${sample}_peptide_${chromosome}_9.xls ${sample}_peptide_${chromosome}_9.xls
    cp neoantigen_result/${sample}/${sample}_peptide_${chromosome}_10.xls ${sample}_peptide_${chromosome}_10.xls
    cp neoantigen_result/${sample}/${sample}_peptide_${chromosome}_11.xls  ${sample}_peptide_${chromosome}_11.xls
    """
}
