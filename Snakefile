configfile: "config.yaml"

from snakemake.utils import min_version
import os
#from os.path import join
#import glob

min_version("3.6.0")

#### Functions

#Here we are collecting the ids of the different sample. the pattern is : xxx_R1_001.fastq.gz. Files MUST be named this way. it will be usefull for the variant calling step, when all of these file will be merged to one file. Without this, snakemake won't be able to do it.
sample_ids, = glob_wildcards("data/Raw_reads/{sample}_R1_001.fastq.gz")

#### Main worwflow rule. Each rule able to run specific parts of the workflow with ease. just launch snakemake with a rule as a paramater in order to produce each of the output files associated. ex : snakemake calling_only

#Run the entire pipeline
rule all:
    input:
        "report.html"
    message:
        "The whole pipeline has been executed"

#Run only the steps necessary to perform a Variant Calling. will not generate the graphs
rule calling_only:
    input:
        "Variant_calling.vcf"
    message:
        "Variant calling done"

#Will only do the mapping step
rule mapping_only:
    input:
        bam=expand("rmduped_reads/{sample}.bam", sample=sample_ids) if config["remove_duplicates"]!=0 else expand("sorted_reads/{sample}.bam", sample=sample_ids),
        bai=expand("rmduped_reads/{sample}.bam.bai", sample=sample_ids) if config["remove_duplicates"]!=0 else expand("sorted_reads/{sample}.bam.bai", sample=sample_ids),
    message:
        "Mapping done"

#WARNING : you !MUST! use the --nt option or all the resulting mapping files will be deleted. They are indeed marked as temp files.
rule raw_mapping_only:
    input:
        bam=expand("mapped_bam/{sample}.bam", sample=sample_ids)
    message:
        "Reads have been mapped against the reference genome and the resulting alignements have been converted to BAM files"

#Will only produce the steps necessary to create the graphs
rule diagnostic:
    input:
        "Mapping_stats.csv",
        "outputDepthBySample.csv",
        "outputDepthByRegion.csv"
    message:
        "The graphs have been created"

#Will create all the different indexes that the analyses needs.
rule required_index:
    input:
        expand("{genome}.bwt", genome=config["genome"]),
        expand("{genome}.fai", genome=config["genome"]),
        expand("{genome_p}.dict", genome_p=config["genome_prefix"])
    message:
        "All the required files related to the ref genome have been created"

#### General messages
# These messages will be print depending on the sucess or the failure of the pipeline
onsuccess:
    print("\nWorkflow finished. If the entire workflow has been runned, a file named report.html has been generated and contains all the output informations.\n")

onerror:
    print("An error occurred. You should read the readme file and see if everything have been setup the right way.")

#### Workflow

#Trimming step, can be set to Off in the config file.
rule cutadapt:
    input:
       read="data/Raw_reads/{sample}_R1_001.fastq.gz",
       read2="data/Raw_reads/{sample}_R2_001.fastq.gz"
    output:
        R1=temp("trimmed_reads/{sample}_R1_001.fastq.gz"),
        R2=temp("trimmed_reads/{sample}_R2_001.fastq.gz") 
    threads:
        50
    message:
        "Cutadapt on {input}..."
    priority:
        20
    log:
        "logs/trimming/{sample}.log"
    shell:
        "cutadapt -q {config[Cutadapt][Quality_value]} -m {config[Cutadapt][min_length]} -a {config[Cutadapt][forward_adapter]} -A  {config[Cutadapt][reverse_adapter]} -o {output.R1} -p {output.R2} {config[Cutadapt][options]} {input.read} {input.read2} > {log}"

#The bwt is necessary for BWA to run, so we build it.
rule bwa_build_bwt:
    input:
        genome=expand("{genome}", genome=config["genome"])
    output:
        "{genome}.bwt"
    log:
        "logs/index/bwa.log"
    message:
        "Building bwt for {input}"
    shell:
        "bwa index {input} > {log}" 

#The genome also need to be indexed by samtools
rule samtools_faidx:
    input:
        genome=expand("{genome}", genome=config["genome"])
    output:
        "{genome}.fai"
    log:
        "logs/index/samtools.log"
    message:
        "Indexing {input}..."
    shell:
       "samtools faidx {input} > {log}" 

#Raw mapping, with conversion to Bam file with samtools. 
rule bwa_map:
    input:
       genome=expand("{genome}", genome=config["genome"]),
       bwt=expand("{genome}.bwt", genome=config["genome"]),
       #read=expand("trimmed_reads/{{sample}}_{pair}_001.fastq.gz", pair=["R1", "R2"])
       read=expand("trimmed_reads/{{sample}}_{pair}_001.fastq.gz", pair=["R1", "R2"]) if config["trimming"]!=0 else expand("data/Raw_reads/{{sample}}_{pair}_001.fastq.gz", pair=["R1", "R2"])
    output:
        temp("mapped_bam/{sample}.bam")
    log:
        "logs/mapping/{sample}.log"
    params:
        rg="@RG\\tID:{sample}\\tPL:ILLUMINA\\tSM:{sample}"
    message:
        "Raw mapping with {input.read}..."
    benchmark:
        "Benchmark/bwa/{sample}-bench.txt"
    shell:
        "bwa mem -t {config[BWA][t]} {config[BWA][options]} -R '{params.rg}' {input.genome} {input.read} 2> {log} | samtools view -Sb - > {output} "

#Sorting reads with Picard. Will be the last step before variant calling if rmdup is set to false.
rule picard_sort:
    input:
        "mapped_bam/{sample}.bam"
    output:
        "sorted_reads/{sample}.bam"
    log:
        "logs/sorting/{sample}.log"
    message:
        "Sorting {input}..."
    shell:
        "java -Xmx4g -jar {config[software_reps][picard]} "
        "SortSam I={input} "
        "O={output} "
        "SO=coordinate "
        "VALIDATION_STRINGENCY=SILENT "
        "2> {log}"

#Remove duplicates from sorted reads. Can be diseable in the config file.
rule picard_rmdup:
    input:
        bam="sorted_reads/{sample}.bam"
    output:
        "rmduped_reads/{sample}.bam",
    log:
        log1="logs/picard_rmdup/{sample}.log",
        log2="logs/picard_rmdup/{sample}.log2"
    message:
        "Removing duplicates in {input}..."
    shell:
        "java -jar -Xmx2g {config[software_reps][picard]} "
        "MarkDuplicates "
        "I={input} "
        "O={output} "
        "VALIDATION_STRINGENCY=SILENT "
        "MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 "
        "REMOVE_DUPLICATES=TRUE "
        "M={log.log1} "
        "2> {log.log2}"

#Index the bam file. Will index the sorted reads or the rmduped reads depending if rmdup is enable or not.
rule samtools_index:
    input:
        "rmduped_reads/{sample}.bam" if config["remove_duplicates"]!=0 else "sorted_reads/{sample}.bam"
    output:
        "rmduped_reads/{sample}.bam.bai" if config["remove_duplicates"]!=0 else "sorted_reads/{sample}.bam.bai"
    message:
        "Indexing {input}..."
    priority:
        10
    shell:
        "samtools index {input}"

#Create genome dictionary with Picard
rule picard_dict:
    input:
        genome=expand("{genome_p}.fasta", genome_p=config["genome_prefix"])
    output:
        "{genome_p}.dict"
    message:
        "Creating dictionary for {input}..."
    shell:
        "java -jar {config[software_reps][picard]} "
        "CreateSequenceDictionary "
        "R={input} "
        "O={output}"

#The variant calling is done in 3 steps.

#VC - Step 1 : Create one gvcf file for each sample.
rule GATK_raw_calling:
    input:
        bam="rmduped_reads/{sample}.bam" if config["remove_duplicates"]!=0 else "sorted_reads/{sample}.bam",
        bai="rmduped_reads/{sample}.bam.bai" if config["remove_duplicates"]!=0 else "sorted_reads/{sample}.bam.bai",
        genome=expand("{genome}", genome=config["genome"]),
        dictionary=expand("{genome_p}.dict", genome_p=config["genome_prefix"]),
        index_genome=expand("{genome}.fai", genome=config["genome"])
    output:
        "Raw_calling/{sample}.g.vcf",
    priority:
        10
    log:
        normal_log="logs/Raw_calling/logs/{sample}.log",
        error_log="logs/Raw_calling/error/{sample}.err"
    message:
        "Proceding raw GVCF calling on {input.bam}"
    benchmark:
        "Benchmark/Raw_calling/{sample}-bench.txt"
    shell:
        "java -Xmx4g -jar {config[software_reps][GATK]} "
        "-ploidy {config[GATK][ploidy]} "
        "--emitRefConfidence GVCF "
        "-T HaplotypeCaller "
        "-R {input.genome} "
        "-I {input.bam} "
        "-L {config[bed_file]} "
        "--genotyping_mode DISCOVERY "
        "{config[GATK][bed_file]} "
        "{config[GATK][option]}"
        "-o {output} "
        "> {log.error_log} "
        "2> {log.normal_log}"

#VC - Step 1.bis : The following steps needs a list of the gvcf files to merge. This list is created here. The dictionary created at the begining 
rule GATK_file_list:
    input:
        expand("Raw_calling/{sample}.g.vcf", sample=sample_ids)
    output:
        temp("fileList.list")
    priority:
        10
    message:
        "Creating the list of input files for GATK"
    shell:
        "ls Raw_calling/*.vcf > {output}"

#VC - Step 2 : Merge all the GVCF files.
rule GATK_merge:
    input:
       listing="fileList.list",
       genome=expand("{genome}", genome=config["genome"]),
       dictionary=expand("{genome_p}.dict", genome_p=config["genome_prefix"]),
       index_genome=expand("{genome}.fai", genome=config["genome"])
    output:
        vcf=temp("merged_raw_calling.g.vcf"),
        vcf_idx=temp("merged_raw_calling.g.vcf.idx")
    priority:
        10
    log:
        normal_log="logs/Merging_GVCF/Merging_GVCF.log",
        error_log="logs/Merging_GVCF/Merging_GVCF.err"
    message:
        "Merging GVCFs files..."
    shell:
        "java -Xmx4g -jar {config[software_reps][GATK]} "
        "-T CombineGVCFs "
        "-R {input.genome} "
        "--variant {input.listing} "
        "-o {output.vcf} "
        "> {log.error_log} "
        "2> {log.normal_log}"

#VC - Step 3 : Convert the resulting gvcf file to a vcf file. End of the variant calling.
rule GATK_final_calling:
    input:
        merged_gvcf="merged_raw_calling.g.vcf",
        merged_gvcf_idx="merged_raw_calling.g.vcf.idx",
        genome=expand("{genome}", genome=config["genome"]),
        dictionary=expand("{genome_p}.dict", genome_p=config["genome_prefix"]),
        index_genome=expand("{genome}.fai", genome=config["genome"])
    output:
        "Variant_calling.vcf"
    priority:
        10
    log:
        normal_log="logs/Final_calling/Variant_calling.log",
        error_log="logs/Final_calling/error/Variant_calling.err"
    message:
        "Proceding variant calling..."
    shell:
        "java -Xmx4g -jar {config[software_reps][GATK]} "
        "-T GenotypeGVCFs "
        "-R {input.genome} "
        "--sample_ploidy 2 "
        "--variant {input.merged_gvcf} "
        "-o {output} "
        "> {log.error_log}"
        "2> {log.normal_log}"

#Create the list of bam files needed for the graphs.
rule BAM_file_list:
    input:
        expand("rmduped_reads/{sample}.bam", sample=sample_ids) if config["remove_duplicates"]!=0 else expand("sorted_reads/{sample}.bam", sample=sample_ids),
    output:
        temp("bamList.list")
    message:
        "Creating the bam files list"
    shell:
        "ls rmduped_reads/*.bam > bamList.list" if config["remove_duplicates"]!=0 else "ls sorted_reads/*.bam > bamList.list"

#Compute the depth of sequencing for each position in the targeted Bed. Is used as an input for the graphs.
rule samtools_depth:
    input:
        "bamList.list"
    output:
        temp("samtools_depth.csv")
    priority:
        0
    message:
        "analyzing the depth of the bam files..."
    shell:
        "samtools depth -aa -b {config[bed_file]} -f {input} > {output}"

rule depth_by_sample:
    input:
        listing="bamList.list",
        depth="samtools_depth.csv"
    output:
        "depth_by_sample_n4_s1_i5.pdf",
        "depth_by_sample_n4_s1_i5.jpg",
        "outputDepthBySample.csv"
    log:
        "logs/graphes/depth_by_sample.log"
    priority:
        0
    message:
        "Creating depth by sample graph..."
    shell:
        "perl {config[depth_by_sample][path]} -l {input.listing} -b {config[bed_file]} -g -n {config[depth_by_sample][interval_size]} -s {config[depth_by_sample][sorting_col]} -i {config[depth_by_sample][number_of_interval]} > {log}"

rule depth_by_region:
    input:
        listing="bamList.list",
        depth="samtools_depth.csv"
    output:
        "depth_by_region_n4_s1_i5.pdf",
        "depth_by_region_n4_s1_i5.jpg",
        "outputDepthByRegion.csv"
    log:
        "logs/graphes/depth_by_region.log"
    priority:
        0
    message:
        "Creating depth by region graph"
    shell:
        "perl {config[depth_by_region][path]} -l {input.listing} -b {config[bed_file]} -g -n {config[depth_by_region][interval_size]} -s {config[depth_by_region][sorting_col]} -i {config[depth_by_region][number_of_interval]} > {log} "

rule mapping_trimming_stats:
    input:
        listing="bamList.list",
        bams=expand("mapped_bam/{sample}.bam", sample=sample_ids)
    output:
        "Mapping_stats.csv",
        "graphe_mapping.pdf",
        "graphe_mapping.jpg"
    log:
        "logs/graphes/mapping_trimming.log"
    priority:
        0
    message:
        "Creating mapping/trimming graph..."
    shell:
        "perl {config[mapping_trimming_stats][path]} -r data/Raw_reads -rb mapped_bam -pb rmduped_reads -b {config[bed_file]} -g" if config["remove_duplicates"]!=0 else "perl {config[mapping_trimming_stats][path]} -r data/Raw_reads -rb mapped_bam -b {config[bed_file]} -g > {log} 2>> {log}"
        
#Create the report...
rule report:
    input:
        VC="Variant_calling.vcf",
        DS="depth_by_sample_n4_s1_i5.jpg",
        DR="depth_by_region_n4_s1_i5.jpg",
        GM="graphe_mapping.jpg",
        file_mapping="Mapping_stats.csv",
        file_DBsample="outputDepthBySample.csv",
        file_DBregion="outputDepthByRegion.csv"
    output:
        "report.html"
    message:
        "Creating the report..."
    run:
        from snakemake.utils import report
        report("""
            Report of the Capture data analyses workflow 
            ============================================

            .. section-numbering::

            This report has been generated at the end of the capture data analyses workflow. You will find here some data related to some steps of the analyses. All of the included graphs are followed by a link that able the user to see the datas that have been used in order to create it.

            .. contents:: Table of Contents

            Parameters of the analyses
            --------------------------

            Cutadapt : trimming
            ~~~~~~~~~~~~~~~~~~~

            - Quality cutoff : {config[Cutadapt][Quality_value]}
            - Minimum read length : {config[Cutadapt][min_length]}
            - Forward adapter : {config[Cutadapt][forward_adapter]}
            - Reverse adapter : {config[Cutadapt][reverse_adapter]}
            - Additional option : {config[Cutadapt][options]}

            Mapping
            ~~~~~~~

            - Additional options : {config[BWA][options]}

            Variant calling : GATK
            ~~~~~~~~~~~~~~~~~~~~~~
    
            - Ploidy : {config[GATK][ploidy]}
            - Additional option : {config[GATK][option]}

            General informations about mapping
            ----------------------------------

            .. |date| date::
            .. |time| date:: %H:%M
            
            .. footer:: This document was generated on |date| at |time|. Author : Alexandre Soriano (INRA) : alex.soriano@laposte.net

            Reads were mapped to *{config[genome]}* using BWA MEM in paired-end mode. After that they have been sorted using picard-tools.

            The remove duplicate option was set to **{config[remove_duplicates]}**

            The trimming option was set to **{config[trimming]}**

            .. figure:: graphe_mapping.jpg
                :alt: Descriptions of mapping results

                Number of reads depending of the different steps of the mapping.

            Detailed data for this figure can be found  `here (file) <{input.file_mapping}>`_ or in the `Detailed mapping stats`_ paragraph.

            Informations about depth of sequencing
            --------------------------------------

            Depth analyses depending of the sample
            ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            .. figure:: depth_by_sample_n4_s1_i5.jpg
                :alt: depth depending of different samples

                Number of regions having a depth included in different intervals for each sample.

            Detailed data for this figure can be found `In this place (file) <{input.file_DBsample}>`_ or in the `Detailed Depth by sample stats`_ part.

            Depth analyses depending of targeted regions
            ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            .. figure:: depth_by_region_n4_s1_i5.jpg
                :alt: depth depending of targeted regions

                Number of targets having a depth included in different intervals for each region.

            Detailed data for this figure can be found `there (file) <{input.file_DBregion}>`_ or in the `Detailed Depth by region stats`_ part.

            Informations about Variant calling
            ----------------------------------

            The file resulting of the Variant calling step can be found `here <{input.VC}>`_ .

            Detailed informations of the different graphs
            ---------------------------------------------

            Detailed mapping stats
            ~~~~~~~~~~~~~~~~~~~~~~

            .. csv-table::
                :file: {input.file_mapping}
                :delim: U+0009

            Detailed Depth by sample stats
            ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            .. csv-table::
                :file: {input.file_DBsample}
                :delim: U+0009

            Detailed Depth by region stats
            ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            .. csv-table::
                :file: {input.file_DBregion}
                :delim: U+0009

            """, output[0])

#Created by Alexandre Soriano - 2017 - alex.soriano@laposte.net
