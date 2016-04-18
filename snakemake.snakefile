#!/usr/bin/python

"""
A pipeline that takes paired end fastq read files from mRNA-Seq,
maps them to the reference transcriptome, calculates abundance,
and calculates differential expression.
"""

include: "config.py"

FINAL = [] #files to be output

#this is default make, the final output files
rule all:
    input:
        #list all the data stuff here
        #transcriptome
        #gene to transcript
        #reads


#remove files to keep stuff clean
rule clean:
    shell: 'rm -f outfiles' #stuff to delete goes here


#install necessary software and set up environment
rule setup:
    input: "setup.sh"
    shell:
        'sudo bash {input[0]}'

#get references and annotation
rule get_genome:
     output:
        "references/dmel-all-transcript-r6.02.fasta"
     shell:
        "cd references; "
        "wget ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r6.02_FB2014_05/fasta/dmel-all-transcript-r6.02.fasta.gz; "
        "gunzip dmel-all-transcript-r6.02.fasta.gz"

rule get_genetotrans:
    output:
        "references/"
        ""
        http://flybase.org/static_pages/downloads/FB2015_05/genes/fbgn_fbtr_fbpp_fb_2015_05.tsv.gz


rule get_reads:
    #need to submit to SRA

#continuing the snakemake file using the following as template for next step:
'''
rule bwt2_idx:
    input:
        "{0}/ref.grp".format(ERCC_RSEM_DIR)
    output:
        "{0}/ref.transcripts.1.bt2".format(ERCC_RSEM_DIR),
        "{0}/ref.transcripts.2.bt2".format(ERCC_RSEM_DIR),
        "{0}/ref.transcripts.3.bt2".format(ERCC_RSEM_DIR),
        "{0}/ref.transcripts.4.bt2".format(ERCC_RSEM_DIR),
        "{0}/ref.transcripts.rev.1.bt2".format(ERCC_RSEM_DIR),
        "{0}/ref.transcripts.rev.2.bt2".format(ERCC_RSEM_DIR)
    shell:
        'bowtie2-build '
        '--seed 42 '
        '--offrate 1 '
        '{ERCC_RSEM_DIR}/ref.transcripts.fa '
        '{ERCC_RSEM_DIR}/ref.transcripts'

rule get_reads
    input: "sample_metadata.txt"
    output: ".fa"
    #SRA download here
    expand
'''
