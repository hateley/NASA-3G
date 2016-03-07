#!/usr/bin/python

"""
A pipeline that takes paired end fastq read files from mRNA-Seq,
maps them to the reference transcriptome, calculates abundance,
and calculates differential expression.
"""


FINAL = [] #files to be output

#this is default make, the final output files
rule all:
    input: FINAL

#remove files to keep stuff clean
rule clean:
    shell: 'rm -f outfiles' #stuff to delete goes here


#install necessary software
rule get_software:
    input: "setup.sh"
    output: "dmel-all-transcript-r6.02.fasta.gz"
    shell:
        'sudo bash {input[0]}'


rule get_genome:
     output:
        "references/dmel-all-transcript-r6.02.fasta"
     shell:
        "mkdir references; "
        "cd references; "
        "wget ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r6.02_FB2014_05/fasta/dmel-all-transcript-r6.02.fasta.gz; "
        "gunzip dmel-all-transcript-r6.02.fasta.gz"


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
