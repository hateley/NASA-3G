#!/bin/bash

"""
Installs the software necessary for running the data analysis.
"""


#get necessary supporting software
apt-get -y install \
unzip fastqc wget screen default-jre samtools


#install packages
cd
#casava illumina filter 0.1	05-Aug-2011
wget http://cancan.cshl.edu/labmembers/gordon/fastq_illumina_filter/release/0.1/fastq_illumina_filter-0.1.tar.gz
tar -xzf fastq_illumina_filter-0.1.tar.gz
rm fastq_illumina_filter-01.tar.gz
cd fastq_illumina_filter-0.1
make
cp fastq_illumina_filter /usr/local/bin
cd
#trimmomatic
wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.30.zip
unzip Trimmomatic-0.30.zip
rm Trimmomatic-0.30.zip
cd Trimmomatic-0.30
chmod +x trimmomatic-0.30.jar
cd
#bowtie2 version2.1.0
wget http://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.1.0/bowtie2-2.1.0-linux-x86_64.zip
unzip bowtie2-2.1.0-linux-x86_64.zip
rm bowtie2-2.1.0-linux-x86_64.zip
#eXpress get version eXpress 1.5.1
wget http://bio.math.berkeley.edu/eXpress/downloads/express-1.5.1/express-1.5.1-linux_x86_64.tgz
tar -zxvf express-1.5.1-linux_x86_64.tgz
mv express-1.5.1-linux_x86_64 /usr/bin/express-1.5.1-linux_x86_64
rm express-1.5.1-linux_x86_64.tgz

#prepare for next step (probably need to make this less jankedy)
chown -R ubuntu:ubuntu /mnt
cd /mnt
mkdir references
mkdir reads
