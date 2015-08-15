"""
A pipeline that takes paired end reads from mRNA-Seq,
maps them to the reference transcriptome, calculates abundance,
and calculates differential expression.
"""

#make it so this can all be done in an Amazon AWS instance
#add commands that output read number and other various stats across processes
#into an output file for review later.
#I need to check versions so that I can get the exact results out that I did before

#update software
sudo bash
apt-get update
apt-get -y upgrade

#get other necessary software
apt-get -y install python-dev unzip python-pip fastqc wget screen default-jre samtools

#mount the hard drive
df -h
mkfs -t ext4 /dev/xvda1
mount /dev/xvda1 /mnt
chown -R ubuntu:ubuntu /mnt
df -h

#install packages
cd
#casava illumina filter 0.1	05-Aug-2011
wget http://cancan.cshl.edu/labmembers/gordon/fastq_illumina_filter/release/0.1/fastq_illumina_filter-0.1.tar.gz
tar -xzf fastq_illumina_filter-0.1.tar.gz
cd fastq_illumina_filter-0.1
make
cp fastq_illumina_filter /usr/local/bin
cd
#trimmomatic
wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.30.zip
unzip Trimmomatic-0.30.zip
cd Trimmomatic-0.30
chmod +x trimmomatic-0.30.jar
cd
#bowtie2 version2.1.0
wget http://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.1.0/bowtie2-2.1.0-linux-x86_64.zip
unzip bowtie2-2.1.0-linux-x86_64.zip
#eXpress get version eXpress 1.5.0
wget http://bio.math.berkeley.edu/eXpress/downloads/express-1.5.1/express-1.5.1-linux_x86_64.tgz
tar -zxvf express-1.5.1-linux_x86_64.tgz
mv express-1.5.1-linux_x86_64 /usr/bin/express-1.5.1-linux_x86_64

#make directories for data
cd /mnt
mkdir references
cd references

#download reference
wget ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r6.02_FB2014_05/fasta/dmel-all-transcript-r6.02.fasta.gz

#get samples reads
#for now this is just one sample, final version will
#will have new links and multiple samples in the pipeline
#will be like: for i,j reads in sample from file that contains sample list, do this:
scp and password
cat *R1_00* > 1GR1_read1.fastq.gz
cat *R2_00* > 1GR1_read2.fastq.gz
#get rid of seperate files
rm LPSH*

#DOING WORK ON THE READS NOW

#Remove reads not passing Casava filter
gunzip 1GR1_read1.fastq.gz
fastq_illumina_filter --keep N -v -v -o 1GR1_read1_filtered.fastq 1GR1_read1.fastq
gunzip 1GR1_read2.fastq.gz
fastq_illumina_filter --keep N -v -v -o 1GR1_read2_filtered.fastq 1GR1_read2.fastq
cd ..

#Trim reads with Trimmomatic
mkdir trimming
cd trimming
#will have custom adapter sequences in github pull 
java -jar $HOME/Trimmomatic-0.33/trimmomatic-0.30.jar PE \
-phred33 ../references/1GR1_read1_filtered.fastq ../references/1GR1_read2_filtered.fastq \
1GR1_read1_trimmed_paired.fastq.gz 1GR1_read1_trimmed_unpaired.fastq.gz \
1GR1_read2_trimmed_paired.fastq.gz 1GR1_read2_trimmed_unpaired.fastq.gz \
ILLUMINACLIP:TruSeq3custom-PE.fa:2:30:10:8:TRUE TRAILING:20 MINLEN:35

#create FastQC reports
fastqc ../trimming/1GR1_read1_trimmed_paired.fastq.gz
fastqc ../trimming/1GR1_read2_trimmed_paired.fastq.gz
cd ..
mkdir fastqc
cd fastqc
mv ../trimming/*fastqc* ./.
cd ..

#khmer and jellyfish QC here

#build Bowtie2 reference
mkdir bt2
cd bt2
gunzip ../references/dmel-all-transcript-r6.02.fasta.gz 
bowtie2-build --offrate 1 ../references/dmel-all-transcript-r6.02.fasta dmel-all-transcript-r6.02

#map reads with Bowtie2
#-q for reads in fastq format
#-a to search for and report all alignments
#-t to show the time
#--un-gz to output unmapped reads for later processing
#-p 8 to use 8 threads
bowtie2 -p 16 -t -q -a --un-gz unmapped.fastq.gz -x dmel-all-transcript-r6.02 -1 ../trimming/1GR1_read1_trimmed_paired.fastq.gz -2 ../trimming/1GR1_read2_trimmed_paired.fastq.gz -U ../trimming/1GR1_read1_trimmed_unpaired.fastq.gz,../trimming/1GR1_read2_trimmed_unpaired.fastq.gz | samtools view -Sb - > bt2out.bam
cd ..

#estimate abundance with eXpress
mkdir express
cd express
/usr/bin/express-1.5.1-linux_x86_64/express /mnt/references/dmel-all-transcript-r6.02.fasta ../bt2/bt2out.bam

#combine expression counts from all samples and get ready for input in DESeq2
#sum transcripts to gene level using TPM * 100(Harold's method)
#build my spreadsheet
#first I'll just use my old spreadsheet for easiness.
DE analysis
GO analysis