"""
A pipeline that takes paired end reads from mRNA-Seq,
maps them to the reference transcriptome, calculates abundance,
and calculates differential expression.
"""

#make it so this can all be done in an Amazon AWS instance



#update software
sudo bash
apt-get update
apt-get -y upgrade

#get other necessary software
apt-get -y install python-dev unzip python-pip fastqc wget screen default-jre

#mount the hard drive
df -h
mkfs -t ext4 /dev/xvdb
mount /dev/xvdb /mnt
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
wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.33.zip
unzip Trimmomatic-0.33.zip
cd Trimmomatic-0.33
chmod +x trimmomatic-0.33.jar
cd
#bowtie2 version2.1.0
apt-get install bowtie2
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
cp ~/Trimmomatic-0.33/adapters/TruSeq3-PE.fa ./.
java -jar $HOME/Trimmomatic-0.33/trimmomatic-0.33.jar PE \
-phred33 ../references/1GR1_read1_filtered.fastq ../references/1GR1_read2_filtered.fastq \
1GR1_read1_trimmed_paired.fastq.gz 1GR1_read1_trimmed_unpaired.fastq.gz \
1GR1_read2_trimmed_paired.fastq.gz 1GR1_read2_trimmed_unpaired.fastq.gz \
ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:8:TRUE TRAILING:20 MINLEN:35

#create FastQC reports
fastqc ../trimming/1GR1_read1_trimmed_paired.fastq.gz
fastqc ../trimming/1GR1_read2_trimmed_paired.fastq.gz
cd ..
mkdir fastqc
cd fastqc
mv ../trimming/*fastqc* ./.
cd ..

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
bowtie2 -p 8 -t -q -a --un-gz unmapped.fastq.gz \
	-x ~/references/bt2build/dmel-all-transcript-r6.02 \
	-1 1GR1_read1_trimmed_paired.fastq.gz -2 1GR1_read2_trimmed_paired.fastq.gz \
	-U 1GR1_read1_trimmed_unpaired.fastq.gz,1GR1_read2_trimmed_unpaired.fastq.gz \
	| samtools view -bS - > bt2out.bam

#estimate abundance with eXpress
express ~/references/dmel-all-transcript-r6.02.fasta bt2out.bam

sum transcripts to gene level using TPM * 100(Harold's method)
DE analysis
GO analysis

do
curl -u single-cell:SanDiegoCA    http://bix.ucsd.edu/projects/singlecell/nbt_data/ecoli_mda_lane${i}.fastq.bz2 |bunzip2 - |head -400000 > ecoli_mda_lane${i}.fastq
bwa index sequence.fasta
bwa aln sequence.fasta ecoli_mda_lane${i}.fastq >ecoli_lane${i}.sai
bwa samse sequence.fasta ecoli_lane${i}.sai ecoli_mda_lane${i}.fastq > ecoli_lane${i}.sam
python ./sam-scan-errhist.py -o sam-scan-errhist_ecoli_lane${i}.out sequence.fasta ecoli_lane${i}.sam
done