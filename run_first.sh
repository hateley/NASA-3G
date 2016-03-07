#!/bin/bash

#Do this first thing.
#Prepares a virtual machine for running the snakefile.
#Tested on 'Amazon EC2 Ubuntu E3.2xl' instance.
#You need about 100G hard drive.

#Run with: sudo bash run_first.sh

echo "Getting your computer ready to run the snakefile."
sleep 5

#update software
apt-get -y update
apt-get -y upgrade

#get snakemake and git
apt-get -y install python3-pip git
pip3 install snakemake

#clone github repo with snakefile, etc.
#git clone ...

echo
echo
echo "Snakemake setup complete."
echo
echo "Type 'snakemake -h' to ensure the install was successful."
echo "Then continue with snakefile."
