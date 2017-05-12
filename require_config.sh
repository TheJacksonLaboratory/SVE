#!/bin/bash

apt-get update \
apt-get install -y \
ant \
build-essential \
cmake \
g++ \
gfortran \
git \
wget \
nano \
libffi6 \
libffi-dev \
libssl1.0.0 \
libssl-dev \
libblas3 \
libblas-dev \
liblapack3 \
liblapack-dev \
libncurses-dev \
libbz2-dev \
zlib1g-dev \
python \
python-dev \
python-pip \
&& apt-get clean

#cd ../
#git clone https://github.com/mysql/mysql-connector-python.git
#cd mysql-connector-python
#python ./setup.py build
#sudo python ./setup.py install
#cd ../


pip install --upgrade pip
pip install --user -I argeparse
pip install --user -I numpy
pip install --user -I scipy
pip install --user -I paramiko
pip install --user -I subprocess32
pip install --user -I pysam
pip install --user -I pysamstats
pip install --user -I HTSeq
pip install --user -I crossmap

#echo "checking SVE executables"
#scripts/variant_processor.py -h

