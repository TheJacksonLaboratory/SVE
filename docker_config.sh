#!/bin/bash
apt-get update && apt-get install -y \
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

cd /software
git clone https://github.com/mysql/mysql-connector-python.git
cd mysql-connector-python
python ./setup.py build
sudo python ./setup.py install
cd /

pip install -I numpy
pip install -I scipy
pip install -I paramiko
pip install -I subprocess32
pip install -I HTSeq

