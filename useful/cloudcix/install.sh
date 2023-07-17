#!/bin/bash

sudo apt-get update && sudo apt-get upgrade

# install basics, for some reason wouldn't work as a multiline command
sudo apt-get -y install vim
sudo apt-get -y install git
sudo apt-get -y install default-jre
sudo apt-get -y install unzip
sudo apt-get -y install tar
sudo apt-get -y install nano
sudo apt-get -y install pip 
sudo apt-get -y install build-essentials
sudo apt-get -y install gcc-multilib
sudo apt-get -y install libncurses-dev
sudo apt-get -y install tmux
sudo apt-get -y install htop
sudo apt-get -y install curl
sudo apt-get -y install ncdu
sudo apt-get -y install sendmail
# install anaconda dependencies
sudo apt-get install libgl1-mesa-glx \
    libegl1-mesa \
    libxrandr2 \
    libxss1 \
    libxcursor1 \
    libxcomposite1 \
    libasound2 \
    libxi6 \
    libxtst6

wget https://repo.anaconda.com/archive/Anaconda3-2021.05-Linux-x86_64.sh
bash Anaconda3-2021.05-Linux-x86_64.sh

# conda create --yes --name pyega python=3.8
# install pyega
#conda activate pyega
#conda config --add channels bioconda
#conda config --add channels conda-forge
#conda install --yes pyega3

# install singularity
## dependencies
sudo apt-get update && sudo apt-get install -y \
    build-essential \
    uuid-dev \
    libgpgme-dev \
    squashfs-tools \
    libseccomp-dev \
    wget \
    pkg-config \
    git \
    cryptsetup-bin

# install GO

export VERSION=1.13.5 OS=linux ARCH=amd64 && \
    wget https://dl.google.com/go/go$VERSION.$OS-$ARCH.tar.gz && \
    sudo tar -C /usr/local -xzvf go$VERSION.$OS-$ARCH.tar.gz && \
    rm go$VERSION.$OS-$ARCH.tar.gz
# add go to path
echo 'export GOPATH=${HOME}/go' >> ~/.bashrc && \
    echo 'export PATH=/usr/local/go/bin:${PATH}:${GOPATH}/bin' >> ~/.bashrc && \
    source ~/.bashrc

# download singularity stable release
export VERSION=3.8.4 && \
    wget https://github.com/apptainer/singularity/releases/download/v${VERSION}/singularity-${VERSION}.tar.gz && \
    tar -xzf singularity-${VERSION}.tar.gz && \
    cd singularity-${VERSION}
# configure singularity
./mconfig && \
    make -C ./builddir && \
    sudo make -C ./builddir install

# install latest nextflow release
wget -qO- https://get.nextflow.io | bash
chmod 777 nextflow
sudo mv ./nextflow /usr/local/bin/
nextflow -v
