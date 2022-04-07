# Need a Pathway-tools installer in the same folder.
FROM ubuntu:18.04

# To fix the issue with locale language.
ENV LANG C.UTF-8
ENV PYTHONIOENCODING=utf-8
# No interactive debian (like tzdata)
ENV DEBIAN_FRONTEND=noninteractive

LABEL maintainer="Meziane AITE, Arnaud BELCOUR & Jeanne GOT"
LABEL Version="0.5.1"
LABEL Description="Metabolic Network comparison dockerfile."

RUN mkdir -p /programs/pathway-tools /programs/diamond/
RUN mkdir -p /home/database/BIOCYC/METACYC/25/

# Install dependencies for Pathway-Tools and Orthofinder.
RUN apt-get -y update && \
    DEBIAN_FRONTEND=noninteractive apt-get install -yq \
    curl \
    make \
    cmake \
    wget \
    csh \
    git \
    ncbi-blast+ \
    libffi-dev \
    libxm4 \ 
    vim \
    python3.7-dev \
    python3-pip \
    iputils-ping \
    exonerate \
    gnome-terminal \
    r-base \
    mcl
RUN apt-get clean \
    apt-get purge

## Install OrthoFinder and its dependencies.
# Install OrthoFinder
RUN wget https://github.com/davidemms/OrthoFinder/releases/download/2.5.4/OrthoFinder.tar.gz;\
    tar xzf OrthoFinder.tar.gz;\
    mv OrthoFinder/ /programs/;\
    rm OrthoFinder.tar.gz
# Install FastME
RUN wget http://www.atgc-montpellier.fr/download/sources/fastme/fastme-2.1.5.tar.gz;\
    tar xzf fastme-2.1.5.tar.gz fastme-2.1.5/binaries/fastme-2.1.5-linux64;\
    mv fastme-2.1.5/binaries/fastme-2.1.5-linux64 /usr/local/bin/fastme;\
    rm -r fastme-2.1.5*
# Install MMseqs2
RUN wget https://mmseqs.com/latest/mmseqs-linux-avx2.tar.gz;\
    tar xzf mmseqs-linux-avx2.tar.gz;\
    cp mmseqs/bin/mmseqs /usr/local/bin/;\
    mv mmseqs/ /programs/;\
    rm mmseqs-linux-avx2.tar.gz
# Install Diamond
RUN wget https://github.com/bbuchfink/diamond/releases/latest/download/diamond-linux64.tar.gz;\
    tar xzf diamond-linux64.tar.gz;\
    mv diamond /programs/diamond/;\
    rm diamond-linux64.tar.gz
 
RUN echo "[ncbi]" > ~/.ncbirc;\
    echo "Data=/usr/bin/data" >> ~/.ncbirc;\
    echo 'export PATH="$PATH:/programs/OrthoFinder:/programs/OrthoFinder/bin:/programs/mmseqs/bin/:/programs/diamond:"' >> ~/.bashrc
    

# Install r upset for pvclust for dendrogram, padmet, mpwt and comparison
# script.
RUN R -e "install.packages('https://cran.r-project.org/src/contrib/Archive/pvclust/pvclust_2.0-0.tar.gz', repos=NULL, type='source', dependencies=TRUE)";\
    R -e "install.packages('ape', dependencies=TRUE)"

#Install Python package
RUN python3 -m pip install --upgrade pip;\
    python3 -m pip install --upgrade setuptools

# Install Padmet
RUN cd /programs/;\
    git clone https://github.com/AuReMe/padmet.git;\
    cd /programs/padmet;\
    pip3 install -r requirements.txt;\
    chmod +x setup.py;\
    python3 setup.py install

# Install aucome
RUN cd /programs/;\
    git clone https://github.com/AuReMe/aucome.git;\
    cd /programs/aucome;\
    pip3 install -r requirements.txt;\
    chmod +x setup.py;\
    python3 setup.py install
