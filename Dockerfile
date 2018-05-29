# Need a Pathway-tools installer in the same folder.
# Use it with the command mpwt -f folder
FROM ubuntu:18.04
MAINTAINER "Meziane AITE & Arnaud BELCOUR"
LABEL Version="0.2"
LABEL Description="Metabolic Network comparison dockerfile."

# Install dependencies for Pathway-Tools and Orthofinder.
RUN apt-get -y update && \
    apt-get install -y \
    curl \
    wget \
    csh \
    git \
    ncbi-blast+ \
    mcl \
    libxm4 \
    gnome-terminal;\
    echo "[ncbi]\nData=/usr/bin/data" > ~/.ncbirc

# Install OrthoFinder.
RUN mkdir /shared;\
    mkdir /programs;\ 
    cd /programs;\
    wget https://github.com/davidemms/OrthoFinder/releases/download/v2.2.6/OrthoFinder-2.2.6.tar.gz;\
    tar xzf OrthoFinder-2.2.6.tar.gz;\
    rm OrthoFinder-2.2.6.tar.gz;\
    echo 'export PATH="$PATH:/programs/OrthoFinder-2.2.6:"' >> ~/.bashrc;\
    wget http://www.atgc-montpellier.fr/download/sources/fastme/fastme-2.1.5.tar.gz;\
    tar xzf fastme-2.1.5.tar.gz fastme-2.1.5/binaries/fastme-2.1.5-linux64;\
    mv fastme-2.1.5/binaries/fastme-2.1.5-linux64 /usr/local/bin/fastme;\
    rm -rf fastme-2.1.5*

# Install padmet, mpwt and comparison script.
RUN curl https://bootstrap.pypa.io/get-pip.py | python2.7;\
    cd /programs;\
    git clone https://gitlab.inria.fr/maite/padmet-utils.git;\
    pip2 install python-libsbml configparser;\ 
    pip2 install padmet mpwt;\
    cd /usr/bin;\
    wget https://gitlab.inria.fr/DYLISS/compare_metabo/raw/master/compare.py;\
    mv compare.py compare;\
    chmod u+x compare