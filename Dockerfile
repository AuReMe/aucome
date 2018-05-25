# Need a Pathway-tools installer in the same folder.
# Use it with the command mpwt -f folder
FROM phusion/baseimage:latest
MAINTAINER "Meziane AITE & Arnaud BELCOUR"
LABEL Version="0.1"
LABEL Description="Comparison dockerfile."

# Install dependencies for Pathway-Tools and Orthofinder.
RUN apt-get -y update && \
    apt-get install -y \
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
    wget http://www.atgc-montpellier.fr/download/sources/fastme/fastme-2.1.5.tar.gz;\
    tar xzf fastme-2.1.5.tar.gz fastme-2.1.5/binaries/fastme-2.1.5-linux64;\
    mv fastme-2.1.5/binaries/fastme-2.1.5-linux64 /usr/local/bin/fastme;\
    rm -rf fastme-2.1.5*

# Install padmet, mpwt and comparison script.
# TODO add:
# cd /usr/bin;\
# wget https://gitlab.inria.fr/DYLISS/compare_metabo/blob/master/compare;\
# chmod u+x compare;\
# But the gitlab repository needs to be public for this.
RUN curl https://bootstrap.pypa.io/get-pip.py | python2.7;\
    cd /programs;\
    git clone https://gitlab.inria.fr/maite/padmet-utils.git;\
    pip2 install python-libsbml configparser;\ 
    pip2 install padmet mpwt