# Need a Pathway-tools installer in the same folder.
FROM ubuntu:18.04

# To fix the issue with locale language.
ENV LANG C.UTF-8
ENV PYTHONIOENCODING=utf-8
# No interactive debian (like tzdata)
ENV DEBIAN_FRONTEND=noninteractive

LABEL maintainer="Meziane AITE & Arnaud BELCOUR"
LABEL Version="0.5"
LABEL Description="Metabolic Network comparison dockerfile."

# Install dependencies for Pathway-Tools and Orthofinder.
RUN apt-get -y update && \
    apt-get install -y \
    curl \
    make \
    wget \
    csh \
    git \
    ncbi-blast+ \
    libxm4 \ 
    vim \
    r-base \
    python3-pip \	
    python3.7 \
    python3.7-dev \
    iputils-ping \
    screen \
    exonerate \
    gnome-terminal;\
    echo "[ncbi]\nData=/usr/bin/data" > ~/.ncbirc

#Install Python packages
RUN python3.7 -m pip install requests \
	mpwt \
	padmet \
	aucome


# Install OrthoFinder.
RUN mkdir /programs/ /programs/diamond/ /shared/;\
    cd /programs;\
    wget https://github.com/davidemms/OrthoFinder/releases/download/2.3.11/OrthoFinder.tar.gz;\
    tar xzf OrthoFinder.tar.gz;\
    rm OrthoFinder.tar.gz;\
    echo 'export PATH="$PATH:/programs/OrthoFinder:"\nexport PATH="$PATH:/programs/OrthoFinder/bin:"\nexport OPENBLAS_NUM_THREADS=1' >> ~/.bashrc

# Install r upset for intervene and pvclust for dendrogram, padmet, mpwt and comparison script.
RUN R -e "install.packages('https://cran.r-project.org/src/contrib/Archive/UpSetR/UpSetR_1.3.3.tar.gz', repos=NULL, type='source', dependencies=TRUE)";\
    R -e "install.packages('https://cran.r-project.org/src/contrib/Archive/pvclust/pvclust_2.0-0.tar.gz', repos=NULL, type='source', dependencies=TRUE)";\
    R -e "install.packages('https://cran.r-project.org/src/contrib/ape_5.3.tar.gz', repos=NULL, type='source', dependencies=TRUE)"


#Change in future by python3.7 -m pip install aucome (need to update aucome pack)
#clone depots
RUN cd /programs;\
    git clone https://github.com/AuReMe/padmet-utils.git
