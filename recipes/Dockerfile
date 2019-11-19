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
    mcl \
    libxm4 \ 
    vim \
    r-base \
    python3-pip \	
    python3.7 \
    python3.7-dev \
    iputils-ping \
    screen \
    gnome-terminal;\
    echo "[ncbi]\nData=/usr/bin/data" > ~/.ncbirc

RUN cd tmp;\
    curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py;\
    python3.7 get-pip.py

#Install Python packages
RUN python3.7 -m pip install requests \
	mpwt \
	eventlet \
	padmet \
	aucome


# Install OrthoFinder.
RUN mkdir /programs/ /programs/diamond/ /shared/;\
    cd /programs;\
    wget https://github.com/davidemms/OrthoFinder/releases/download/2.3.3/OrthoFinder-2.3.3.tar.gz;\
    tar xzf OrthoFinder-2.3.3.tar.gz;\
    rm OrthoFinder-2.3.3.tar.gz;\
    wget http://www.atgc-montpellier.fr/download/sources/fastme/fastme-2.1.5.tar.gz;\
    tar xzf fastme-2.1.5.tar.gz fastme-2.1.5/binaries/fastme-2.1.5-linux64;\
    mv fastme-2.1.5/binaries/fastme-2.1.5-linux64 /usr/local/bin/fastme;\
    rm -rf fastme-2.1.5*;\
    cd diamond;\
    wget https://github.com/bbuchfink/diamond/releases/download/v0.9.22/diamond-linux64.tar.gz;\
    tar xzf diamond-linux64.tar.gz;\
    echo 'export PATH="$PATH:/programs/OrthoFinder-2.3.3:"\nexport PATH="$PATH:/programs/diamond:"' >> ~/.bashrc

#Change in future by curl only version, for now exonerate servers are down.
#Install exonerate
COPY exonerate-2.2.0 /bin/
RUN echo 'export PATH="$PATH:/bin/exonerate-2.2.0"' >> ~/.bashrc
 
# Install r upset for intervene and pvclust for dendrogram, padmet, mpwt and comparison script.
RUN R -e "install.packages('UpSetR',,dependencies=TRUE, repos='http://cran.rstudio.com/' )";\
    R -e "install.packages('pvclust',,dependencies=TRUE, repos='http://cran.rstudio.com/' )"

#Copy metacyc database
RUN mkdir -p /home/data/database
COPY database/* /home/data/database/

#Change in future by python3.7 -m pip install aucome (need to update aucome pack)
#clone depots
RUN cd /programs;\
    git clone https://github.com/AuReMe/padmet-utils.git

 

	

