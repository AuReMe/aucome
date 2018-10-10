# Need a Pathway-tools installer in the same folder.
FROM ubuntu:16.04

# To fix the issue with locale language.
ENV LANG C.UTF-8
ENV PYTHONIOENCODING=utf-8

LABEL maintainer="Meziane AITE & Arnaud BELCOUR"
LABEL Version="0.4"
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
    python3.5-dev \
    iputils-ping \
    gnome-terminal;\
    echo "[ncbi]\nData=/usr/bin/data" > ~/.ncbirc

# Install OrthoFinder.
RUN mkdir /programs/ /shared/;\
    cd /programs;\
    wget https://github.com/davidemms/OrthoFinder/releases/download/v2.2.7/OrthoFinder-2.2.7.tar.gz;\
    tar xzf OrthoFinder-2.2.7.tar.gz;\
    rm OrthoFinder-2.2.7.tar.gz;\
    echo 'export PATH="$PATH:/programs/OrthoFinder-2.2.7:"' >> ~/.bashrc;\
    wget http://www.atgc-montpellier.fr/download/sources/fastme/fastme-2.1.5.tar.gz;\
    tar xzf fastme-2.1.5.tar.gz fastme-2.1.5/binaries/fastme-2.1.5-linux64;\
    mv fastme-2.1.5/binaries/fastme-2.1.5-linux64 /usr/local/bin/fastme;\
    rm -rf fastme-2.1.5*;\
    wget https://mmseqs.com/latest/mmseqs-static_avx2.tar.gz;\
    tar xvzf mmseqs-static_avx2.tar.gz;\
    rm mmseqs-static_avx2.tar.gz;\
    echo 'export PATH="$PATH:/programs/mmseqs2/bin/:"' >> ~/.bashrc


# Install r upset for intervene, padmet, mpwt and comparison script.
RUN R -e "install.packages('UpSetR',,dependencies=TRUE, repos='http://cran.rstudio.com/' )";\
    curl https://bootstrap.pypa.io/get-pip.py | python3;\
    cd /programs;\
    git clone https://gitlab.inria.fr/maite/padmet-utils.git;\
    pip3 install configparser padmet mpwt eventlet requests seaborn sklearn fastcluster intervene;\
    cd /usr/bin;\
    wget https://gitlab.inria.fr/DYLISS/compare_metabo/raw/master/compare.py;\
    mv compare.py compare;\
    chmod u+x compare

COPY database /home/database/
