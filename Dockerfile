# Need a Pathway-tools installer in the same folder.
# Use it with the command mpwt -f folder
FROM ubuntu:16.04
LABEL maintainer="Meziane AITE & Arnaud BELCOUR"
LABEL Version="0.2"
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
    gnome-terminal;\
    echo "[ncbi]\nData=/usr/bin/data" > ~/.ncbirc

# Install OrthoFinder.
RUN mkdir /programs/ /shared/ -p /data/database/BIOCYC/Metacyc/22.0_enhanced/ /data/database/MNX/;\
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
    echo 'export PATH="$PATH:/programs/mmseqs2/bin/:"' >> ~/.bashrc


# Install padmet, mpwt and comparison script.
RUN curl https://bootstrap.pypa.io/get-pip.py | python3;\
    cd /programs;\
    wget https://gitlab.inria.fr/DYLISS/compare_metabo/raw/master/ptools_installer;\
    git clone https://gitlab.inria.fr/maite/padmet-utils.git;\
    pip2 install python-libsbml configparser;\ 
    pip2 install padmet mpwt eventlet;\
    cd /usr/bin;\
    wget https://gitlab.inria.fr/DYLISS/compare_metabo/raw/master/compare.py;\
    mv compare.py compare;\
    chmod u+x compare