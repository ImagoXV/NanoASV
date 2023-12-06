# Use a base image with Ubuntu and necessary dependencies
FROM ubuntu:latest

#CMD ["-e", "TZ=$(cat /etc/timezone)"]

RUN touch /etc/localtime
RUN touch /etc/timezone

# Set the timezone
RUN ln -snf /usr/share/zoneinfo/$(cat /etc/timezone) /etc/localtime && echo $(cat /etc/timezone) > /etc/timezone

ENV LANG C.UTF-8
ENV LC_ALL C.UTF-8

# Install required packages
RUN apt-get update && \
    apt-get install -y \
    gcc \
    python3 \
    bwa \
    samtools \
    vsearch \
    wget \
    build-essential \
    zlib1g-dev \
    python3-pip \
    r-base \
    parallel \
    cowsay \
    uuid-runtime

# Install NanoFilt
RUN pip3 install NanoFilt #A changer

# Install Porechop
RUN wget https://github.com/rrwick/Porechop/archive/refs/tags/v0.2.4.tar.gz && \
    tar -zxvf v0.2.4.tar.gz && \
    rm v0.2.4.tar.gz && \
    cd Porechop-0.2.4 && \
    python3 setup.py install

# Install MAAFT
RUN wget https://mafft.cbrc.jp/alignment/software/mafft-7.475-with-extensions-src.tgz && \
    tar -zxvf mafft-7.475-with-extensions-src.tgz && \
    rm mafft-7.475-with-extensions-src.tgz && \
    cd mafft-7.475-with-extensions/core && \
    make && \
    mv mafft /usr/local/bin

# Install FastTree
RUN wget http://www.microbesonline.org/fasttree/FastTree && \
    chmod +x FastTree && \
    mv FastTree /usr/local/bin/

# Download SILVA138.1
RUN wget https://www.arb-silva.de/fileadmin/silva_databases/release_138_1/Exports/SILVA_138.1_SSURef_tax_silva.fasta.gz && \
    gunzip SILVA_138.1_SSURef_tax_silva.fasta.gz

    #mv Silva /data/database/Silva
    echo Ca va etre long

    #bwa-mem index SILVA

# # Install R packages 
# RUN R -e "install.packages('dplyr', repos='http://cran.rstudio.com/')"
# RUN R -e "BiocManager::install('phyloseq')"


# Copy the script into the container
COPY script.sh /script.sh
COPY script.r /script.r

# Set the script as the entry point
ENTRYPOINT ["/script.sh"]
