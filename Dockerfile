#Command to run with Docker 
############"sudo docker run -v $(pwd)/Minimal:/data/Minimal -it nanoasv_dev -d /data/Minimal -o OUTPUT

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

# # Use the official Miniconda image as a base
# FROM continuumio/miniconda3:latest
# RUN  conda update -n base -c defaults conda
# # Switch to root user for installation
# USER root

# #Install Chopper
# RUN conda install -c bioconda chopper

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
    mkdir database && \
    mv SILVA_138.1_SSURef_tax_silva.fasta.gz database/SILVA_138.1_SSURef_tax_silva.fasta.gz

    Check if the index exists
RUN echo && \
    if [[ $(ls database/SILVA_138.1_SSURef_tax_silva.amb 2>/dev/null | wc -l) -eq 0 ]]; then && \
    (cd database && \
    # Create the index
    echo Indexing SILVA && \
    date && \
    bwa index SILVA_138.1_SSURef_tax_silva.fasta.gz && \
    zcat SILVA_138.1_SSURef_tax_silva.fasta.gz | grep ">"  | sed 's/.//' > Taxonomy_SILVA138.1.csv) && \
    fi

RUN echo "Building the index, grab a cup of coffe, it's the longest part" 

RUN cd database && \
    bwa index -p SILVA_IDX SILVA_138.1_SSURef_tax_silva.fasta.gz && \
    zcat SILVA_138.1_SSURef_tax_silva.fasta.gz | grep ">"  | sed 's/.//' > Taxonomy_SILVA138.1.csv

# # Install R packages 
# RUN R -e "install.packages('dplyr', repos='http://cran.rstudio.com/')"
# RUN R -e "BiocManager::install('phyloseq')"


# Copy the script into the container
COPY script.sh /script.sh
COPY script.r /script.r

# Set the script as the entry point
ENTRYPOINT ["/script.sh"]
