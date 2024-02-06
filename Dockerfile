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

# # Install required packages
 RUN apt-get update && \
    apt-get install -y \
    fasttree \
    gcc \
    mafft \
    python3 \
    samtools \
    vsearch \
    wget \
    bwa \
    build-essential \
    zlib1g-dev \
    python3-pip \
    r-base \
    parallel \
    cowsay \
    uuid-runtime \
    libcurl4-openssl-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libssl-dev \
    libharfbuzz-dev \
    libfribidi-dev

# Install Porechop
RUN wget https://github.com/rrwick/Porechop/archive/refs/tags/v0.2.4.tar.gz && \
    tar -zxvf v0.2.4.tar.gz && \
    rm v0.2.4.tar.gz && \
    cd Porechop-0.2.4 && \
    python3 setup.py install

#Install Chopper
RUN wget https://github.com/wdecoster/chopper/releases/download/v0.7.0/chopper-linux.zip &&  unzip chopper-linux.zip
RUN chmod ugo+rwx chopper && mv chopper /opt/

#Install bwa-mem2
RUN wget https://github.com/bwa-mem2/bwa-mem2/releases/download/v2.2.1/bwa-mem2-2.2.1_x64-linux.tar.bz2 && tar -xf bwa-mem2-2.2.1_x64-linux.tar.bz2
RUN chmod ugo+rwx bwa-mem2-2.2.1_x64-linux/* && mv bwa-mem2-2.2.1_x64-linux /opt/
ENV PATH="${PATH}:/opt/bwa-mem2-2.2.1_x64-linux/"

# Download SILVA138.1
RUN wget https://www.arb-silva.de/fileadmin/silva_databases/release_138_1/Exports/SILVA_138.1_SSURef_tax_silva.fasta.gz && \
    mkdir database && \
    mv SILVA_138.1_SSURef_tax_silva.fasta.gz ./database/ && \
    zcat SILVA_138.1_SSURef_tax_silva.fasta.gz | awk '/^>/ {printf("%s%s\n",(NR==1)?"":RS,$0);next;} {printf("%s",$0);} END {printf("\n");}' | \
     gzip > SILVA_138.1_SSURef_tax_silva.fasta.gz

RUN echo "Building the index, grab a cup of coffe, it's the longest part" 

RUN cd database && \
    #bwa-mem2 index -p SILVA_IDX SILVA_138.1_SSURef_tax_silva.fasta.gz && \
    bwa index -p SILVA_IDX SILVA_138.1_SSURef_tax_silva.fasta.gz && \
    zgrep "^>" SILVA_138.1_SSURef_tax_silva.fasta.gz | tr -d ">" > Taxonomy_SILVA138.1.csv

RUN apt-get autoremove

# # Install R packages 
RUN R -e 'install.packages("dplyr")'
RUN R -e 'install.packages("BiocManager")'
#RUN R -e 'install.packages("tidyverse")'
RUN R -e 'BiocManager::install("phyloseq", dependencies = TRUE)'

RUN mkdir Rdata

# Copy the script into the container
COPY script.sh /script.sh
COPY script.r /script.r

# Set the script as the entry pointbwa-mem2-2.2.1_x64-linux
ENTRYPOINT ["/script.sh"]
