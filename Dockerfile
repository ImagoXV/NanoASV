#Command to run with Docker 
############"sudo docker run -v $(pwd)/Minimal:/data/Minimal -it nanoasv_dev -d /data/Minimal -o OUTPUT

# Use a base image with Ubuntu and necessary dependencies
FROM ubuntu:22.04

#CMD ["-e", "TZ=$(cat /etc/timezone)"]

RUN touch /etc/localtime
RUN touch /etc/timezone

# Set the timezone
RUN ln -snf /usr/share/zoneinfo/$(cat /etc/timezone) /etc/localtime && echo $(cat /etc/timezone) > /etc/timezone

ENV LANG=C.UTF-8
ENV LC_ALL=C.UTF-8

# # # Install required packages
#  RUN apt-get update && \
#     apt-get install -y \
#     fasttree=2.1.11-2 \
#     gcc=4:13.2.0-7ubuntu1 \
#     mafft=7.505-1 \ 
#     python3=3.12.3-0ubuntu1 \
#     samtools=1.19.2-1build2 \
#     vsearch=2.27.0-1 \
#     wget=1.21.4-1ubuntu4 \
#     bwa=0.7.17-7 \
#     build-essential=12.10ubuntu1 \
#     zlib1g-dev=1:1.3.dfsg-3.1ubuntu2 \
#     python3-pip=24.0+dfsg-1ubuntu1 \
#     r-base=4.3.3-2build2 \
#     parallel=20231122+ds-1 \
#     cowsay=3.03+dfsg2-8 \
#     uuid-runtime=2.39.3-9ubuntu6 \
#     libcurl4-openssl-dev=8.5.0-2ubuntu10.1 \
#     libxml2-dev=2.9.14+dfsg-1.3ubuntu3 \
#     libfontconfig1-dev=2.15.0-1.1ubuntu2 \
#     libharfbuzz-dev=8.3.0-2build2 \
#     libfribidi-dev=1.0.13-3build1
#     #libssl-dev=3.0.13-0ubuntu3 

    # # Install required packages
 RUN apt-get update && \
    apt-get install -y \
    git \
    fasttree \
    gcc \
    mafft \ 
    python3 \
    samtools \
    vsearch \
    wget \
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
    libharfbuzz-dev \
    libfribidi-dev
 #libssl-dev=3.0.13-0ubuntu3 

# Instal minimap2

RUN git clone https://github.com/lh3/minimap2 && \
    cd minimap2/ && make &&\
    mv minimap2 /bin/minimap2

# Install Porechop
RUN wget https://github.com/rrwick/Porechop/archive/refs/tags/v0.2.4.tar.gz && \
    tar -zxvf v0.2.4.tar.gz && \
    rm v0.2.4.tar.gz && \
    cd Porechop-0.2.4 && \
    python3 setup.py install

#Install Chopper
RUN wget https://github.com/wdecoster/chopper/releases/download/v0.7.0/chopper-linux.zip &&  unzip chopper-linux.zip
RUN chmod ugo+rwx chopper && mv chopper /bin/chopper

# Download SILVA138.1
RUN wget https://www.arb-silva.de/fileadmin/silva_databases/release_138_1/Exports/SILVA_138.1_SSURef_tax_silva.fasta.gz && \
    mkdir database && \
    mv SILVA_138.1_SSURef_tax_silva.fasta.gz ./database/ && \
    zcat SILVA_138.1_SSURef_tax_silva.fasta.gz | awk '/^>/ {printf("%s%s\n",(NR==1)?"":RS,$0);next;} {printf("%s",$0);} END {printf("\n");}' | \
     gzip > SILVA_138.1_SSURef_tax_silva.fasta.gz

#RUN echo "Building the index, grab a cup of coffe, it's the longest part" 

#RUN cd database && \
#Changing the indexing step from installation to first run 
#    bwa index -p SILVA_IDX SILVA_138.1_SSURef_tax_silva.fasta.gz && \
#   zgrep "^>" SILVA_138.1_SSURef_tax_silva.fasta.gz | tr -d ">" > Taxonomy_SILVA138.1.csv

RUN apt-get autoremove

# # Install R packages 
RUN R -e 'install.packages("dplyr")'
RUN R -e 'install.packages("BiocManager")'
RUN R -e 'BiocManager::install("phyloseq", dependencies = TRUE)'
#RUN R -e 'BiocManager::install("phyloseq")'

RUN mkdir Rdata

# Copy the script into the container
COPY script.sh /script.sh
COPY script.r /script.r
COPY help.txt /help.txt


# Set the script as the entry pointbwa-mem2-2.2.1_x64-linux
ENTRYPOINT ["/script.sh"]
