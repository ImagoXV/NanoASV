# Use a base image with Ubuntu and necessary dependencies
FROM ubuntu:latest

#CMD ["-e", "TZ=$(cat /etc/timezone)"]

RUN touch /etc/localtime
RUN touch /etc/timezone

# Set the timezone
RUN ln -snf /usr/share/zoneinfo/$(cat /etc/timezone) /etc/localtime && echo $(cat /etc/timezone) > /etc/timezone


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
    r-base

# Install NanoFilt
RUN pip3 install NanoFilt

# Install Porechop
RUN wget https://github.com/rrwick/Porechop/archive/refs/tags/v0.2.4.tar.gz && \
    tar -zxvf v0.2.4.tar.gz && \
    rm v0.2.4.tar.gz && \
    cd Porechop-0.2.4 && \
    python3 setup.py install

# Download SILVA138.1
RUN wget https://www.arb-silva.de/fileadmin/silva_databases/release_138_1/Exports/SILVA_138.1_SSURef_tax_silva.fasta.gz && \
    gunzip SILVA_138.1_SSURef_tax_silva.fasta.gz

# Copy the script into the container
COPY script.sh /script.sh

# Set the script as the entry point
ENTRYPOINT ["/script.sh"]
