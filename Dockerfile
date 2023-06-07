# Use a base image with Ubuntu and necessary dependencies
FROM ubuntu:latest

# Install required packages
RUN apt-get update && \
    apt-get install -y \
    gcc \
    python3 \
    bwa \
    samtools \
    vsearch \
    r-base

# Install R packages
RUN R -e "install.packages(c('ggplot2', 'phyloseq', 'vegan'), repos = 'https://cloud.r-project.org/')"

# Copy the script into the container
COPY script.sh /script.sh

# Set the script as the entry point
ENTRYPOINT ["/script.sh"]
