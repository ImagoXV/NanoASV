# NanoASV
NanoASv is a docker based Nanopore 1500bp 16 Metabarcoding amplicon data analysis workflow. 


# Installation
## Build from source with Docker
```
git clone --filter=NanoASV.tar:none https://github.com/ImagoXV/NanoASV
docker build -t nanoasv NanoASV/.
```


## Build image with Singularity
```
wget PATH/TO/ARCHIVE
singularity build --sandbox nanoasv docker-archive://NanoASV.tar
```
# Usage
## With Singularity
```
singularity exec nanoasv workflow -d path/to/sequences [--options]

| Option    | Description                            |
| --------- | -------------------------------------- |
| `-h`, `--help` | Show help message                 |
| `-v`, `--version` | Show version information       |
| `-d`, `--dir` | Description of the option          |
| `-q`, `--quality | Quality threshold for NanoFilt, default 8
| `-l`, `--minlength` | Minimum amplicon length for Nanofilt, default 1300
| `-L`, `--maxlength` | Maximum amplicon lmength for Nanofilt, default 1700
| `-i`, `--id_vsearch | Identity threshold for vsearch unknown sequences clustering step, default 0.7
| `-p`, '--num_process` | Number of core for parallelization, default = 6



```
# How it works 
## Data preparation
Directly input your /path/to/sequence/data/fastq_pass directory 
4000 sequences fastq.gz files are uncompressed and concatenated by barcode identity

## Filtering
NanoFilt will filter for inappropriater sequences. 
Default parameters will filter for sequences with quality>8 1300bp<length<1700bp

## Chimera detection
Ongoing work

## Adapter trimming
Porechop will trimm known adapters 
