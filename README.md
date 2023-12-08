# NanoASV
NanoASV is a docker based Nanopore 1500bp 16 Metabarcoding amplicon data analysis workflow. 


# Installation
## Build from source with Docker
```
git clone https://github.com/ImagoXV/NanoASV
docker build -it nanoasv NanoASV/.
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
```
## With Docker
```
sudo docker run -v $(pwd)/Minimal:/data/Minimal -it nanoasv_dev -d /data/Minimal -o OUTPUT
```
## Options
```
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

## Building from source

Building from source is pretty long at the moment.
The main time bottle neck is bwa-meme2 SILVA138.1 indexing step (~30min)
It is way faster if you download the archive and build with Singularity. However, the archive is pretty heavy. 

## Data preparation
Directly input your /path/to/sequence/data/fastq_pass directory 
4000 sequences fastq.gz files are concatenated by barcode identity to make one barcodeXX.fastq.gz file.

## Filtering
NanoFilt will filter for inappropriate sequences. 
Default parameters will filter for sequences with quality>8 1300bp<length<1700bp

## Chimera detection
Ongoing work

## Adapter trimming
Porechop will trimm known adapters 
