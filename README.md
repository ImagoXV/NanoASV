![Logo](NanoASV_logo.png)

# NanoASV
NanoASV is a docker based Nanopore 1500bp 16S Metabarcoding amplicon data analysis workflow. 

# Installation
## Download for Singularity
Just download and uncompress the archive then run with singularity accordingly
```sh
wget path/to/archive
singularity run nanoasv -d path/to/sequences -o out [--options]
```

## Build from source with Docker
Takes 75 min on my computer (32Gb RAM - 12 cores).
The longest part is SILVA indexing step.
Avoid this step by downloading the (heavy) NanoASV.tar archive
```sh
git clone https://github.com/ImagoXV/NanoASV
docker build -t nanoasv NanoASV/.
```
### Create Docker archive to build with Singularity

```sh
docker save NanoASV.tar nanoasv
```
## Build image with Singularity

```sh
wget PATH/TO/ARCHIVE
singularity build nanoasv docker-archive://NanoASV.tar
```
# Usage
## With Singularity
```sh
singularity run nanoasv -d path/to/sequences -o out [--options]
```
## With Docker
I highly recommand you not to run it with docker because of root privileges.
Plus, the workflow is Singularity oriented, which means it might not work if running with docker
Don't forget the --docker flag
```sh
docker run -v $(pwd)/Minimal:/data/Minimal -it nanoasv -d /data/Minimal -o out --docker 1
```
## Options

```
| Option    | Description                            |
| --------- | -------------------------------------- |
| `-h`, `--help` | Show help message                 |
| `-v`, `--version` | Show version information       |
| `-d`, `--dir` | path/to/fastq_pass/                |
| `-q`, `--quality | Quality threshold for NanoFilt, default 8
| `-l`, `--minlength` | Minimum amplicon length for Nanofilt, default 1300
| `-L`, `--maxlength` | Maximum amplicon lmength for Nanofilt, default 1700
| `-i`, `--id_vsearch | Identity threshold for vsearch unknown sequences clustering step, default 0.7
| `-p`, `--num_process` | Number of core for parallelization, default = 6
| `--subsampling`, | Max number of sequences per barcode, default 4.10^7
| `--r_cleaning` | logical 0-1 to remove Eukaryota, Chloroplast and Mitochondria sequences from phyloseq object, default 1 (TRUE)
| `--notree` | Flag - To remove phylogeny step and subsequent tree from phyloseq object
| `--docker` | Flag - To run NanoASV with Docker
```

# How it works 

## Building from source

Building from source is pretty long at the moment.
The main time bottle neck is bwa-meme2 SILVA138.1 indexing step (~60min on my computer)
It is way faster if you download the archive and build with Singularity. However, the archive is pretty heavy. 

## Data preparation
Directly input your /path/to/sequence/data/fastq_pass directory 
4000 sequences fastq.gz files are concatenated by barcode identity to make one barcodeXX.fastq.gz file.

## Filtering
Chopper will filter for inappropriate sequences.
Is executed in parrallel (default --num-process = 6 )
Default parameters will filter for sequences with quality>8 1300bp<length<1700bp

## Chimera detection
Chimera detection is performed with vsearch --uchime_denovo.
Is executed in parrallel (default --num-process = 6 )

## Adapter trimming
Porechop will trimm known adapters 
Is executed in parrallel (default --num-process = 6 )

## Subsampling
50 000 sequences per barcode is enough for most common questions.
Default is set to 2.5 millions sequences per barcode. 
Can be modified with --subsampling int

## Alignment
bwa-mem2 will align previously filtered sequences against SILVA 138.1
Is executed in parrallel (default --num-process = 6 )
In the future, I will add the possibility to use another database than SILVA
barcode*_abundance.tsv, Taxonomy_barcode*.csv and barcode*_exact_affiliations.tsv like files are produced.
Those files can be found in Resumlts directory.

## Unknown sequences clustering
Non matching sequences fastq are extracted then clustered with vsearch -id 0.7.
Outputs into Results/Unknown_clusters

## Phylogenetic tree generation
Reference ASV sequence from SILVA138.1 are extracted accordingly to detected references. 
Unknown OTUs seed sequence are added. The final file is fed to FastTree to produce a tree file
Tree file is then implemented into the final phyloseq object.
This allows for phylogeny of unknown OTUs and 16S based phylogeny taxonomical estimation of the entity.

## Phylosequization
Alignements results, taxonomy and clustered unknown entities are used to produce a phyloseq opbject: NanoASV.rdata
In the future, I will add the phylogeny and the tree. But at the moment, the efficient way would be to produce a tree from SILVA database reference sequence. 
Please refer to the metadata.csv file in Minimal dataset to be sure to input the correct file format for phyloseq to produce a correct phyloseq object.
In the future, I will add a possibility to just start from output results if metadata.csv is bad format
You can choose not to remove Eukaryota, Chloroplasta and Mitochondria sequences (pruned by default) using --r_cleaning 0

## Citation
Please don't forget to cite NanoASV and dependencies if it helped you treat your Nanopore data
Thank you !







