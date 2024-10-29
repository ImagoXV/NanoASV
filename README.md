![Logo](NanoASV_logo.png)

# NanoASV
 NanoASV is a container based workflow using state of the art bioinformatic software to process full-length SSU rRNA (16S/18S) amplicons acquired with Oxford Nanopore Sequencing technology. Its strength lies in reproducibility, portability and the possibility to run offline. It can be installed on the Nanopore MK1C  sequencing device and process data locally. 

# Installation
At the moment, the only way to install NanoASV is building it from source with Docker. 
At this point you can whether run it with Docker or to build a Singularity image file (SIF) from the docker version to run with Singularity. 

## ADVANCED - Build from source with Docker
Takes 75 min on my computer (32Gb RAM - 12 cores).
The longest part is SILVA indexing step.
Avoid this step by downloading the (heavy) NanoASV.tar archive
```sh
git clone https://github.com/ImagoXV/NanoASV
docker build -t nanoasv NanoASV/.
```
### Create Docker archive to build with Singularity

```sh
docker save nanoasv -o NanoASV.tar
```
## ADVANCED - Build image with Singularity
I recommend building the sif file from the docker archive 
```sh
singularity build nanoasv docker-archive://NanoASV.tar
```

## NOT WORKING ATM - EASY - Download for Singularity
Archive is too big at the moment to be a GitHub release. You have to build from source.
```sh
wget path/to/archive
tar -xvzf nanoasv.tar.gz 
sudo mv nanoasv /opt/
echo 'export PATH=$PATH:/opt/' >> ~/.bashrc && source ~/.bashrc

```
Then test if everything is working properly
The low vsearch clustering identity threshold allows to successfully recover OTUs from a small number of sequences. 
You should not use such a low identity threshold for analysis. -i 0.7 works fine.

```sh singularity
singularity run nanoasv -d Minimal -o Out_test -i 0.3 [--options]
```
```sh docker
docker run -v $(pwd)/Minimal:/data/Minimal -it nanoasv -d /data/Minimal -o out --docker -i 0.3 [--options] 
```

## ADVANCED - Install on MK1C sequencing device

All previous steps can be used to install on MK1C, but be sure to use the aarch64 version. **IT WILL NOT RUN IF IT'S NOT AARCH64 VERSION**

# Usage
## RECOMMENDED - With Singularity
If added to the path

```sh
nanoasv -d path/to/sequences -o out [--options]
```
Or 
```sh
singularity run nanoasv -d path/to/sequences -o out [--options]
```
Or if installed elsewhere 
```
/path/to/installation/nanoasv -d path/to/sequences -o out [--options] 
```
## ADVANCED - With Docker
I recommand you not to run it with docker because of root privileges.
**Don't forget the --docker flag**
```sh
docker run -v $(pwd)/Minimal:/data/Minimal -it nanoasv -d /data/Minimal -o out --docker
```
You can mount your sequences directory anywhere in the container, but I recommand you to mount in /data/

## Technical recommandations
If running on a PC, I suggest to not use more than two threads with 32Gb of RAM. Otherwise, you might crash your system. 
I highly suggest you to run it on a cluster. 
96 samples (--subsampling 50000) took 4h (without tree) with 150Gb and 8 threads. The tree is highly computer intensive.

## Options

```
| Option               | Description                                                                    |
| -------------------- | ------------------------------------------------------------------------------ |
| `-h`, `--help`       | Show help message                                                              |
| `-v`, `--version`    | Show version information                                                       |
| `-d`, `--dir`        | Path to fastq_pass/                                                            |
| `-db`, `--database`  | Path to reference fasta file                                                   |
| `-q`, `--quality`    | Quality threshold for Chopper, default: 8                                      |
| `-l`, `--minlength`  | Minimum amplicon length for Chopper, default: 1300                             |
| `-L`, `--maxlength`  | Maximum amplicon length for Chopper, default: 1700                             |
| `-i`, `--id-vsearch` | Identity threshold for vsearch unknown sequences clustering step, default: 0.7 |
| `-ab`, `--minab`     | Minimum unknown cluster total abundance to be kept                             |
| `-p`, `--num-process`| Number of cores for parallelization, default: 1                                |
| `--subsampling`      | Max number of sequences per barcode, default: 50,000                           |
| `--no-r-cleaning`    | Flag - to keep Eukaryota, Chloroplast, and Mitochondria sequences              |
|                      | from phyloseq object                                                           |
| `--metadata`         | Specify metadata.csv file directory, default is demultiplexed directory (--dir)|
| `--notree`           | Flag - To remove phylogeny step and subsequent tree from phyloseq object       |
| `--docker`           | Flag - To run NanoASV with Docker                                              |
| `--ronly`            | Flag - To run only the R phyloseq step                                         |
| `--requirements`     | Flag - To display personal reference fasta requirements                        |
```

# How it works 

## Building from source

Building from source is pretty long at the moment.
The main time bottle neck is bwa SILVA138.1 indexing step (~60min on 32Gb RAM PC)
It is way faster if you download the archive and build with Singularity. However, the archive is pretty heavy and not available for download at the moment.

## Data preparation
Directly input your /path/to/sequence/data/fastq_pass directory 
4000 sequences fastq.gz files are concatenated by barcode identity to make one barcodeXX.fastq.gz file.

## Filtering
Chopper will filter for inappropriate sequences.
Is executed in parrallel (default --num-process = 1 )
Default parameters will filter for sequences with quality>8 and 1300bp<length<1700bp

## Chimera detection
<!-- Chimera detection is performed with vsearch --uchime_denovo.
Is executed in parrallel (default --num-process = 6 ) -->
There is no efficient chimera detection step at the moment

## Adapter trimming
Porechop will trimm known adapters 
Is executed in parrallel (default --num-process = 1 )

## Subsampling
50 000 sequences per barcode is enough for most common questions.
Default is set to 50 000 sequences per barcode. 
Can be modified with --subsampling int

## Alignment
bwa will align previously filtered sequences against SILVA 138.1
Can be executed in parrallel (default --num-process = 1 )
50 000 sequences take 14Gb of RAM per thread. Don't overload your computer. With two threads, some RAM will use more than 2x14Gb of RAM. It is important to have whether SWAP partition or to parallelize your computation function of your computer capacities.
In the future, I will add the possibility to use another database than SILVA
barcode*_abundance.tsv, Taxonomy_barcode*.csv and barcode*_exact_affiliations.tsv like files are produced.
Those files can be found in Results directory.

## Unknown sequences clustering
Non matching sequences fastq are extracted then clustered with vsearch (default --id 0.7).
Clusters with abundance under 5 are discarded to avoid useless heavy computing.
Outputs into Results/Unknown_clusters

## Phylogenetic tree generation
Reference ASV sequence from SILVA138.1 are extracted accordingly to detected references. 
Unknown OTUs seed sequence are added. The final file is fed to FastTree to produce a tree file
Tree file is then implemented into the final phyloseq object.
This allows for phylogeny of unknown OTUs and 16S based phylogeny taxonomical estimation of the entity.

## Phylosequization
Alignements results, taxonomy, clustered unknown entities and 16S based phylogeny tree are used to produce a phyloseq opbject: NanoASV.rdata
Please refer to the metadata.csv file in Minimal dataset to be sure to input the correct file format for phyloseq to produce a correct phyloseq object.
You can choose not to remove Eukaryota, Chloroplasta and Mitochondria sequences (pruned by default) using --r_cleaning 0
### --ronly option
Sometimes, your metadata.csv file will not meet phyloseq standards. 
To avoid you recomputing all the previous steps, a --ronly flag can be added. 
Just precise --dir and --out as in your first treatment. NanoASV will find final datasets and run only the r script. 
This will save you time.

## Acknowledgments

We thank Antoine Cousson, Fiona Elmaleh and Meren for their time and energy with NanoASV beta testing !

## Citation
Please don't forget to cite NanoASV and dependencies if it helped you treat your Nanopore data
Thank you !

Dependencies citations :

Danecek, Petr, James K Bonfield, Jennifer Liddle, John Marshall, Valeriu Ohan, Martin O Pollard,
Andrew Whitwham, et al. 2021. “Twelve Years of SAMtools and BCFtools.” GigaScience
10 (2): giab008. https://doi.org/10.1093/gigascience/giab008.

De Coster, Wouter, and Rosa Rademakers. 2023. “NanoPack2: Population-Scale Evaluation of
Long-Read Sequencing Data.” Edited by Can Alkan. Bioinformatics 39 (5): btad311.
https://doi.org/10.1093/bioinformatics/btad311.

Katoh, K., and D. M. Standley. 2013. “MAFFT Multiple Sequence Alignment Software Version 7:
Improvements in Performance and Usability.” Molecular Biology and Evolution 30 (4):
772–80. https://doi.org/10.1093/molbev/mst010.

McMurdie, Paul J., and Susan Holmes. 2013. “Phyloseq: An R Package for Reproducible
Interactive Analysis and Graphics of Microbiome Census Data.” Edited by Michael Watson.
PLoS ONE 8 (4): e61217. https://doi.org/10.1371/journal.pone.0061217.

Md, Vasimuddin, Sanchit Misra, Heng Li, and Srinivas Aluru. 2019. “Efficient Architecture-Aware
Acceleration of BWA-MEM for Multicore Systems.” arXiv. http://arxiv.org/abs/1907.12931.

Nygaard, Anders B., Hege S. Tunsjø, Roger Meisal, and Colin Charnock. 2020. “A Preliminary
Study on the Potential of Nanopore MinION and Illumina MiSeq 16S rRNA Gene
Sequencing to Characterize Building-Dust Microbiomes.” Scientific Reports 10 (1): 3209.
https://doi.org/10.1038/s41598-020-59771-0.

Price, M. N., P. S. Dehal, and A. P. Arkin. 2009. “FastTree: Computing Large Minimum Evolution
Trees with Profiles Instead of a Distance Matrix.” Molecular Biology and Evolution 26 (7):
1641–50. https://doi.org/10.1093/molbev/msp077.

Quast, Christian, Elmar Pruesse, Pelin Yilmaz, Jan Gerken, Timmy Schweer, Pablo Yarza, Jörg
Peplies, and Frank Oliver Glöckner. 2012. “The SILVA Ribosomal RNA Gene Database
Project: Improved Data Processing and Web-Based Tools.” Nucleic Acids Research 41 (D1): D590–96. https://doi.org/10.1093/nar/gks1219.

Rodríguez-Pérez, Héctor, Laura Ciuffreda, and Carlos Flores. 2021. “NanoCLUST: A Species-Level Analysis of 16S rRNA Nanopore Sequencing Data.” Edited by Birol Inanc. Bioinformatics 37 (11): 1600–1601. https://doi.org/10.1093/bioinformatics/btaa900.

Rognes, Torbjørn, Tomáš Flouri, Ben Nichols, Christopher Quince, and Frédéric Mahé. 2016. “VSEARCH: A Versatile Open Source Tool for Metagenomics.” PeerJ 4 (October): e2584. https://doi.org/10.7717/peerj.2584. 

Santos, Andres, Ronny van Aerle, Leticia Barrientos, and Jaime Martinez-Urtaza. 2020. “Computational Methods for 16S Metabarcoding Studies Using Nanopore Sequencing Data.” Computational and Structural Biotechnology Journal 18: 296–305. https://doi.org/10.1016/j.csbj.2020.01.005.





