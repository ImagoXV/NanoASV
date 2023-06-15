# NanoASV
NanoASv is a docker based Nanopore 1500bp 16 Metabarcoding amplicon data analysis workflow. 


# Installation
## Build from source with Docker
```
git clone
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
| `-h`, `--help` | Show help message                  |
| `-v`, `--version` | Show version information          |
| `-d`, `--dir` | Description of the option           |


```
# How it works 
## Filtering
NanoFilt 
