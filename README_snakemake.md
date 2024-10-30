
### Installation

Clone the repository from [github](https://github.com/ImagoXV/NanoASV.git)

```git clone https://github.com/ImagoXV/NanoASV.git```

Install the required dependencies with Conda:
```
cd NanoASV
conda env create -f environment.yml
```

Then activate the environment:

```conda activate nanoASV```

### Usage

First download the SILVA database (if needed):

```wget https://www.arb-silva.de/fileadmin/silva_databases/release_138_2/Exports/SILVA_138.2_SSURef_tax_silva.fasta.gz -P resources/```


Ensure the conda environment is active:
```
conda activate nanoASV
```

Test your configuration by performing a dry-run:

```
snakemake -np -s workflow/snakefile --configfile config/config.yaml
```

Execute the workflow (change the number of cores if needed):
```
snakemake -p --cores 4 -s workflow/snakefile --configfile config/config.yaml
```


The workflow can be executed on a cluster using snakemake cluster configuration. Install a [profile](https://github.com/Snakemake-Profiles) for your cluster's job submission system. Edit the defaults in the file `cluster.json` and run the workflow. For example:

```
snakemake -p --jobs 100 --profile slurm --cluster-config cluster.json -s workflow/snakefile --configfile config/config.yaml
```



