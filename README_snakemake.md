
### Installation

Clone the repository from [github](https://github.com/ImagoXV/NanoASV.git)

```git clone https://github.com/ImagoXV/NanoASV```

Install the required dependencies with Conda:
```
cd NanoASV
conda env create -f environment.yml
ACTIVATE_DIR=$(conda info --base)/envs/NanoASV/etc/conda/activate.d
cp config/alias.sh $ACTIVATE_DIR/
cp config/paths.sh $ACTIVATE_DIR/
echo "export NANOASV_PATH=$(pwd)" >> $ACTIVATE_DIR/paths.sh
DEACTIVATE_DIR=$(conda info --base)/envs/NanoASV/etc/conda/deactivate.d
cp config/unalias.sh $DEACTIVATE_DIR/
chmod ugo+x workflow/run.sh

```

Then activate the environment:

```conda activate nanoASV```

### Usage

First download the SILVA database (if needed):

```wget https://www.arb-silva.de/fileadmin/silva_databases/release_138_2/Exports/SILVA_138.2_SSURef_tax_silva.fasta.gz -P resources/```


Test your configuration by performing a dry-run:

```
snakemake -np -s workflow/snakefile --configfile config/config.yaml
```

Execute the workflow (change the number of cores if needed):
```
snakemake -p --cores 4 -s workflow/snakefile --configfile config/config.yaml
```




