# Changelog

All notable changes to NanoASV are documented in this file. The
format is based on [Keep a Changelog](https://keepachangelog.com/),
and NanoASV follows [Semantic Versioning](https://semver.org/).
Earlier versions are described in the
[GitHub releases](https://github.com/ImagoXV/NanoASV/releases).


## [Unreleased]

### Upgrade notes

- Results change: reads aligned to the reverse strand are now counted,
  so abundances roughly double. On the mock dataset, assigned reads go
  from 3,000 to 6,043 and taxa from 1,700 to 2,714.
- Runs take longer: the extra taxa (above) make the phylogeny larger,
  and FastTree runs on a single thread whatever `--num-process`. On the
  mock dataset (12-core computer), FastTree goes from 44 s to 85 s, and
  the whole run from 174 s to 224 s with `nanoasv --mock` (2 cores),
  from 122 s to 168 s with `--num-process 12` (#213). `--notree` skips
  the phylogeny

### Added

- `Results/Read_counts.tsv`: number of reads of each sample after each
  step, including the reads discarded because their mapping quality is
  below `--sam-qual` (#145)

### Changed

- Messages of the pipeline tools (minimap2, chopper, vsearch, mafft,
  FastTree...) go to log files in `tmp_files/logs/` instead of the
  terminal; when a step fails, Snakemake prints the path of its log.
  The R step still writes to the terminal (and to its log)
- Error and warning messages of `nanoasv` (and the help or requirements
  text printed with them) go to stderr instead of stdout
- Interleaved reference files are converted to one line per sequence
  with vsearch: same result for valid files, but alignment characters
  (`-`, `.`) now stop the run with an error, and entries without a
  sequence are dropped
- `--fasttree-model` no longer accepts `wag` and `lg`: these protein
  models had no effect on nucleotide sequences (same tree as the
  default model)

### Removed

- Support for Oxford Nanopore MK1C devices, discontinued by Oxford
  Nanopore: README installation instructions removed.

### Fixed

- Reads aligned to the reverse strand (about half of the reads) were
  dropped after alignment (#154)
- Reference IDs are matched exactly: an ID that was a prefix of another
  (`AB1.1.15` and `AB1.1.1500`) added references without reads to
  `Results/Tax/` files and to the tree, which could change its
  topology. Abundances were not affected.
- A run started in the same directory as an earlier run picked up the
  earlier run's samples, tree references and unknown clusters from
  `tmp_files/` (the final R step then failed on samples missing from
  `metadata.csv`)
- Unmatched reads of samples not named `barcode*` were left out of the
  unknown-sequence clustering
- Uncompressed fastq files given directly in `-d` (no barcode
  directories) were not processed, and `nanoasv` exited with status 0
- `--dry-run` replaced `-d`, `-o` and `-db` given before it with the
  mock dataset values
- `--subsampling` accepted any value: a non-number or `0` kept no read.
  It must now be a positive integer
- `metadata.csv` errors were followed by a generic "Something went
  wrong" message suggesting a NanoASV bug
- Reference headers with a space after `>` (`> ID taxonomy`) made
  minimap2 crash; they are now rejected before the run starts, and
  `nanoasv --requirements` no longer shows one as a valid example
- Input, output, reference and installation paths containing spaces
  made the run fail
- When some fastq files were uncompressed, samples whose directory is
  not named `barcode*` were lost ("sequences not found")
- In scripts, cluster jobs and CI (stdin not a terminal), the error
  messages for a missing `-d`, `-o` or `metadata.csv` were lost, and
  `nanoasv` could wait as long as stdin stayed open
- A sample name containing an underscore (e.g. `soil_A`) made the final
  R step fail
- A sample name with a dash or a leading digit (e.g. `soil-A`, `2soil`)
  lost its unknown-cluster counts, and was renamed in the CSV header
  (`soil.A`, `X2soil`)
- Input made only of sample directories not named `barcode*` was
  rejected ("No barcode directories or fastq files found")


## [1.2.3] - 2026-10-27

Bug-fix release: no new features. Results are unchanged for runs that
worked with v1.2.2 (checked against the mock dataset).

### Upgrade notes

- The conda environment is now created inside the repository and
  activated by path:
  `conda activate "${HOME}/NanoASV/NanoASV-v1.2.3"`
- `-p/--num-process` is now honoured, with its documented default of 1.
  Runs without `-p` use a single core (v1.2.2 used all CPUs): pass
  `-p N` to go faster.
- `metadata.csv` must list every sample directory. A missing one used
  to make the final R step fail; it is now rejected upfront.

### Fixed

- Installation: `unalias.sh` installed in `deactivate.d`, paths with
  spaces, stop on the first error (#157, #159)
- `nanoasv` no longer runs `conda init` (which modified shell startup
  files) nor `conda activate NanoASV` (#169)
- Gzipped reference files: headers are read after decompression; IDs
  are added in about a second instead of about 20 minutes on SILVA,
  without `uuidgen` (#164)
- Uncompressed interleaved reference files were replaced by an empty
  file (#165)
- `--fasttree-model` accepted any value; error messages now show the
  default and available models (#166)
- `metadata.csv` barcodes were looked for in the `--metadata` directory
  instead of the input directory (#167)
- Porechop, chopper, minimap2 and Snakemake ignored `-p` (#66, #168)
- Python 3.12 `SyntaxWarning: invalid escape sequence '\+'` (#107, #112)
- R errors are no longer hidden; R libraries are looked for in the
  conda environment

### Documentation

- `help.txt` listed `--sam-quality`: the option is `--sam-qual`.
  `--model` and `--fasttree-model` are now documented (#170)
- README: activation path, SILVA 138.2 download, mock output location,
  stale `--r_cleaning` option, broken link, typos (#171, #172)
- `CITATION.cff` cites the Bioinformatics paper (#156)

### Maintenance

- CI: installation, command-line tests, mock run compared with expected
  outputs, lint report (#160, #161, #162, #163)
- Removed the unused `config/config.yaml` and `UV_transition/` (#173)

[Unreleased]: https://github.com/ImagoXV/NanoASV/compare/v1.2.3...dev
[1.2.3]: https://github.com/ImagoXV/NanoASV/compare/v1.2.2...v1.2.3
