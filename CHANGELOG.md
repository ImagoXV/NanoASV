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
