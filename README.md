# `ma_mapper`
Multiple Alignment Mapper


A python package for mapping genome-wide data onto multiple alignment of transposable elements.

This multiple alignment of THE1C for example:

<img src="docs/img/the1c_aln.png" alt="THE1C alignment" style="width:50%;" />


can be overlayed by genome-wide data such as TF motif prediction, conservation socre, ChIP-exo mapped reads, common base frequencies, alternater allele frequencies, and more:

<p>
  <img src="docs/img/motif.png" alt="THE1C alignment" style="width:30%; display:inline-block;"/>
  <img src="docs/img/bigwig.png" alt="THE1C signal" style="width:30%; display:inline-block;"/>
  <img src="docs/img/bam.png" alt="THE1C signal" style="width:30%; display:inline-block;"/>
  <img src="docs/img/maf.png" alt="THE1C signal" style="width:30%; display:inline-block;"/>
  <img src="docs/img/vcf.png" alt="THE1C signal" style="width:30%; display:inline-block;"/>
</p>


## Quick links
- [Dependencies](#dependencies)
- [Installation](#installation)
- [Package structure](#package-structure)
- [Usage](#usage)
- [Documentaiton quick links](#documentation-quick-links)

## Dependencies
- `python` 3.10
- `biopython` 1.83
- `compress-pickle` 2.1.0
- `cyvcf2` 0.30.28
- `h5py` 3.7.0
- `logomaker` 0.8
- `matplotlib` 3.8.4
- `numpy` 1.21.5
- `pandas` 1.3.5
- `pybedtools` 0.10.0
- `pybigwig` 0.3.23
- `pysam` 0.22.0
- `scipy` 1.7.3

## Installation
To install this package, use the following command to install on a python environment:
```bash
pip install git+https://github.com/pkncsk/ma_mapper@experimental
```
## Package structure
`ma_mapper` is ultimately a wrapper of various external packages, guiding them to work under under the multiple alignment mapping framework. The package itself is separated into modules and submodules based on their main functions on the framework. Notably, the data extraction and data mapping/overlay are streamlined into one module as shown below:

<img src="docs/img/packagestructure.png" alt="THE1C alignment" style="width:70%;" />

## Usage

## Documentation quick links
- [Getting started](docs/gettingstarted.md)
- [Package documentations](docs/package_docs.md)
- [Aligning large TE families](docs/aligning_large_families.md)
## Under construction
- [Tutorials]()
- [API references]()

## License
This package is licensed under the MIT License. 

