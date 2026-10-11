# Getting started
## Quick links
- [About `ma_mapper`](#about-ma-mapper)
- [Dependencies](#dependencies)
- [Installation](#installation)
- [Package structure](#package-structure)
- [Package documentations](#package-documentation)
## Under construction
- [Tutorials](script/documentation/tutorial_nb/00overview.md)
- [API references](script/documentation/api_references/01overview.md)


## About `ma_mapper`
`ma_mapper` is a python package for mapping genome-wide data onto multiple alignment of transposable elements. The main function of this package is to work as a wrapper for other biological data file handlers such as biopython to steamline multiple alignment data overlay. The output of this package is a numerical numpy matrix, which can be easily used for visualization by matplotlib or seaborn. 

## Dependencies
### General package
- `python` 3.10
- `compress-pickle` 2.1.0
- `h5py` 3.7.0
- `numpy` 1.21.5
- `pandas` 1.3.5
- `scipy` 1.7.3

### Genome-wide data file handler
- `biopython` 1.83
- `pybedtools` 0.10.0
- `pybigwig` 0.3.23
- `pysam` 0.22.0
- `cyvcf2` 0.30.28

### Visualization 
- `logomaker` 0.8
- `matplotlib` 3.8.4

## Installation
To install this package, use the following command to install on a python environment:
```bash
pip install git+https://github.com/pkncsk/ma_mapper@experimental
```
## Package structure
`ma_mapper` is ultimately a wrapper of various external packages, guiding them to work under under the multiple alignment mapping framework. The package itself is separated into modules and submodules based on their main functions on the framework. Notably, the data extraction and data mapping/overlay are streamlined into one module as shown below:

<img src="img/packagestructure.png" alt="THE1C alignment" style="width:70%;" />

## Package documentation
[Documentaion](/docs/package_docs.md) for general usages and function parameters.

[Aligning large TE families](/docs/aligning_large_families.md) for families above
~50,000 copies, where the default MAFFT settings run out of memory.

