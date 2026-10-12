# API references

> **Status note (0.9.0, 2026-10-12).** These pages were written for 0.8.x. Every
> parameter they describe still exists and still means what they say — in two
> places the code has been *fixed to match this documentation*, notably
> `save_to_file`, which is described here as accepting a path and until 0.9.0
> silently ignored one. What is new in 0.9.0 is listed in
> [new in 0.9.0](10new_in_0.9.0.md); what changed in behaviour is in
> [../../CHANGELOG.md](../../CHANGELOG.md) and
> [../CHANGES_FROM_UPSTREAM.md](../CHANGES_FROM_UPSTREAM.md).

## Overview

ma_mapper is a practical realization of multiple alignment mapping framework in TE analysis. The main function is to extract data from genome-wide dataset and map them on TE multiple alignment and visualize the output. 

## Main APIs

Main APIs for this package are: 

- mapper.map_and_overlay(): for data extraction and mapping 
- plots.plot_experimental(): for visualization 

## ma_mapper.mapper.map_and_overlay()

Use TEs coordinate as inputs to extract data from genome=wide datasets and map the extrated data onto multiple alignment, resultig in a data matrix with the same gap structure as multiple alignemnt.

```
from ma_mapper import mapper
alignment_filepath = '/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/output/ma_mapper/hg38_main/alignment/THE1C.fasta.aligned'
genomewide_data_filepath = '/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/resource/annotation/homer_known_motif_hg38/AP-1(bZIP).bed'
output_matrix=mapper.map_and_overlay(alignment_filepath, genomewide_data_filepath,data_format='bed')
```

## ma_mapper.plots.plot()

A matplotlib wrapper that takes a data matrix as an input and plot it in a specific layout, mainly a heatmap of data points with possible annotation and colorbar, accompanied by an aggregated plot under the heatmap.

```
from ma_mapper import plots
plots.plot(data = [output_matrix], heatmap_color=['Greens'], vlim = [[0,0.1]], opacity = 1)
```

## Modules

List of modules:

- ma_mapper
- ma_mapper.mapper
- ma_mapper.plots
- ma_mapper.sequence_alignment
- ma_mapper.extract_bed
- ma_mapper.extract_bam
- ma_mapper.extract_bigwig
- ma_mapper.extract_maf
- ma_mapper.extract_vcf
