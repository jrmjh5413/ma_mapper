# API references

> **Verified against 0.9.0 on 2026-10-12.** Pages 02 to 09 were written for
> 0.8.x. Every signature on them is now compared with the code by
> `tests/test_docs_match_code.py`, so a parameter cannot be added, removed or
> given a new default without a test failing. That check found, and this
> revision fixes: two functions that had gained parameters, one documented
> parameter that does not exist (`map_data(nested_data=...)`), four BAM
> filtering parameters that were never documented, a counting mode list that
> advertised four modes and a default returning no data, and two function
> names that could not be imported (`seqeunce_io`,
> `ma_mapper.coordinate_table_to_bed`).
>
> In two places the code was changed to match what these pages already said,
> rather than the other way round: `save_to_file` accepting a path, which until
> 0.9.0 was silently ignored, and the default output directory, which was the
> installed package directory and is now the current one.
>
> Parameters added since 0.8.x are in [new in 0.9.0](10new_in_0.9.0.md).
> Behaviour that changed is in [../../CHANGELOG.md](../../CHANGELOG.md) and
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
