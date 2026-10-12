# `ma_mapper`
Multiple Alignment Mapper

A Python package for mapping genome-wide data onto a multiple alignment of
transposable elements.

A TE family is hundreds to thousands of copies of the same ancestral sequence,
each mutated, truncated and interrupted differently. Aligning signal to the
family's **consensus** squashes them all onto one idealised coordinate; aligning
it to a **multiple alignment of the copies themselves** keeps every copy as its
own row, and puts each measured value in the column of the base it was measured
on. What comes out is a copies × columns matrix you can slice, sort, split and
plot.

This multiple alignment of THE1C, for example:

<img src="docs/img/the1c_aln.png" alt="THE1C alignment" style="width:50%;" />

can be overlaid with transcription factor motif predictions, conservation
scores, ChIP-exo read pile-ups, cross-species base frequencies, population
allele frequencies, and anything else that comes as BED, BAM, bigWig, VCF or
MAF:

<p>
  <img src="docs/img/motif.png" alt="motif overlay" style="width:30%; display:inline-block;"/>
  <img src="docs/img/bigwig.png" alt="bigWig overlay" style="width:30%; display:inline-block;"/>
  <img src="docs/img/bam.png" alt="BAM overlay" style="width:30%; display:inline-block;"/>
  <img src="docs/img/maf.png" alt="MAF overlay" style="width:30%; display:inline-block;"/>
  <img src="docs/img/vcf.png" alt="VCF overlay" style="width:30%; display:inline-block;"/>
</p>

## Contents

- [Installation](#installation)
- [Quick start](#quick-start)
- [What the matrix means](#what-the-matrix-means)
- [Data formats](#data-formats)
- [Defaults and semantics](#defaults-and-semantics)
- [Package structure](#package-structure)
- [Documentation](#documentation)
- [Citation](#citation)
- [Provenance and licence](#provenance-and-licence)

## Installation

```bash
pip install ma-mapper
```

From a checkout, for development:

```bash
pip install -e ".[dev]"
pytest
```

Python 3.10 or newer. Two external binaries are used when the formats that need
them are used, and neither is a Python dependency:

| binary | needed for | install |
|---|---|---|
| `bedtools` | BED overlays (via pybedtools) | `conda install -c bioconda bedtools` |
| `mafft` | building alignments with `ma-mapper align` | `conda install -c bioconda mafft` |

`ma-mapper info` prints the version of everything installed and says which of
these two it can find.

## Quick start

### From the command line

```bash
# 1. the genomic sequence of every copy, with its coordinates in the header
ma-mapper sequence --repeatmasker hg38.rmsk.tsv --subfamily THE1C \
                   --genome hg38.fa --output THE1C.fasta

# 2. align them
ma-mapper align --input THE1C.fasta --threads 8 --output THE1C.aligned

# 3. map a dataset onto the alignment
ma-mapper overlay --alignment THE1C.aligned \
                  --data hg38.phyloP447way.bw --format bigwig \
                  --output THE1C_phyloP.npy

# 4. look at it
ma-mapper plot --matrix THE1C_phyloP.npy --kind profile --output THE1C_phyloP.png
```

Step 3 writes three files: the matrix, `THE1C_phyloP.rows.tsv` (which copy each
row is) and `THE1C_phyloP.alignment.npy` (which cells are gaps). A matrix
without those two is not interpretable later.

### From Python

```python
from ma_mapper import mapper

matrix, copies, alignment = mapper.map_and_overlay(
    alignment='THE1C.aligned',
    data_file='hg38.phyloP447way.bw',
    data_format='bigwig',
    return_details=True,
)

# mean conservation per alignment column, over the copies that have a base there
profile = mapper.normalise(alignment, matrix, method='average')
```

Each step is also usable on its own: `sequence_alignment.sequence_io` to fetch
sequence, `mapper.parse_and_filter` to get an alignment matrix and its copies,
one of the `extract_*` modules to read a dataset, `mapper.map_data` to put the
two together, and `plots` to draw the result.

## What the matrix means

Row *i* is one copy. Column *j* is one column of the multiple alignment. Cell
(*i*, *j*) is the value measured at the base copy *i* has in column *j*, in the
copy's own 5'->3' orientation. A minus-strand copy is reverse complemented and
its fragments are read right to left, so column 0 is the start of the element
and not the leftmost base in the genome.

A cell where the copy has a gap is 0. A cell where the copy has a base the
dataset does not cover is whatever the dataset returns there, which for bigWig
is `nan` and for BED is 0. That difference matters: see below.

`parse_alignment` codes the alignment itself as 0 gap, 1 A, 2 C, 3 T, 4 G, 5 N,
so `alignment > 0` is the non-gap mask.

## Data formats

| `data_format` | reads | value per base |
|---|---|---|
| `bigwig` | bigWig | the track's value, `nan` where it has no data |
| `bed` | BED (needs `bedtools`) | sum of the scores of the features covering the base; a `.` score counts as 1 |
| `vcf` | tabix-indexed VCF | sum of an INFO field over the SNPs at that base, 0 elsewhere |
| `read_sum`, `read_forward`, `read_reverse`, `read_max`, `read_min` | indexed BAM | number of read 5' ends, forward/reverse relative to the **copy's** strand |
| `normal_*` | indexed BAM | the same ends smoothed with a Gaussian |
| `maf` | MAF + Biopython index | a per-column summary across species: `coverage`, `common_freq`, `ref_freq`, `base_count`, ... |

## Defaults and semantics

Some defaults and definitions decide what a number means. The two that catch
people most often: the BAM `offset` is **5**, not 0, so every read is shifted
five bases downstream, and `row_threshold` is a non-gap **density** within each
copy's own span, which correlates with element age rather than with length.

The full list, with the measured effect of each, is in
[docs/semantics.md](docs/semantics.md). It is worth reading once before
publishing a number from this package.

## Package structure

`ma_mapper` wraps the file-format readers (pysam, pyBigWig, pybedtools, cyvcf2,
Biopython) and gives them a common framework: extract per base, then map onto
alignment columns.

<img src="docs/img/packagestructure.png" alt="package structure" style="width:70%;" />

| module | does |
|---|---|
| `mapper` | parse and filter an alignment, map extracted data onto it, `map_and_overlay` |
| `sequence_alignment` | coordinates to sequence to alignment |
| `extract_bed`, `extract_bam`, `extract_bigwig`, `extract_vcf`, `extract_maf` | one dataset format each |
| `plots` | heatmaps, profiles, annotations, base counts, sequence logos |
| `cli` | the `ma-mapper` command |
| `utility`, `custom_cmap` | RepeatMasker and age tables, colour maps |

## Documentation

- [Getting started](docs/gettingstarted.md)
- [Package documentation](docs/package_docs.md): every function and parameter
- [Tutorial](docs/tutorial.md): a toy family in a minute, then a real one
- [Defaults and semantics](docs/semantics.md): what the defaults mean for a
  number, in one place
- [API reference](docs/api_reference/01overview.md): every function and every
  parameter ([what 0.9.0 added](docs/api_reference/10new_in_0.9.0.md))
- [Notebook tutorials](docs/tutorials/00overview.md): six notebooks written
  against lab paths, which need substituting
- [Aligning large TE families](docs/aligning_large_families.md)
- [Releasing](docs/RELEASING.md)
- [Changes from the original `ma_mapper`](docs/CHANGES_FROM_UPSTREAM.md): what
  this fork changes and what it does not, with the provenance of each layer
- [CHANGELOG](CHANGELOG.md). Read the 0.9.0 entry if you have results from an
  earlier version. Three of its fixes change numbers.

## Citation

See [`CITATION.cff`](CITATION.cff).

## Provenance and licence

`ma_mapper` was written by Pakkanan Chansongkrow and is described in their PhD
thesis (University of Cambridge, 2024). The multiple-alignment overlay approach
it implements comes from the Imbeault lab (Ecco *et al.*, 2016). Later
maintenance, the large-family alignment route, the command line, the test suite
and the 0.9.0 fixes are by Jeremi Jaksina.

MIT licence; see [`LICENSE`](LICENSE).
