# Tutorial

Two walkthroughs. The first needs nothing but the installed package and takes
about a minute; the second is the same thing on a real family, and needs a
genome, a RepeatMasker table and a signal track.

- [1. A family built from nothing](#1-a-family-built-from-nothing)
- [2. A real family](#2-a-real-family)
- [3. Reading the matrix](#3-reading-the-matrix)
- [4. Splitting the copies](#4-splitting-the-copies)
- [5. When the numbers look wrong](#5-when-the-numbers-look-wrong)

## 1. A family built from nothing

```bash
python examples/toy_end_to_end.py /tmp/ma_mapper_tutorial
```

That script writes a 400 bp chromosome, puts four copies of a short element on
it — one on the minus strand, one split into two fragments — and builds a
bigWig whose value at position *p* is *p*. Then it runs the whole chain and
checks the result against arithmetic:

```
fetched 4 copies, lengths [20, 20, 20, 15]
overlay matrix (4, 20), 75 non-gap cells
every value in the column of its own base: True
```

That last line is the point. Because the track's value at each position is the
position itself, a copy's row can be read directly: `copyB` is on the minus
strand at chrT:200-220, so its row must read 219, 218, … 200. If the
orientation, the fragment order or the offset were wrong, the check would fail
rather than produce a plausible-looking figure.

Read the script — it is 120 lines and it is the whole API.

## 2. A real family

Four commands. Substitute your own paths.

```bash
GENOME=hg38.fa
RMSK=hg38.rmsk.tsv          # UCSC rmsk table, with genoName/genoStart/repName
TRACK=hg38.phyloP447way.bw

# 1. the genomic sequence of every copy, with its coordinates in the header
ma-mapper sequence --repeatmasker $RMSK --subfamily THE1C \
                   --genome $GENOME --output THE1C.fasta

# 2. align them. Above ~50,000 copies, read docs/aligning_large_families.md
#    instead of using this step.
ma-mapper align --input THE1C.fasta --threads 8 --output THE1C.aligned

# 3. map the track onto the alignment
ma-mapper overlay --alignment THE1C.aligned --data $TRACK --format bigwig \
                  --row-threshold 0.10 --output THE1C_phyloP.npy

# 4. one value per column
ma-mapper plot --matrix THE1C_phyloP.npy --kind profile \
               --output THE1C_phyloP.png
```

Note `--row-threshold 0.10` rather than the default 0.50. The row threshold is
a *density* — the share of non-gap bases within a copy's own span — so at 0.50
an old, deletion-riddled family loses most of its copies, and the ones it keeps
are the unusually intact ones. That is a selection on exactly the covariate
most TE analyses condition on. For a young, compact family the default is fine.

Step 3 writes four files:

| file | what it is |
|---|---|
| `THE1C_phyloP.npy` | the matrix: one row per copy, one column per alignment column |
| `THE1C_phyloP.rows.tsv` | which copy each row is, with its coordinates |
| `THE1C_phyloP.alignment.npy` | the alignment itself, so a gap is distinguishable from a measured zero |
| `THE1C_phyloP.run.json` | the inputs, their checksums, the parameters and the versions |

## 3. Reading the matrix

```python
import numpy as np, pandas as pd
from ma_mapper import mapper, plots

matrix = np.load('THE1C_phyloP.npy')
alignment = np.load('THE1C_phyloP.alignment.npy')
copies = pd.read_csv('THE1C_phyloP.rows.tsv', sep='\t')

matrix.shape                      # (copies, columns)
alignment > 0                     # True where that copy has a base
np.isnan(matrix).sum()            # bigWig returns nan where it has no data
```

Three distinctions that matter, and which the matrix keeps separate:

- **gap** — the copy has no base in this column: `alignment == 0`, and the
  overlay cell is 0;
- **not measured** — the copy has a base but the track does not cover it:
  `nan` (bigWig) or 0 (BED, where absence genuinely means no feature);
- **measured as zero** — a real 0.0.

Collapsing a column therefore needs a decision about the second case:

```python
# divides by the copies with a base here, nan values included
mapper.normalise(alignment, matrix, method='average')

# divides by the copies that actually have a value here
mapper.normalise(alignment, matrix, method='average', denominator='measured')
```

## 4. Splitting the copies

The reason for keeping one row per copy is to be able to split them. The rows
table is in matrix row order, so any per-copy annotation can be joined onto it
and used as a mask:

```python
bound = copies['name'].isin(copies_with_a_peak).to_numpy()

import matplotlib.pyplot as plt
fig, (top, bottom) = plt.subplots(2, 1, figsize=(8, 4), sharex=True)
plots.plot_heatmap(data=matrix[bound], matplot_axes=top)
plots.plot_heatmap(data=matrix[~bound], matplot_axes=bottom)
```

and the two profiles compared column by column:

```python
import numpy as np
with_peak = mapper.normalise(alignment[bound], matrix[bound], method='average')
without  = mapper.normalise(alignment[~bound], matrix[~bound], method='average')
difference = with_peak - without
print(f'largest difference at column {int(np.nanargmax(np.abs(difference)))}')
```

A column index is not a consensus position. If you need to speak about the
consensus, map it explicitly — do not assume the two frames agree, because for
fragmented or offset copies they do not.

## 5. When the numbers look wrong

| symptom | likely cause |
|---|---|
| `copy N has X non-gap bases in the alignment but Y extracted values` | the alignment and the coordinate table were built from different copy sets, or the aligner trimmed columns (`hmmalign --trim`) or wrote `.` for gaps without translating them to `-` |
| every value is 0 | for BED, no feature overlaps these copies; for a bigWig, check the chromosome naming (`chr1` against `1`) — an unknown name raises rather than returning zeros |
| the whole matrix is gaps | the alignment was passed as an exhausted iterator (fixed in 0.9.0, but check the version with `ma-mapper info`) |
| a profile looks flat where a footprint is expected | the row filter may have removed the copies that carry it: compare `rows.tsv` against the family size |
| read counts stop at 255 | an overlay built before 0.9.0 — the profile was `uint8` and wrapped |
| memory blows up | the overlay canvas is copies × columns × 4 bytes; a 100,000-copy family with 20,000 columns is 8 GB. Filter columns first, or work in blocks |

Two checks worth running before building on a new dataset:

```bash
# re-read the track independently and compare base by base
python examples/verify_overlay.py --coordinates copies.bed \
    --genome hg38.fa --bigwig track.bw --sample 10

# what produced this matrix
python -c "import json;print(json.dumps(json.load(open('THE1C_phyloP.run.json')),indent=1))"
```
