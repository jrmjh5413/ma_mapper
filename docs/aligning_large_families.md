# Aligning large TE families

## Quick links
- [The problem](#the-problem)
    - [Matrix footprint](#matrix-footprint)
- [What ma_mapper requires of an alignment](#what-ma_mapper-requires-of-an-alignment)
- [Choosing an approach](#choosing-an-approach)
    - [Alternatives, and when they are worth it](#alternatives-and-when-they-are-worth-it)
- [Profile alignment with `hmmalign`](#profile-alignment-with-hmmalign)
- [De novo alignment with FAMSA](#de-novo-alignment-with-famsa)
- [Row ordering (the second O(n^2) wall)](#row-ordering-the-second-on2-wall)
- [What not to use](#what-not-to-use)
- [Overlay output sizes](#overlay-output-sizes)
- [Benchmark results](#benchmark-results)
- [How far you can trust the accuracy numbers](#how-far-you-can-trust-the-accuracy-numbers)
- [Provenance](#provenance)

## The problem

`sequence_alignment.mafft_align()` calls MAFFT with `--treeout --reorder` and no
algorithm flags, which selects FFT-NS-2. That is fine up to a few tens of
thousands of copies and fails above roughly 65,000.

The failure is **memory**, not time. Peak RSS scales as about `n^1.19`:

| copies | peak RSS |
|---|---|
| 5,000 | 1.60 GB |
| 20,000 | 8.41 GB |
| 223,029 (extrapolated) | ~150 GB |

An `icelake-himem` node provides ~107 GB at 16 cores, so the largest human TE
families cannot complete. This is the sole cause of the families listed in
`missing_alignments.txt` on the genome-wide run: `AluJb`, `AluJr`, `AluSx`,
`AluSx1`, `AluSz`, `AluY`, `L1MB7`, `L1ME1`, `L2`, `L2a`, `L2b`, `L2c`, `L2d2`,
`MIR`, `MIR3`, and `MIRb`.

A second, subtler issue affects every de novo aligner: **alignment width grows
with copy count**, because each additional copy contributes private insertions
that add columns without adding homologous positions. Measured on real outputs:

| family | consensus | copies | columns | inflation |
|---|---|---|---|---|
| MER41B | 634 bp | 2,824 | 6,558 | 10.3x |
| MIR1_Amn | 230 bp | 11,727 | 4,623 | 20.1x |
| AluSq2 | 313 bp | 50,953 | 10,091 | 32.2x |
| AluJo | 312 bp | 45,217 | 19,209 | 61.6x |
| L3 | 4,099 bp | 45,083 | 19,815 | 4.8x |
| L1MB3 | 6,183 bp | 17,223 | 78,918 | 12.8x |

Inflation is worst for short families, mildest for long ones, which is what you
expect if each copy contributes a roughly fixed number of private columns.
Growth is sub-linear in copy count, about `n^0.27` for FAMSA and `n^0.44` for
MAFFT, so it is a nuisance rather than a wall. It does decide how much RAM the
overlay needs. See [Matrix footprint](#matrix-footprint) below.

### Matrix footprint

`map_data()` holds three arrays over the full alignment simultaneously: the
parsed alignment, the coded matrix, and the overlay canvas. Peak memory is
therefore several times the cell count.

```python
parsed_alignment = np.zeros(shape, dtype=np.uint8)     # 1 byte/cell
alignment_matrix = np.zeros(shape, dtype=np.uint8)     # 1 byte/cell
canvas           = np.zeros(shape, dtype=np.float32)   # 4 bytes/cell
```

Peak for one overlay track over the full MIRb alignment:

| aligner | columns | cells | matrix (uint8) | **peak in `map_data`** |
|---|---|---|---|---|
| `mafft --retree 1 --parttree` | 25,365 | 5.66 G | 5.66 GB | **33.9 GB** |
| `mafft --dpparttree --partsize 50` | 22,313 | 4.98 G | 4.98 GB | **29.9 GB** |
| `famsa -medoidtree` | 13,201 | 2.94 G | 2.94 GB | 17.6 GB |
| `kalign --mode fast` | 7,634 | 1.70 G | 1.70 GB | 10.2 GB |
| `hmmalign --outformat afa` | 4,596 | 1.03 G | 1.03 GB | **6.1 GB** |
| `famsa` (default) | 3,750 | 0.84 G | 0.84 GB | **5.0 GB** |

`icelake-himem` provides ~107 GB at 16 cores and plain `icelake` ~54 GB, so every
option now fits on an ordinary node.

> **Changed 2026-08-18.** `alignment_matrix` was `np.ndarray(shape)`, float64
> for data whose only values are 0-5, and `canvas` was float64. Peak was
> therefore 96.2 GB for MAFFT PartTree and 84.6 GB for `--dpparttree`, i.e.
> himem-only with no headroom. Casting them to `uint8` and `float32` cut peak by
> ~65%. This also repaired `base_count()`, which calls
> `np.bincount(col, minlength=5)` and raises `TypeError` on float input.

After the standard `>=10%` occupancy filter the surviving matrix is 244 to 287
columns, i.e. 54 to 64 MB as `uint8`, for every method. The filter runs
*after* `map_data()`, so peak memory is set by the unfiltered width.

## What ma_mapper requires of an alignment

Two constraints decide which aligners are usable. Both fail silently rather than
loudly, so they are worth stating explicitly.

### 1. The per-base invariant

`mapper.map_data()` performs a positional broadcast:

```python
canvas[idx, row > 0] = extracted_data[idx]
```

`row > 0` selects the non-gap positions of alignment row `idx`, and
`extracted_data[idx]` holds the per-base signal extracted from the genome across
that copy's interval. Therefore:

> **For every copy, the number of non-gap characters in its alignment row must
> equal the number of genomic bases in that copy.**

Any aligner, flag, or post-processing step that deletes residues breaks the
overlay. Measured on 5,000 MIRb copies:

| method | bases lost | usable |
|---|---|---|
| `famsa` (default) | 0.00% | yes |
| `mafft --parttree` | 0.00% | yes |
| `hmmalign --outformat afa` (no `--trim`) | 0.00% | yes |
| `famsa -remove-rare-columns 0.05` | 2.10% | **no** |
| `hmmalign --trim`, match-columns only | 10.77% | **no** |

**Corollary:** width reduction must happen *after* mapping, never inside the
aligner. `mapper.create_filter()` is safe precisely because it runs downstream of
`map_data()`, because it hides columns rather than destroying the coordinate
correspondence. Aligner-side column trimming (`-remove-rare-columns`,
`hmmalign --trim`) is not equivalent and must not be substituted for it.

### 2. The parser character set

`mapper.parse_alignment()` maps characters through:

```python
hash_table = {45: 0, 46: 0, 65: 1, 67: 2, 84: 3, 71: 4, 78: 5}  # '-' '.' A C T G N
alignment_matrix = np.zeros(parsed_alignment.shape, dtype=np.uint8)
```

Anything outside that set becomes a gap, and a warning naming the offending
characters is logged. Sequence case is safe, because `.upper()` is applied
before lookup. IUPAC ambiguity codes (R, Y, S, W, ...) are *not* in the table and will
be counted as gaps, with a warning; add them explicitly if a data source emits
them.

> **Changed 2026-08-18.** `alignment_matrix` was allocated with
> `np.ndarray(shape)`, which does not initialise memory, and `.` (ASCII 46) was
> not in `hash_table`. Any unmapped character therefore retained whatever was in
> that memory, making the `row > 0` test in `map_data()` non-deterministic.
> That silently corrupted overlays rather than raising. This mattered for HMMER
> output specifically, where `.` is ~73% of all gap characters. FAMSA and MAFFT
> emit only `-`, which is why it never surfaced.

## Choosing an approach

**Default: MAFFT PartTree.** MAFFT is the reference implementation for this
project. It has proper nucleotide scoring parameters (FAMSA does not; see
below), and it is what published work in this area uses.

```bash
mafft --retree 1 --maxiterate 0 --nofft --parttree --thread 16 \
      all_copies.fasta > family.afa
```

Full MIRb (223,029 copies): 36 m 26 s, 14.7 GB, 25,365 columns, zero bases lost.

The MAFFT modes, all lossless, all on the full 223,029 copies:

| mode | wall | RSS | width | SP | overlay peak* |
|---|---|---|---|---|---|
| `--retree 1 --parttree` | 36 m | 14.7 GB | 25,365 | 0.199 | 33.9 GB |
| **`--dpparttree --partsize 50`** | 1 h 11 m | 12.4 GB | 22,313 | **0.281** | 29.9 GB |
| `--parttree --partsize 1000` | 4 h 43 m | 26.5 GB | **20,101** | 0.194 | 26.9 GB |
| `--dpparttree --partsize 1000` | 5 h 35 m | 13.0 GB | 24,448 | 0.215 | 32.7 GB |

\* after the dtype fix in [Matrix footprint](#matrix-footprint).

`--dpparttree --partsize 50` is the best balance: the most accurate MAFFT mode by
a clear margin (0.281 vs 0.199), narrower than `--retree 1`, at double the
runtime. Prefer it unless 36 minutes vs 71 minutes matters.

This needs the dtype fix in [Matrix footprint](#matrix-footprint). Without it
the overlay peaks at 96 GB, with it 34 GB.

Two things MAFFT does **not** solve, both covered below:

- **Row ordering.** MAFFT's `--treeout` guide tree scores 9% of scipy UPGMA and
  cannot drive the clustering. See
  [Row ordering](#row-ordering-the-second-on2-wall).
- **Width.** 25,365 columns is 6.8x FAMSA's and 5.5x `hmmalign`'s. That is a
  memory cost, not a correctness one. The alignment is lossless and the
  post-filter matrix is the same ~250 columns either way.

### Alternatives, and when they are worth it

| | MAFFT PartTree | FAMSA | `hmmalign` + tree |
|---|---|---|---|
| nucleotide scoring model | **yes** | **no** (protein matrices) | yes |
| wall (223,029) | 36 m | 6 m 41 s | 2 m 25 s |
| width | 25,365 | 3,750 | 4,596 |
| overlay peak (after dtype fix) | 33.9 GB | 5.0 GB | 6.1 GB |
| alignment accuracy (SP) | 0.199 | 0.171 | **0.452-0.513** |
| columns mean consensus positions | no | no | **yes** |
| produces a usable guide tree | no (9%) | **yes (48%)** | no |

**FAMSA** is 5x faster and 6.8x narrower but scores DNA through protein
substitution matrices (`-sm` offers only MIQS and PFASUM31/43/60). That is a
real methodological weakness and a reason to prefer MAFFT for anything
reported. FAMSA remains the best source of a *guide tree* (see below), a role in
which the scoring-matrix objection does not apply. The tree only decides the
order rows are drawn in and enters no reported number.

**`hmmalign`** is the most accurate and gives consensus-anchored columns, which
is what makes matrices comparable across families. It cannot produce a row
ordering. Consider it when cross-family comparability is the point.

**MAFFT's own profile modes do not work here**. see
[What not to use](#what-not-to-use).

## Profile alignment with `hmmalign`

Align every copy independently against a profile HMM of the family rather than
aligning copies to each other. This is `O(n)` in both time and memory, has no
guide tree and no all-vs-all distance matrix, is chunkable across array tasks,
and yields match columns that correspond to consensus positions, so matrices
become comparable *between* families, which no de novo MSA provides.

```bash
# 1. seed: ~1,000 held-out copies through a high-accuracy aligner   (~5 min)
mafft --thread 16 --maxiterate 1000 --localpair seed1000.fasta > seed.afa

# 2. build the profile
hmmbuild --dna --cpu 16 --informat afa family.hmm seed.afa

# 3. align all copies. NO --trim. afa format.                       (~2.5 min)
hmmalign --dna --outformat afa family.hmm all_copies.fasta > family.afa
```

Full MIRb (223,029 copies): **2 min 25 s, 3.13 GB, 4,596 columns, zero bases
lost**; 287 columns after the `>=10%` occupancy filter (64 MB matrix).

**`--trim` must not be used, and the match-column-only projection must not be
used.** They produce a fixed-width alignment (286 to 334 columns regardless of copy
count) which looks ideal but achieves it by discarding 10.77% of bases, exactly
the data ma_mapper exists to map. The fixed-width form is only appropriate for
analyses that do not need per-base fidelity.

HMMER's `.` gap character is handled by `parse_alignment()` since 2026-08-18; no
`tr '.' '-'` preprocessing is needed (it remains harmless if present).

Recovering consensus coordinates from the `afa` output is straightforward: HMMER
applies its own case convention regardless of input case, so match-state columns
are entirely uppercase and insert columns entirely lowercase, with no mixed
columns. On 500 MIRb copies this gives 271 pure-uppercase and 549 pure-lowercase
columns.

Seed choice matters more than seed size:

| seed | model length | SP vs L-INS-i |
|---|---|---|
| L-INS-i, 2,000 copies | 305 | 0.5132 |
| L-INS-i, 1,000 copies | 286 | n/a |
| L-INS-i, 500 copies | 276 | 0.4960 |
| FAMSA, 5,000 copies | 334 | 0.4519 |
| FAMSA, 2,000 copies | 283 | 0.4447 |
| FAMSA, 500 copies | 259 | 0.3351 |

An L-INS-i seed of 500 copies beats a FAMSA seed of 5,000. Returns saturate by
~1,000 copies, which costs about five minutes per family.

## De novo alignment with FAMSA

The simplest complete route: one tool for both the matrix and the guide tree.
Also the right choice when a profile is unwanted: subfamily structure,
phylogeny, or any analysis where imposing a single consensus model would beg the
question. Use default settings:

```bash
famsa -t 16 all_copies.fasta family.afa
```

Full MIRb: 6 min 41 s, 13.4 GB, 3,750 columns, zero bases lost. Request ~25 GB.

Do **not** use `-medoidtree` at this scale. It is faster on subsamples but at
223,029 copies it is 2.6x slower and 3.5x wider than plain FAMSA; its only
advantage is memory (1.21 GB vs 13.4 GB). Extrapolating from 5k/20k subsamples
alone would select it wrongly.

**Caveat:** FAMSA 2.4.1 has no nucleotide scoring model. `-sm` offers only MIQS
and PFASUM31/43/60, all protein substitution matrices. It performs well
empirically on TE DNA (mean pairwise identity 0.459 against L-INS-i's 0.466) but
this should be disclosed in any write-up.

## Row ordering (the second O(n^2) wall)

Producing the alignment is only half the problem. The overlay pipeline orders
rows by hierarchical clustering (`bulk_produce_mapper_overlays.smk`,
`_parse_filter_cluster`):

```python
dist_matrix    = pdist(alignment_matrix, metric="hamming")
linkage_matrix = linkage(dist_matrix, method="average")
row_order      = leaves_list(linkage_matrix)
```

`pdist` is O(n^2) and builds the full condensed distance matrix in memory:

| rows | condensed distance matrix (float64) |
|---|---|
| 20,000 | 1.6 GB |
| 50,000 | 10.0 GB |
| 65,000 | 16.9 GB |
| 100,000 | 40.0 GB |
| 223,029 (MIRb) | **199 GB** |

So this step fails at almost exactly the copy count where the old MAFFT
alignments failed. Fixing the aligner alone does not make a large family
processable. Both walls sit at ~65,000 copies and both must be addressed.

Note also that `row_order` is computed on the **filtered** matrix and applied as
`data[mask][row_order]`. An externally supplied ordering covers all copies
including those `mask` drops, so it must be subset to the surviving rows before
use, or the ordering and the data silently desynchronise.

### Scalable alternatives

Measured on 8,000 MIRb rows (244 columns after filtering) as mean Hamming
distance between adjacent rows. Lower means similar copies are adjacent, which
is what makes the heatmap band. Scored as the percentage of scipy UPGMA's
improvement over a random ordering that each method recovers:

| ordering | % of UPGMA | cost at 223,029 copies |
|---|---|---|
| random baseline | 0% | - |
| `mafft --parttree --treeout` | **9%** | free with the alignment |
| sort by similarity to consensus | 12% | 0.1 s |
| `famsa` default tree (`-gt sl`) | 21% | 3 m 02 s |
| `famsa -gt upgma` | 42% | - |
| `famsa -gt upgma -medoidtree` | 46% | ~10 s |
| `famsa -gt sl -medoidtree` | 46% | 9.3 s |
| **`famsa -gt nj -medoidtree`** | **48%** | **~10 s** |
| scipy UPGMA | 100% | O(n^2) - infeasible above ~65,000 |

```bash
famsa -t 16 -gt nj -medoidtree -gt_export all_copies.fasta family.nwk
```

Leaf labels are truncated at the first `:`, so they come out as the bare
internal ID (`MIRb_181684`), which is unique per copy and matches the key
ma_mapper uses. Map tree order onto matrix rows by ID, never by position.

**Findings that matter for the choice:**

- **MAFFT's PartTree tree is not usable for clustering (9%)**, worse than
  sorting by distance to consensus. `--treeout` does work with `--parttree` and
  costs nothing, but PartTree builds its topology precisely *by not* computing
  pairwise distances, so it partitions well enough to guide a progressive
  alignment and not well enough to order rows.
- **Sorting by similarity to consensus recovers only 12%.** It is a
  one-dimensional sort, so it separates copies by *how much* they diverged, not
  by *how*: subfamilies interleave instead of banding. It is a legitimate
  ordering if an age/divergence gradient is what the figure should show, but it
  is not a clustering.
- **`-medoidtree` is essential and, for ordering, it is also better.** Single
  linkage chains, giving a poor leaf order (21%); medoid trees give compact
  clusters (46-48%) 20x faster.
- **Never use `-gt nj` without `-medoidtree`**. Neighbour-joining is O(n^3) and
  took 96 minutes on 20,000 copies against 1.15 s with `-medoidtree`.
- **`hmmalign` cannot produce a tree at all.** Sequences never interact, so
  there is nothing to emit. If you align with `hmmalign`, source the ordering
  from a separate `famsa -gt_export` run; the two are independent and the tree
  costs ~10 s.

**No scalable method exceeds ~48% of scipy UPGMA's banding.** Below ~65,000
copies `pdist` still fits and UPGMA is strictly better, so use it. Above that
the choice is between a half-quality ordering and none.

## What not to use

**`--partsize 1000`**: poor value rather than useless. `--parttree
--partsize 1000` gives the narrowest MAFFT alignment (20,101 columns, 21% below
`--retree 1`) but takes 4 h 43 m against 36 m, and buys no accuracy (SP 0.194 vs
0.199). `--dpparttree --partsize 1000` is worse on every axis than
`--dpparttree --partsize 50`. If width matters more than time it is defensible;
otherwise use the default partsize.

**MAFFT's profile modes, `--add` and `--addfragments`.** The natural idea, to
align every copy to a seed alignment as `hmmalign` does, does not work in MAFFT:

- `--addfragments` **segfaults** in `pairlocalalign` (4 m 50 s in, 20,000 copies,
  1,000-sequence L-INS-i seed).
- `--add` completes but produces **25,867 columns from only 20,000 copies** in
  1 h 01 m, more than the full 223,029-copy de novo PartTree alignment. Insertions
  relative to the seed each open their own columns, and the seed's own width
  (1,511 columns for an L-INS-i alignment of 1,000 copies) propagates.

Use `hmmalign` if a profile alignment is wanted; MAFFT has no working equivalent
at this scale.

**Kalign 3.6.** Disqualified on quality, not speed. Mean pairwise identity within
aligned columns is 0.276, against 0.466 for L-INS-i and 0.459 for FAMSA, so it
stacks non-homologous bases. Its `default` mode emits 25,169 columns for 500
sequences at 99.4% gaps, i.e. effectively unaligned. Its speed (2 min 10 s at
full scale) comes from not doing the work.

## Overlay output sizes

`map_data()` returns the matrix already filtered (`mapped_data[np.ix_(row_filter,
col_filter)]`), so overlay outputs are ~250 columns rather than the full
alignment width. The width problem is therefore an in-RAM problem only; it does
not reach disk. What does reach disk is one matrix per data track, written by the
overlay pipeline as CSV.

Measured on MER41B (1,589 rows x 661 columns = 1.05 M cells after filtering; the
row filter drops 44% of the 2,824 copies), gzip level 4:

| overlay | CSV | bytes/cell | gzipped | ratio |
|---|---|---|---|---|
| TFBS (one JASPAR motif) | 4.2 MB | 4.0 | 0.02 MB | **200x** |
| cCREs | 4.2 MB | 4.0 | 0.03 MB | **161x** |
| MAF coverage | 5.9 MB | 5.6 | 0.10 MB | 61x |
| MAF common_freq | 16.8 MB | 16.0 | 1.1 MB | 15.6x |
| MAF ref_freq | 17.0 MB | 16.2 | 1.8 MB | 9.5x |
| alignment | 4.2 MB | 4.0 | 0.5 MB | 9.2x |
| MAF base_count | 30.1 MB | 28.7 | 3.9 MB | 7.7x |
| phyloP | 18.0 MB | 17.1 | 2.9 MB | 6.2x |
| gnomAD | 7.2 MB | 6.8 | 1.4 MB | 5.3x |
| CADD | 16.9 MB | 16.0 | 4.4 MB | 3.9x |
| **MAF raw** | **9,669 MB** | **5,667** | - | - |

Two regimes. Binary or sparse tracks (TFBS, cCREs, coverage) compress 60-200x
and are effectively free. Continuous tracks (phyloP, CADD, gnomAD) compress only
4-6x, because real-valued data carries real entropy. No storage choice will
make them small, only smaller.

Projected to MIRb (~125,000 rows x 250 columns = 31 M cells, 30x MER41B; the
per-cell costs should hold, the row count assumes the same 56% survival under the
0.5 row threshold and has not been measured at this scale):

| overlay | CSV | gzipped |
|---|---|---|
| phyloP | 536 MB | 87 MB |
| CADD | 502 MB | 129 MB |
| MAF base_count | 900 MB | 117 MB |
| gnomAD | 213 MB | 40 MB |
| cCREs, or each TFBS motif | 125 MB | ~1 MB |
| **MAF raw** | **178 GB** | - |

A complete overlay set for one large family is roughly **2.5 GB as CSV, ~370 MB
gzipped**, with one exception.

### Do not run `MAF raw` on large families

At 5,667 bytes per cell it is three orders of magnitude heavier than any other
track. It is the nested-object path in `map_data()`: `canvas = np.full(shape,
None)`, every cell holding a per-species list, serialised as Python list reprs.
At MIRb scale that is ~178 GB for a single output. The aggregated MAF variants
(`base_count`, `coverage`, `common_freq`, `ref_freq`) carry the same information
at 6-29 bytes/cell. `bulk_produce_mapper_overlays.smk` already requests
`mem_mb = 65536` for `raw` against 8192 for every other track.

### Storage format

Numeric matrices are written with `pd.DataFrame(...).to_csv(...)`. Two cheap
improvements, in increasing order of disruption:

1. **Compress on write.** `to_csv("....csv.gz")` is transparent in pandas and
   takes a full MIRb overlay set from ~2.5 GB to ~370 MB. One-line change.
2. **Stop using text for numeric matrices.** phyloP costs 17.1 bytes/cell as
   decimal text against 4 as float32, a 4x penalty on disk, plus a pandas parse
   on every read. `np.save` or an h5py dataset would cut phyloP from 536 MB to
   125 MB before compression and load instantly. Binary tracks (TFBS, cCREs) are
   0/1 and belong in `uint8`, at 1 byte/cell.

For the parsed alignment cache, `parse_alignment()` uses `compression="lzf"`.
On the full MIRb matrix (4.98 GB uint8, 99.32% gaps) `lzf` is strictly dominated:

| method | size | ratio | time |
|---|---|---|---|
| h5py lzf (current default) | 150 MB | 33x | 13 s |
| h5py gzip-1 | 120 MB | 41x | 11 s |
| **h5py gzip-4** | **85 MB** | **58x** | 21 s |
| h5py gzip-9 | 56 MB | 89x | 208 s |
| scipy CSR sparse | 170 MB | 29x | 15 s |

`gzip-1` beats `lzf` on both size and time; `gzip-4` is the best balance. Sparse
storage is *worse* than generic compression here. int32 indices cost 4 bytes
per non-gap cell (136 MB) against 34 MB of actual data.

Note that compression is a disk-side concern only. `map_data()` needs the array
live and uncompressed, so none of this affects the peak memory in
[Matrix footprint](#matrix-footprint).

## Benchmark results

MIRb, 223,029 copies, 268 bp consensus, 29.3% mean divergence, median copy length
153 bp. 16 threads on `icelake-himem`.

Aligner RSS is the memory the aligner itself needs; overlay peak is what
`map_data()` then needs to use the result (see [Matrix footprint](#matrix-footprint)).

| method | wall | aligner RSS | width | overlay peak | bases lost |
|---|---|---|---|---|---|
| `hmmalign --outformat afa` | **2:25** | **3.13 GB** | 4,596 | 17.4 GB | **0.00%** |
| `famsa` (default) | 6:41 | 13.4 GB | **3,750** | **14.2 GB** | **0.00%** |
| `famsa -medoidtree` | 17:08 | 1.21 GB | 13,201 | 50.1 GB | 0.00% |
| `mafft --retree 1 --parttree` | 36:26 | 14.7 GB | 25,365 | 96.2 GB | 0.00% |
| `mafft --dpparttree --partsize 50` | 1:10:38 | 12.4 GB | 22,313 | 84.6 GB | 0.00% |
| `hmmalign --trim` (match only) | 2:19 | 2.98 GB | 334 | 1.3 GB | 10.77% |
| `famsa -medoidtree -remove-rare-columns` | 8:38 | 1.24 GB | 255 | 1.0 GB | 2.10% |
| `kalign --mode fast` | 2:10 | 4.76 GB | 7,634 | 28.9 GB | 0.00% |
| `mafft --treeout --reorder` (current default) | n/a | ~150 GB est. | n/a | n/a | n/a |

Alignment accuracy, scored as sum-of-pairs against an L-INS-i reference on 500
copies:

| method | SP |
|---|---|
| *ceiling: L-INS-i vs itself, shuffled input* | *0.9177* |
| `hmmalign`, L-INS-i seed | 0.5132 |
| `hmmalign`, FAMSA seed | 0.4519 |
| *G-INS-i (reference-class method)* | *0.4366* |
| `mafft --dpparttree --partsize 50` | 0.2809 |
| `mafft --treeout --reorder` | 0.2266 |
| `mafft --retree 1 --parttree` | 0.1993 |
| `famsa` | 0.1710 |
| `kalign --mode fast` | 0.0323 |

Occupancy structure is stable across scale. For FAMSA on MIRb:

| threshold | 5,000 copies | 223,029 copies |
|---|---|---|
| >=1% | 369 | 378 |
| >=5% | 260 | 258 |
| >=10% | 246 | 244 |
| >=50% | 157 | 158 |

The homologous core is ~250 columns whether 5,000 or 223,029 copies are aligned.
Width grew 2.7x across a 45x increase in sequences; useful width did not grow at
all.

## How far you can trust the accuracy numbers

The SP column has a low ceiling and should not be over-read.

- L-INS-i and G-INS-i, both reference-class aligners, agree with each other only
  **SP = 0.4366** on this family. At ~29% divergence with 153 bp copies there is
  no single defensible MSA of MIRb.
- L-INS-i against itself on shuffled input scores **0.9177**, so ~8% of pairs
  move under nothing but input order.
- The reference *is* L-INS-i, so L-INS-i-seeded profiles share its biases and
  gain some unearned credit. The FAMSA-seeded figure (0.4519) is the bias-free
  lower bound, and it already sits at the inter-reference level.
- FAMSA's low SP (0.1710) reflects **different gap placement, not worse
  alignment**. Reference-free mean pairwise identity within aligned columns puts
  it at 0.459 against L-INS-i's 0.466. SP cannot discriminate between them at
  this resolution; it can and does discriminate Kalign (0.276), which is a real
  quality failure.

This low ceiling is itself the strongest argument for the profile approach: it
replaces an ill-posed all-vs-all question with a well-posed one: where does
each copy sit relative to a shared model.

## Provenance

Benchmark harness, inputs, alignments and logs:
`/rds/project/rds-XrHDlpCeVDg/users/jeremi/02-ABCDom/01-data/00-ma_mapper/aligner_benchmark/`

- `methods.tsv`, `bench.sh`: scaling array (10 methods x 5,000 / 20,000 / 223,029)
- `accuracy.sh`, `sp_score.py`: SP/TC scoring against the L-INS-i reference
- `control.sh`: reference-stability controls (shuffled-input L-INS-i, G-INS-i)
- `seedtest.sh`: profile seed quality, on held-out copies
- `hmm_full.sh`, `hmm_afa_full.sh`: full-scale `hmmalign` variants
- `input/MIRb_linsi.hmm`, `input/MIRb_seed1000_linsi.afa`: reusable MIRb profile
- `out/*.result`: one row per run with method, size, n, exit, wall, max RSS, width, n_out

Versions: MAFFT v7.525, FAMSA 2.4.1, HMMER 3.4, Kalign 3.6.0. Consensus lengths
were derived as the mode of `repEnd + repLeft` over all copies in
`hg38.fa.out.tsv`; mode agreement is >=0.78 for every family quoted.

Run 2026-08-17.
