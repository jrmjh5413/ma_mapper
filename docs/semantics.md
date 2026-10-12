# Defaults and semantics

The defaults and definitions that decide what a number means. None of these is
a bug. Each one has caught someone out, and each one is easier to find here
than in a parameter list.

This is the only copy. The README, the package documentation and the tutorial
link here rather than repeating it.

| thing | what it actually does |
|---|---|
| `extract_bam(offset=...)` | defaults to **5**, shifting every read 5 bp downstream |
| `row_threshold` | non-gap *density within a copy's own span*, which correlates with element age. Not a length filter |
| `normalise(denominator='alignment')` | divides by the copies with a base in the column, `nan` values included. `'measured'` divides by the copies with a finite value |
| `normalise(method='perc_coverage')` | counts non-zero cells, so a genuine 0.0 reads as absent |
| MAF `count_arg='coverage'` | counts the species present in the block, deletions included |
| `generate_new_id=True` | one copy per row, so a fragmented element is not reassembled |
| gap, not measured, and zero | three different things in the matrix, and only two of them are 0 |
| alignments above ~50,000 copies | MAFFT's defaults run out of memory, and so does row clustering |

## The BAM offset is 5, not 0

Every read's 5' end is shifted five bases downstream before it is counted. That
suits a Tn5-style correction and nothing else. Pass `offset=0` for raw 5' ends.
The default is inherited and has not been changed, because changing it would
silently move every existing BAM result by five bases.

## `row_threshold` is a density, not a length

A copy survives if its non-gap bases are dense *within its own span*. A short
intact fragment passes; a full-length, deletion-riddled copy fails. The
threshold therefore correlates with element age: at 0.5 an ancient SINE loses
most of its copies, and the ones it keeps are the unusually intact ones, which
is a selection on exactly the covariate most TE analyses condition on.

For anything stratified by age use 0.10, or no row filter, and say which.
Measured on real families: at 0.5, LFSINE_Vert keeps 235 of 2,301 copies and
MamRep434 keeps 130 of 2,102, against 2,297 and 2,095 at 0.10.

## `normalise` divides by the copies with a base, not the copies with a value

With `denominator='alignment'`, the default and the only behaviour before
0.9.0, a column's mean is the sum of the data over the number of copies that
have a base there, **including copies whose value is `nan`**. A track that does
not cover part of the family therefore reads as a low value rather than as
missing. `denominator='measured'` divides by the copies that have a finite
value.

Which one is right depends on whether a gap in the data means zero or means
not measured. For a bigWig it means not measured.

## `perc_coverage` treats a real zero as absent

It counts non-zero cells, so it is a coverage measure only for data where 0
means "nothing here", such as peak overlap. For signal that can genuinely be
zero, it undercounts.

## MAF `coverage` counts species, including deletions

It answers "how many species align here", not "how many have a base here". A
species whose row exists in the block but whose base at that column is a
deletion is counted. If the number is going to be the denominator of a
conservation fraction, that distinction matters.

## Gap, not measured, and zero are three different things

| in the matrix | means |
|---|---|
| `alignment == 0`, overlay cell 0 | the copy has no base in this column |
| `nan` | the copy has a base, the track does not cover it (bigWig) |
| 0.0 | measured, and zero |

BED input has no third case: absence of a feature is a real zero. This is why
`overlay` writes the alignment alongside the matrix, and why collapsing a
column needs the `denominator` decision above.

## Large families need a different aligner

Above roughly 50,000 copies MAFFT's default settings run out of memory, and so
does the row clustering in any downstream plot. See
[Aligning large TE families](aligning_large_families.md) for the profile
alignment route and the measured scaling.

Never use `hmmalign --trim`. It deletes bases, which breaks the invariant the
whole package rests on: a copy's extracted values must number exactly as many
as its non-gap bases. Always translate HMMER's `.` to `-`, or rely on 0.9.0,
which reads `.` as a gap.

## The alignment and the data must describe the same copies

In the same order. `map_data` checks that each copy's value count equals its
non-gap base count and names the copy that disagrees. The error is worth taking
seriously: it usually means an alignment and a coordinate table built from
different copy sets, which no amount of filtering will fix.

`examples/verify_overlay.py` re-reads any bigWig independently and compares a
sample of copies base by base. It is the check worth running on a new dataset
before building anything on it.
