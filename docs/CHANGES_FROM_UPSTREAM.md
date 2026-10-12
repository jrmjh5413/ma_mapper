# Changes from the original `ma_mapper`

What this fork changes, and what it does not, relative to the original package.
Read it alongside [`../CHANGELOG.md`](../CHANGELOG.md), which is organised by
release; this document is organised by **what came from whom**, and is the
basis of the contribution statement in any paper describing the tool.

## Provenance

| Layer | Origin |
|---|---|
| The multiple-alignment overlay approach | Imbeault lab (Ecco *et al.*, 2016) |
| The `ma_mapper` package that implements it | **Pakkanan Chansongkrow**, 2023–2025; described in their PhD thesis (University of Cambridge, 2024) |
| Everything in this document | Jeremi Jaksina |

The original lives at **`github.com/pkncsk/ma_mapper`**. This fork is
**`github.com/jrmjh5413/ma_mapper`**.

The engine is not claimed here. What is claimed is the correctness work, the
packaging and the tests below.

## Where the fork diverges

Both branches share history up to **`8f7d870`** (2025-07-21, "remove history
folder"), which is still the tip of the original's `main`. Nothing from
upstream is missing from this fork. Three commits sit on top of it:

| commit | date | what it is |
|---|---|---|
| `761c3b1` | 2026-06-04 | MAF extraction: species-level deduplication and three smaller fixes |
| `38e603b` | 2026-10-11 | the VCF strand fix, the mapper gap/dtype fixes, the large-family alignment guide |
| `5c699e0` | 2026-10-11 | the audit: the remaining fixes, a command line, 105 tests, CI, packaging |

In total, against `8f7d870`: **43 files, +5,508 / −1,741 lines**. Per module:

| module | upstream lines | changed here |
|---|---|---|
| `mapper.py` | 489 | +570 / −152 |
| `sequence_alignment.py` | 186 | +328 / −186 |
| `extract_bam.py` | 152 | +276 / −135 |
| `extract_bed.py` | 143 | +220 / −129 |
| `utility.py` | 74 | +196 / −66 |
| `extract_maf.py` | 522 | +174 / −77 |
| `extract_vcf.py` | 83 | +167 / −67 |
| `extract_bigwig.py` | 71 | +117 / −68 |
| `__init__.py` | 17 | +84 / −16 |
| `plots.py` | 722 | +38 / −32 |
| `custom_cmap.py`, `gzmaf.py` | 524, 194 | untouched |
| `_io_common.py`, `cli.py` | — | new, 181 and 378 lines |

---

## A. Changes that alter results

Three. Anything computed with the original version should be checked against
this list.

### A1 · Minus-strand VCF extraction (`extract_vcf.py`, `38e603b`)

The original reversed the window **once per variant** instead of once per
window, because the reversal sat inside the per-variant loop. With *n* variants
in a window, every other value counting back from the last was mirrored about
the window centre, and with an even *n* the window was returned in the wrong
orientation entirely. A window holding one variant was correct, which is why
sparse tracks looked right and dense ones (gnomAD) did not. Minus-strand copies
are about half of any family.

Measured, old against new, on a window of ten bases with values 0.1, 0.2 and
0.4 at offsets 0, 2 and 6:

| | expected | original |
|---|---|---|
| plus strand | `0.1 · 0.2 · · · 0.4 · · ·` | same |
| minus strand, 3 variants | `· · · 0.4 · · · 0.2 · 0.1` | `· · 0.2 0.4 · · · · · 0.1` |
| minus strand, 2 variants | `· · · · · · · 0.2 · 0.1` | `0.1 · · · · · · 0.2 · ·` |
| minus strand, 1 variant | `· 0.1` | same |

Regression test: `tests/test_extract_vcf.py`, which checks the output against an
independently constructed array rather than against a recorded result. No
reproduction switch: the old behaviour was simply wrong.

### A2 · Read counts above 255 (`extract_bam.py`, `5c699e0`)

The read profile was `uint8`, which wraps silently: 300 reads at one base were
reported as 44. The smoothed profile was `float16` — three decimal digits,
saturating at 65504. Now `uint32` and `float32`. Regression test:
`tests/test_extractors.py::TestBam::test_read_counts_above_255_do_not_wrap`.

### A3 · The second column pass of `create_filter` (`mapper.py`, `5c699e0`)

The filter runs columns → rows → columns. The second column pass was meant to
recount occupancy among the copies that survived the row filter. Instead it

1. **added** those counts to the counts already accumulated from *all* copies,
   so the final mask came from a mixture of the two populations;
2. judged `col_content_threshold` against the original copy count, which makes
   the threshold stricter in proportion to how many copies were removed;
3. excluded each copy's own last non-gap column from the span denominator
   (`col_counts_max[nonzeros[0]:nonzeros[-1]]`), so a column could be non-gap
   in more copies than were counted as spanning it.

**The old masks remain reproducible exactly**: `create_filter(...,
legacy_second_pass=True)`, or `ma-mapper overlay --legacy-second-pass`. That
path was compared against the original code on five real families at three
threshold settings — every row and column mask identical.

Effect of the correction on the surviving column count (row counts are
unchanged except MamRep434 at the package default, 1796 → 1795):

| family | 0.05/0.05/0.5 | 0.10/0.10/0.10 | 0.50/0.10/0.50 |
|---|---|---|---|
| MER11A | 1233 → 1226 | 1215 → 1208 | 1029 → 1029 |
| MER11C | 1143 → 1141 | 1123 → 1111 | 1077 → 1077 |
| L1PA2 | 6053 → 6035 | 6035 → 6031 | 6034 → 6031 |
| LFSINE_Vert | 603 → 737 | 516 → 363 | 168 → 146 |
| MamRep434 | 789 → 682 | 673 → 523 | 169 → 171 |

Under 1 % on young compact families, 15–30 % on the two ancient ones.

---

## B. Crashes and silent wrong answers

Each was a case where the original either stopped with a misleading error or
returned something wrong without saying so. Grouped by module; every one has a
test.

### `mapper.py`

- **A one-shot iterator produced an all-gap matrix.** `parse_alignment` walked
  its input twice; given a generator the second pass saw nothing and every cell
  stayed zero, with no error.
- **Uninitialised memory in the gap test.** The matrix was allocated with
  `np.ndarray`, so any character outside `{-,A,C,T,G,N}` left whatever was in
  memory, making the `row > 0` non-gap test non-deterministic.
- **`.` was read as an unknown base**, not a gap. HMMER writes `.` in insert
  columns, so an `hmmalign` alignment that was not translated to `-` silently
  broke the per-base invariant.
- **Empty, malformed and ragged alignments** raised `UnboundLocalError` from a
  later line instead of saying what was wrong. Fifteen high-copy families were
  skipped in the 2026-08-08 matrix build for exactly this reason, because
  MAFFT had written zero-byte output.
- **No saved matrix could be read back**: `parse_alignment` wrote the HDF5
  dataset `parsed_alignment`, `import_parsed_alignment` looked for
  `result_array`.
- **Any family containing a fragmented copy crashed.** `parse_and_filter`
  applied the row mask to the coordinate table positionally; a copy split into
  blocks contributes one coordinate row per block, so the mask was the wrong
  length (`IndexError: Boolean index has wrong length`). Both FASTA conventions
  are now handled — one record per copy with the fragments in the header, and
  one record per fragment.
- **`map_and_overlay(coordinate_file=…)` paired the wrong copy with the wrong
  data.** The permutation was built as the position of each alignment row in
  the coordinate table and used directly as an index, which is the inverse of
  the permutation needed. It is correct for an identity or a single swap and
  wrong for anything else; a three-cycle demonstrates it.
- **A `coordinate_file` path raised `NameError`**: the existence check tested a
  name that is only assigned afterwards.
- **The documented `data_format` values were rejected by the dispatch.** The
  annotation advertised `bam_max` and friends; those calls logged "dataformat
  field unspecified" and then failed with `UnboundLocalError`. The accepted
  values are now `mapper.DATA_FORMATS` and an unknown one raises at once.
- **`base_count` crashed on an N outside the first column**, because
  `np.bincount(..., minlength=5)` sized the output from the first column alone.
- `map_data` now states which copy has how many values against how many
  non-gap bases, instead of raising a bare numpy broadcast error — or, when the
  counts coincide, quietly mis-assigning rows.
- The iterative filter could not be reached through `map_and_overlay` at all
  (`filter_mode` was never forwarded).

### `extract_maf.py`

- **Three arguments were delivered to the wrong parameters.** `maf_io` passed
  `target_species, count_arg, e_value_df, internal_id_df, species_list`
  positionally while `extract_maf` declared `..., species_list, e_value_df,
  internal_id_df`. A `species_list` therefore arrived as `internal_id_df` and
  hit an undefined name (`NameError` in the worker); an `e_value_table` arrived
  as `species_list` and became a regex of its own column names. The species
  filter has never worked through `maf_io`.
- **The target's gap-column removal was discarded** by every branch of the
  species filter, which re-derived from the unfiltered frame: a `species_list`
  or `e_value_table` call returned more values than the window has bases.
- An absent target species printed an `IndexError` and carried on, failing on
  the next line with `UnboundLocalError`.
- `common_freq_nogap` divided the commonest non-gap count by a total that still
  included the gaps it had just discarded, and crashed on an all-gap column.
- `get_maf_filepath`'s gzip branch was unreachable: `…chr1.maf.gz` does not end
  with `…chr1.maf`.

### `extract_bed.py`

- The overlapping feature's score was read with `int()`: 0.75 became 0 and `.`
  raised `ValueError`. Scores are floats; a missing score counts as 1.
- A DataFrame of features raised `TypeError`, from a leftover second conversion
  of the same argument.
- `generate_new_id=True` renamed the copies for grouping but not in the
  intersect query, so they could not be found again.
- The intersect output is now read by column position; pybedtools' BED12
  auto-names only line up for a feature file of exactly six columns.

### `extract_bigwig.py`, `extract_bam.py`

- `extract_bigwig` left its handle unbound for an unexpected argument type and
  never closed a file it opened. An unknown chromosome surfaced as a pyBigWig
  `RuntimeError`; it now names the mismatch and suggests the `chr` prefix.
- An unknown `bam_format` raised `UnboundLocalError` from the next line.

### `sequence_alignment.py`, `utility.py`, `plots.py`

- `extract_coord_from_repeatmasker_table` ignored the subfamily when an
  internal-id table was given: a list of *N* subfamilies returned the whole
  table *N* times.
- `parse_alignment_to_fasta` wrote gaps as `N` and raised `KeyError` on a
  genuine N — its mapping had 0 as N and no entry for 5.
- `mafft_align` filled `--threadtb` and `--threadit` with the value of
  `nthread`, producing `--threadb None` when they were set alone.
- `sequence_io(save_to_file=True)` built the output path by appending
  `.fasta` to the coordinate table's **parent directory**.
- `utility.divergence_table_prep` computed its filter and returned the
  unfiltered table.
- `plot_colorbar`, `plot_basecount` and `plot_logos` tested an
  `output_filepath` that was never assigned, so saving raised `NameError`; the
  other plotting functions wrote PNGs **into the installed package directory**,
  read-only in a normal install.
- Every `*_io` function compared `save_to_file == True`, which is `False` for a
  string, so asking to save to a path wrote nothing and reported nothing. The
  reported path is also now the path written: `compress_pickle` replaces the
  extension with the compression's own, so `vcf_out.p` was really
  `vcf_out.lzma`.
- `plots` appended a hard-coded path in another user's home directory to
  `sys.path` at import, and imported `tabnanny` and a non-existent Biopython
  module.

---

## C. Behaviour changed on purpose

Not bug fixes; decisions, listed so nothing is a surprise.

| what | was | now |
|---|---|---|
| Logging | importing the package attached a `FileHandler` writing `ma_mapper.log` into the current directory — it failed where that directory is read-only, and littered every directory it ran in | no handler, no file; WARNING and above still reach stderr through the standard library fallback. `configure_logging(level, logfile)` or `MA_MAPPER_LOG` opt in |
| Worker processes | 40, hardcoded in three extractors | the CPUs the process may use (`os.sched_getaffinity`, i.e. the Slurm allocation), overridable with `MA_MAPPER_WORKERS`, `MA_MAPPER_VCF_WORKERS`, `MA_MAPPER_MAF_WORKERS` |
| Genome access | `SeqIO.to_dict` loaded the whole FASTA into memory — gigabytes of hg38 to fetch a few hundred bases | read through the `.fai` index with pysam, falling back to the old route with a warning |
| `normalise` | always divided by the copies with a base in the column, `nan` values included | same default, plus `denominator='measured'` to divide by the copies with a finite value |
| `repeatmasker_prep` | a hardcoded human chromosome list, silently dropping everything in any other assembly | `main_chr` and `discard_class` parameters, human as the default |
| Overlay canvas / alignment matrix | `float64` / `float64` | `float32` / `uint8`, on the two largest arrays in the pipeline |
| Fragment order | taken as given | sorted into ascending genomic order, with a warning when that was necessary — the minus-strand handling reverses fragment order and is only correct for sorted input |
| `*_io` side effects | added a `meta_id` column to the caller's DataFrame | the input frame is left alone |
| VCF handles | one `cyvcf2.VCF` open per window, i.e. per copy | one per process, per file. The old pattern re-read the tabix index and the gnomAD header once per copy; CSD3 support reported the resulting load on 2026-08-02 |
| `map_data` type sniffing | every element of every copy, in pure Python | one value per copy |
| `os.listdir` in `maf_io` | once per copy | once per chromosome |
| `sequence_io` grouping | rescanned the whole table per copy (quadratic) | one pass |

---

## D. Added

- **A command line**, which the package did not have at all: `ma-mapper
  sequence | align | overlay | plot | index-maf | info`, with `--help` on every
  option. `overlay` also writes the row names and the filtered alignment,
  without which a matrix cannot be interpreted later.
- **105 tests** (`tests/`) on fixtures built at test time: a four-copy family
  with one minus-strand copy and one fragmented copy, and a dataset in every
  supported format whose value at position *p* is *p*, so a wrong orientation,
  fragment order or offset cannot pass unnoticed. Every item in sections A and
  B has a test.
- **Continuous integration** (`.github/workflows/ci.yml`): Linux on Python
  3.10–3.12 plus macOS, ruff, pytest, and a check that the command line runs.
- **Packaging**: PEP 621 `pyproject.toml` with a console entry point and
  per-format extras; `setup.py` removed, because the two disagreed about the
  package name. The wheel builds clean.
- **`extract_maf.build_maf_index`**, so the Biopython index is built once
  rather than by every worker process at the same time, plus a warning when it
  is missing.
- **`_io_common.py`**, the shared coordinate handling: five copies of the same
  reader, grouper and save logic — and therefore five copies of the same three
  defects — became one.
- **Documentation**: `CHANGELOG.md`, `CITATION.cff`, this file,
  [`aligning_large_families.md`](aligning_large_families.md) (measured memory
  scaling of MAFFT's defaults, the profile-alignment route above ~50,000
  copies, the second quadratic wall in row ordering), a rewritten `README.md`
  with a "things worth knowing before you trust a number" section, and new
  sections in [`package_docs.md`](package_docs.md) on the command line,
  parallelism, logging and semantics.
- **`examples/`**: `toy_end_to_end.py`, which builds a chromosome, a family and
  a track from nothing, runs the whole chain and checks its own output; and
  `verify_overlay.py`, which rebuilds an overlay from any BED, genome and
  bigWig and compares a sample of copies against a direct reading of the same
  file.

## E. Deliberately not changed

- **The overlay method itself.** Per-base extraction, placement into alignment
  columns, and the three filter thresholds are the original design.
- **The public API.** Every function keeps its name and signature; new
  behaviour sits behind new keyword arguments with the old defaults. The
  dependent analysis code in `02-ABCDom/02-code` needed no edits.
- **`custom_cmap.py` and `gzmaf.py`**, untouched.
- **The plotting functions' appearance and parameters.** Only the save paths
  and the import-time path manipulation changed.
- **Pakkanan's `MafIO.get_spliced` replacement**, including the IUPAC
  collapsing of duplicated species rows and the fragment-overlap trimming: the
  logic is theirs and is unchanged. Note that importing `extract_maf` still
  monkey-patches Biopython process-wide.
- **The analysis scripts under `script/`**, kept for provenance and excluded
  from the installed package.

## F. How to check any of this

```bash
pytest                                   # 105 tests, ~8 s
python examples/toy_end_to_end.py /tmp/x # self-checking, needs no data
ma-mapper info                           # versions, and which binaries are present
```

The claims in section A were checked against the original code rather than
against recorded output:

- **A1** against an independently constructed array
  (`tests/test_extract_vcf.py`).
- **A3** by running the original `create_filter` and the new one with
  `legacy_second_pass=True` over five real families at three settings and
  requiring every mask to be identical.
- The whole chain end to end on real data: 312 LTR10A copies (175 on the minus
  strand) from
  `/rds/project/rds-XrHDlpCeVDg/users/jeremi/02-ABCDom/01-data/00-ma_mapper/subfam_matrices_000/LTR10A/metadata.tsv`,
  against `hg38.fa` and `hg38.phyloP447way.bw`, with sampled copies compared
  base by base to a direct pyBigWig read — which is what `examples/verify_overlay.py`
  does for any dataset.
