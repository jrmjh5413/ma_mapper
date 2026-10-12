# New in 0.9.0

Parameters and functions that the rest of this reference does not cover. Pages
01–09 document 0.8.x and remain accurate for everything they describe.

## `mapper.create_filter`

`ma_mapper.mapper.create_filter(parsed_alignment, col_threshold=0.50, col_content_threshold=0.10, row_threshold=0.50, save_to_file=False, legacy_second_pass=False)`

**legacy_second_pass : bool, default: False**
Reproduce the pre-0.9.0 arithmetic of the second column pass, which added the
surviving copies' occupancy to the counts already taken from all copies rather
than recounting, judged `col_content_threshold` against the original copy
count, and left each copy's last non-gap column out of the span denominator.
Set it to rebuild a matrix built before 0.9.0 bit-for-bit; leave it alone for
new work. The difference is under 1 % of columns on young compact families and
15–30 % on ancient ones (`../CHANGES_FROM_UPSTREAM.md` section A3).

## `mapper.map_and_overlay`

`ma_mapper.mapper.map_and_overlay(alignment, data_file, data_format, filter=True, coordinate_file=None, extension_length=None, source_fasta=None, filter_mode='one_pass', return_details=False, **kwargs)`

**data_format : str**
One of `mapper.DATA_FORMATS`: `bigwig`, `bed`, `maf`, `vcf`, or a BAM readout —
`read_sum`, `read_max`, `read_min`, `read_forward`, `read_reverse`,
`normal_max`, `normal_min`, `normal_forward`, `normal_reverse`. An unknown
value raises `ValueError` immediately. (Before 0.9.0 the type annotation
advertised `bam_max` and similar, which the dispatch did not accept.)

**filter_mode : {'one_pass', 'iterative'}, default: 'one_pass'**
Which filter to apply when `filter` is true. The iterative filter could not be
reached through this function before 0.9.0.

**return_details : bool, default: False**
Return `(matrix, coordinate_table, alignment_matrix)` instead of just the
matrix: the coordinates of the surviving copies in matrix row order, and the
filtered alignment. Without them a matrix has no row labels and no way to tell
a gap from a measured zero.

## `mapper.normalise`

`ma_mapper.mapper.normalise(alignment_matrix, data_matrix, method='average', denominator='alignment')`

**denominator : {'alignment', 'measured'}, default: 'alignment'**
What `average` and `perc_coverage` divide by. `alignment` — the only behaviour
before 0.9.0 — counts the copies with a base in the column, *including those
whose value is `nan`*, so a track that does not cover part of the family reads
as a low value rather than as missing. `measured` counts the copies that have a
finite value there. A column no copy occupies returns `nan` rather than `inf`.

## `mapper.base_count`

Returns six rows, indexed by the parse codes 0 gap, 1 A, 2 C, 3 T, 4 G, 5 N.
Before 0.9.0 it sized its output from the first column alone and failed on an
alignment whose later columns held an N.

## `utility.repeatmasker_prep`

`ma_mapper.utility.repeatmasker_prep(repeatmasker_filepath, main_chr=None, discard_class=None)`

**main_chr : list or None, default: None**
Chromosomes to keep; `None` means `utility.HUMAN_MAIN_CHROMOSOMES`. **Pass this
for any assembly other than human** — the human list keeps mouse `chr1..chr19`
by coincidence and silently drops everything else. `[]` keeps every chromosome.

**discard_class : list or None, default: None**
Repeat classes to drop; `None` means `utility.DEFAULT_DISCARD_CLASS` (simple
repeats, satellites, low complexity, structural RNAs and the question-marked
classes). `[]` keeps every class.

## `utility.divergence_table_prep`

**drop_unaligned : bool, default: True**
Drop species marked `3. Zoonomia (not in alignment)`. Before 0.9.0 this filter
was computed and then discarded, so the function returned the unfiltered table.

## `extract_maf.build_maf_index`

`ma_mapper.extract_maf.build_maf_index(maf_file, target_species='Homo_sapiens', chrom=None)`

Build Biopython's sqlite index for a MAF once, before any parallel extraction.
`extract_maf` constructs a `MafIndex` per copy; if the index does not exist yet,
every worker process starts building the same one into the same file at the same
time. Returns the index path. `maf_io` warns when the index is missing.

## `extract_maf.maf_io`

**species_list : list or None**
Keep only these species. It must include `target_species`, or the reference row
is dropped and the call raises. Before 0.9.0 this argument was delivered to the
wrong parameter of the worker and never took effect.

## `extract_bam.extract_bam` and `bam_io`

**offset : int, default: 5**
Bases to shift each read's 5' end by. **The default is not zero**: it suits a
Tn5-style correction and nothing else. Pass `offset=0` for raw 5' ends.

**drop_dup, drop_unmap, drop_2, drop_multimap : bool, default: False**
Skip reads flagged duplicate, unmapped, secondary or supplementary, or carrying
`NH` > 1. `drop_multimap` only works for aligners that write the `NH` tag.

`read_*` readouts return `uint32` counts and `normal_*` return `float32`.
Before 0.9.0 they were `uint8` and `float16`, so a count above 255 wrapped.

## Every `*_io` function

**save_to_file : bool or str, default: False**
A path is written to as given; `True` writes a default name beside the
coordinate table, or in the current directory when the table was passed as a
DataFrame. The function returns the output path when it saves, and the payload
otherwise. Before 0.9.0 the test was `save_to_file == True`, which is `False`
for a string, so passing a path wrote nothing and reported nothing.

The default names end in `.lzma`, which is what the files actually are:
`compress_pickle` replaces the extension with the compression's own, so the old
`*_out.p` names were never the names on disk.

**generate_new_id : bool, default: False**
Number the copies `entry_0, entry_1, …` instead of grouping rows by column 3.
Each row becomes its own copy, so a fragmented element is **not** reassembled.
It cannot be combined with `coordinate_file` in `map_and_overlay`, where rows
are matched to the alignment by name.

## `ma_mapper.configure_logging`

`ma_mapper.configure_logging(level=logging.INFO, logfile=None, stream=True)`

Importing the package no longer configures logging or creates a file; it used to
attach a handler writing `ma_mapper.log` into the current working directory.
Messages of WARNING and above still reach stderr through the standard library's
fallback. Call this for progress output, or set `MA_MAPPER_LOG` before import.

## `ma_mapper.provenance`

`write_run_record(path, inputs, parameters, outputs=None, command=None, extra=None)`

Write a JSON record of how a matrix was made: the inputs with their sizes and
checksums (files above `MAX_HASH_BYTES`, 256 MB, are identified by size and
mtime instead), the parameters, the outputs, the package versions, the host and
any Slurm job variables. `ma-mapper overlay` writes one automatically;
`read_run_record(path)` reads it back.

## Environment variables

| variable | effect |
|---|---|
| `MA_MAPPER_LOG` | a path to append log records to, from import time |
| `MA_MAPPER_WORKERS` | worker processes for the parallel extractors |
| `MA_MAPPER_VCF_WORKERS` | the same, for `extract_vcf` only |
| `MA_MAPPER_MAF_WORKERS` | the same, for `extract_maf` only |

Without them the extractors use the CPUs the process is allowed
(`os.sched_getaffinity`, i.e. the Slurm allocation). Before 0.9.0 they
hardcoded 40.

## The command line

`ma-mapper sequence | align | overlay | plot | index-maf | info`, documented in
[`../package_docs.md`](../package_docs.md#command-line) with `--help` on every
option.
