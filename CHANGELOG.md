# Changelog

## 0.9.0 (2026-10-11)

First release with a test suite, a command line and a package build. The
headline is a set of correctness fixes: **three of them change numbers that
earlier analyses depended on**, and those three are listed first.

### Corrected results

Anything built with an earlier version should be checked against this list.

- **VCF, minus strand.** `extract_vcf` reversed the window once per variant
  instead of once per window. With *n* variants in a window, every other value
  counting back from the last was mirrored about the window centre, and with an
  even *n* the window was returned in the wrong orientation entirely. A window
  with a single variant was correct, which is why sparse tracks looked fine and
  dense ones (gnomAD) did not. Minus-strand copies are about half of any
  family. Verified against an independent construction in
  `tests/test_extract_vcf.py`.
- **BAM, deep pile-ups.** The read profile was `uint8`, which wraps silently at
  256: 300 reads at one base were reported as 44. The smoothed profile was
  `float16` (three decimal digits, saturating at 65504). Now `uint32` and
  `float32`.
- **The second column pass of `create_filter`.** The pass that was meant to
  recount column occupancy among the surviving copies instead *added* those
  counts to the counts already accumulated from all copies, and judged the
  content threshold against the original copy count. It also excluded each
  copy's last non-gap column from the span denominator, so a column could be
  non-gap in more copies than were counted as spanning it. Pass
  `legacy_second_pass=True` to reproduce the old masks exactly; that was
  checked bit-for-bit on five real families (MER11A, MER11C, LFSINE_Vert,
  MamRep434, L1PA2) at three threshold settings. The corrected filter changes
  the number of surviving columns by under 1 % on young compact families
  (MER11A at 0.05/0.05/0.5: 1233 to 1226) and by more on ancient ones
  (MamRep434 at the same setting: 789 to 682).

### Fixed: crashes and silently wrong behaviour

- `parse_alignment` walked its input twice; given a one-shot iterator the
  second pass saw nothing and the matrix came back **entirely gaps**, with no
  error. It also allocated with `np.ndarray`, leaving the cells uninitialised,
  so any character outside `{-,A,C,T,G,N}` produced non-deterministic values in
  the `row > 0` gap test. Now allocated with `np.zeros`, `uint8`, and `'.'`
  (HMMER's insert-column gap) is read as a gap rather than an unknown base.
  Unmapped characters are counted and reported.
- `parse_alignment` raises on an empty or malformed alignment instead of
  failing later with `UnboundLocalError`, and on a ragged one it names the
  record and both lengths.
- `import_parsed_alignment` looked for an HDF5 dataset called `result_array`
  while `parse_alignment` wrote `parsed_alignment`, so no saved matrix could
  ever be read back.
- `parse_and_filter` applied the row mask to the coordinate table positionally.
  A copy split into several fragments has one coordinate row per fragment, so
  any family containing one died with `IndexError: Boolean index has wrong
  length`. Both FASTA conventions are now handled: one record per copy (names
  unique, fragments in the header) and one record per fragment.
- `map_and_overlay(coordinate_file=...)` built the row permutation backwards —
  the position of each alignment row in the coordinate table, used directly as
  an index, which is the inverse permutation. It is correct for an identity or
  a single swap and pairs the wrong copy with the wrong data for anything else.
- `map_and_overlay` tested `os.path.isfile(coordinate_table)` before that name
  existed, so passing a `coordinate_file` path raised `NameError`.
- `map_and_overlay`'s type annotation advertised formats (`bam_max`, ...) that
  its dispatch did not accept; those calls logged "dataformat field
  unspecified" and then failed with `UnboundLocalError`. The accepted values
  are in `mapper.DATA_FORMATS` and an unknown one raises immediately.
- `map_data` checks that each copy's value count matches its non-gap base count
  and says which copy disagrees, instead of raising a bare numpy broadcast
  error — or, when the counts happen to match, silently mis-assigning rows.
- `maf_io` passed five fixed arguments positionally in an order that did not
  match `extract_maf`'s signature: `species_list`, `e_value_df` and
  `internal_id_df` were rotated. A `species_list` arrived as `internal_id_df`
  and hit an undefined name (`NameError` inside the worker); an `e_value_table`
  arrived as `species_list` and became a regex of its own column names. The
  worker is now called with keywords.
- `extract_maf` dropped the target species' gap columns and then went back to
  the *unfiltered* frame in every branch of the species filter, so a
  `species_list` or `e_value_table` call returned more values than the window
  has bases.
- `extract_maf` printed an `IndexError` and carried on when the target species
  was absent, failing on the next line with `UnboundLocalError`; it now says
  which copy, which file and that `species_list` has to include the target.
- `common_freq_nogap` divided the commonest non-gap count by a total that still
  included the gaps it had just discarded, and crashed on an all-gap column.
- `get_maf_filepath`'s gzip branch was unreachable: `...chr1.maf.gz` does not
  end with `...chr1.maf`.
- `base_count` sized its output from the first column alone
  (`np.bincount(..., minlength=5)`), so an alignment whose first column held no
  N but whose later columns did failed with a broadcast error. It now returns
  six rows (gap, A, C, T, G, N) and is vectorised.
- `extract_bed` read the overlapping feature's score with `int()`: a score of
  0.75 became 0 and `.` raised `ValueError`. Scores are floats and a missing
  score counts as 1. The intersect output is read by column position rather
  than by pybedtools' BED12 auto-names.
- `extract_bed.bed_io` crashed on a DataFrame of features, from a leftover
  second conversion of the same argument, and ignored `generate_new_id` when
  building the query, so the renamed copies could not be found again.
- `sequence_alignment.extract_coord_from_repeatmasker_table` ignored the
  subfamily when an internal-id table was given: a list of N subfamilies
  returned the whole table N times.
- `parse_alignment_to_fasta` wrote gaps as `N` and raised `KeyError` on a
  genuine N: its mapping had 0 as N and no entry for 5.
- `mafft_align` filled `--threadtb` and `--threadit` with the value of
  `nthread`, producing `--threadb None` if they were set alone. It no longer
  runs through a shell.
- `utility.divergence_table_prep` computed its filter and returned the
  unfiltered table.
- `plot_colorbar`, `plot_basecount` and `plot_logos` tested an
  `output_filepath` that was never assigned, so saving raised `NameError`; the
  other plotting functions defaulted to writing a PNG **into the installed
  package directory**, which is read-only in a normal install. Figures now
  default to the current directory, and a path passed as `save_to_file` is
  honoured.
- Every `*_io` function compared `save_to_file == True`, which is `False` for a
  string: asking to save to a path wrote nothing and reported nothing. The
  reported path is also the path actually written — `compress_pickle` replaces
  the extension with the compression's own, so `vcf_out.p` was really
  `vcf_out.lzma`. Default names now end in `.lzma`.
- A missing coordinate table was logged and then crashed with `NameError`
  further down; it raises `FileNotFoundError` with the path.
- `extract_bigwig` left `bigwig` unbound for an unexpected argument type, and
  never closed a file it opened. An unknown chromosome name now says so and
  suggests the `chr`-prefix mismatch, rather than surfacing as a pyBigWig
  `RuntimeError`.
- `plots` appended a hard-coded path in another user's home directory to
  `sys.path` at import, and imported a module (`tabnanny`) and a dead
  Biopython path that do not belong to this package.

### Changed

- **Logging.** Importing the package no longer configures logging or creates a
  file. It used to attach a `FileHandler` writing `ma_mapper.log` into the
  current working directory, which failed outright where that directory is
  read-only and littered every directory it ran in (six such files were
  committed to this repository). Messages of WARNING and above still reach
  stderr through the standard library's fallback. Call
  `ma_mapper.configure_logging(level, logfile)` for progress output, or set
  `MA_MAPPER_LOG` before import.
- **Worker counts.** The three parallel extractors each hardcoded 40 worker
  processes. They now default to the CPUs the process is actually allowed to
  use (`os.sched_getaffinity`, i.e. the Slurm allocation) and respect
  `MA_MAPPER_WORKERS`, `MA_MAPPER_VCF_WORKERS` and `MA_MAPPER_MAF_WORKERS`.
- `normalise` takes `denominator='alignment'|'measured'`. The default is the
  old behaviour, which divides by the number of copies with a base in the
  column **including those whose value is `nan`** — so a track that does not
  cover part of the family reads as a low value rather than as missing.
  `'measured'` divides by the copies that have a finite value. A column no copy
  occupies yields `nan` rather than `inf`.
- `map_and_overlay` takes `filter_mode` (the iterative filter was unreachable
  through it) and `return_details`, which also returns the surviving copies'
  coordinates and the filtered alignment — without them a matrix has no row
  labels.
- `utility.repeatmasker_prep` takes `main_chr` and `discard_class`. The default
  chromosome list is human, which silently dropped everything in any other
  assembly.
- `*_io` functions no longer add a `meta_id` column to the caller's DataFrame.
- Fragments of a copy are sorted into ascending genomic order before
  extraction, with a warning when that was necessary: the minus-strand handling
  reverses fragment order and is only correct for sorted input.
- Overlay canvases are `float32` rather than `float64`, and parsed alignments
  `uint8`: half and an eighth of the memory, on the two largest arrays in the
  pipeline.
- `extract_maf` no longer calls `os.listdir` once per copy, and
  `sequence_io` no longer rescans the whole coordinate table once per copy
  (quadratic in the number of copies).
- `sequence_io` reads the genome through its `.fai` index (pysam) instead of
  loading the whole FASTA into memory with `SeqIO.to_dict` — several gigabytes
  for hg38, to fetch a few hundred bases. Its `save_to_file=True` path built the
  output name by appending `.fasta` to the coordinate table's *parent
  directory*.
- `map_data`'s decision between a numeric and an object canvas inspects one
  value per copy instead of every value of every copy in pure Python.

### Added

- A command line: `ma-mapper sequence | align | overlay | plot | index-maf |
  info`, with `--help` on every option.
- `extract_maf.build_maf_index`, so the Biopython index is built once instead
  of by every worker process at the same time, and a warning when it is missing.
- A test suite (105 tests) built on synthetic fixtures: a four-copy family with
  one minus-strand copy, one fragmented copy, and a dataset in each supported
  format whose value at position *p* is *p*, so a wrong orientation, fragment
  order or offset cannot pass.
- `pyproject.toml` with PEP 621 metadata, a single build backend, optional
  extras per format, and a console entry point. `setup.py` is gone; the two
  files disagreed about the package name.
- `CITATION.cff`.
