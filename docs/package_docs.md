# Package Documentation

## Quick links
- [Input preparation](#input-preparation)
    - [TE coordinates](#te-coordinates)
    - [Transposable element multiple sequeunce alignment (TE MSA)](#transposable-element-multiple-sequeunce-alignment-te-msa)
        - [Large families](#large-families)
    - [Genome-wide datasets](#genome-wide-datasets)
        - [`BED` file](#bed-file)
        - [`BAM` file](#bam-file)
        - [`BIGWIG` file](#bigwig-file)
        - [`VCF` file](#vcf-file)
        - [`MAF` file](#maf-file)
- [TE MSA conversion](#te-msa-conversion)
- [Data extraction and mapping](#data-extraction-and-mapping)
    - [`BED` file](#bed-file-1)
    - [`BAM` file](#bam-file-1)
    - [`BIGWIG` file](#bigwig-file-1)
    - [`VCF` file](#vcf-file-1)
    - [`MAF` file](#maf-file-1)
- [Visualization](#visualization)
- [Command line](#command-line)
- [Parallelism](#parallelism)
- [Logging](#logging)
- [Semantics worth checking before publishing a number](#semantics-worth-checking-before-publishing-a-number)
- [Reproducing results from before 0.9.0](#reproducing-results-from-before-090)

See also [Changes from the original `ma_mapper`](CHANGES_FROM_UPSTREAM.md).

## Input preparation
TE coordinates, TE MSA, and genome-wide dataset are required for the pipeline. Although there are different methods to prepare these files, the format requirement for each input file, and the preparation method and used including parameters and options in each showcase are outlined below:

### TE coordinates
The requirements for this input is to be in BED-like format (`chromosome`, `start`, `end`, and `strand` must be present) and optionally provide unique ID for TE entries. It is possible for one ID to have multiple entries in the table as the package supports fragment joining. Critically, the ID from TE coordinate should be the same as TE sequence ID in TE MSA. 

For the testing sample, the coordinates for THE1C TEs were obtained from [RepeatMasker output](https://www.repeatmasker.org/genomes/hg38/RepeatMasker-rm405-db20140131/hg38.fa.out.gz). The header rows were edited from two rows into one rows to facilitate importing by `pandas`. After importing, the THE1C entries were extracted from the table using `ma_mapper.sequence_alignment`. This module is a wrapper for `biopython.SeqIO` and `MAFFT` to streamline the sequence alignment process.

```python
	from ma_mapper import sequence_alignment
	
	coord_table=sequence_alignment.extract_coord_from_repeatmasker_table(
	subfamily = 'THE1C',
	repeatmasker_table = '/path/to/repeatmasker/output/file/',
	save_to_file = False
	)
```
The resulting `coord_table` is a `pandas.DataFrame` object that is in BED format and compatible with the pipeline. To save the table to a file, users can provide a destination path via `save_to_file`.

### Transposable element multiple sequeunce alignment (TE MSA)
The TEA MSA file is required to be in FASTA format with a header specifically generated using their coordinate:     `>ID::chromosome:start-end(strand)` to ensure proper row syncing across matrices in the downstream analyses. For this package, the `ma_mapper.sequence_io()` was used to extract sequence from the [human reference genome (`hg38`) obtained from UCSC](https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.fa.gz).

```python
	from ma_mapper import sequence_alignment
	
	te_seqeunces=sequence_alignment.sequence_io(
	coordinate_table=coord_table,
	source_fasta='/path/to/hg38/fasta',
	save_to_file='/target/path/for/TE/sequences')
```

`coordinate_table` can accept both a dataframe (as produced in the previous step) and a path to BED file. Although the function returns a `SeqIO` object, it is recommended to save the result as a FASTA file for use in the next step.

Then, the multiple sequence alignment was performed on TE sequences using `ma_mapper.sequence_alignment.mafft_align`.

```python
	from ma_mapper import sequence_alignment
	
	sequence_alignment.mafft_align(
	input_filepath='/path/to/TE/sequences',
	nthread=6,
	output_filepath='/target/path/for/TE/alignment'
	)
```

The output was saved in FASTA format. Users can adjust multithread parameters using `nthread, nthreadb, nthreadit` and provide additional flags using `mafft_arg` argument.

#### Large families

`mafft_align()` uses MAFFT's default algorithm (FFT-NS-2), whose memory grows as
about `n^1.19`. It is reliable up to a few tens of thousands of copies and fails
above roughly 65,000 — for the largest human families it would need ~150 GB.

For families above ~50,000 copies, use a profile alignment instead:

```bash
mafft --thread 16 --maxiterate 1000 --localpair seed1000.fasta > seed.afa
hmmbuild --dna --cpu 16 --informat afa family.hmm seed.afa
hmmalign --dna --outformat afa family.hmm all_copies.fasta | tr '.' '-' > family.afa
```

This aligns all 223,029 copies of MIRb in 2.5 minutes on 3.1 GB. Two rules:
never use `hmmalign --trim` or a match-column-only projection (they delete ~11%
of bases and silently break `map_data()`), and always translate `.` to `-`
(HMMER's insert-column gap character is not in `parse_alignment()`'s hash table).

Note that the aligner is only half the problem: the overlay pipeline's row
clustering (`pdist` + `linkage`) is also O(n^2) and needs 199 GB at 223,029
copies, so it fails at the same ~65,000-copy threshold. Use an aligner guide
tree (`famsa -gt nj -medoidtree -gt_export`, ~10 s) above that size.

Overlay outputs are written filtered (~250 columns), so they stay small — but
`MAF raw` costs ~5,700 bytes per cell against 4-29 for every other track, which
projects to ~178 GB for a single large family. Use the aggregated MAF variants
(`base_count`, `coverage`, `common_freq`, `ref_freq`) instead at that scale.

Benchmarks for ten aligner configurations, the constraints any alignment must
satisfy for `map_data()`, memory footprints per aligner, the row-ordering
comparison, and overlay output sizes and compressibility are in
[Aligning large TE families](/docs/aligning_large_families.md).

### Genome-wide datasets
To show case the capability of the framework, genome-wide datasets in this thesis were obtained from various publicly available sources as follows:

#### `BED` file
An example of `BED` file was obtained from [HOMER transcription factor motif prediction](http://homer.ucsd.edu/homer/data/motifs/homer.KnownMotifs.hg38.191020.bed.gz). The table was then separated by motif names into smaller table for efficient data retrieval.

#### `BAM` file 
An example `BAM` file was generated from ZNF267 ChIP-exo derived from [KZFP overexpression experiments](https://www.ncbi.nlm.nih.gov/sra?term=SRX2512744). The `FASTQ` file, retrieved from Sequence Read Archive (SRA) was mapped to [`GRCh38_noalt_as` reference genome](https://genome-idx.s3.amazonaws.com/bt/GRCh38\_noalt\_as.zip) using `bowtie2`. Subsequently, the mapped reads in `SAM` format was converted into `BAM`, then sorted and indexed by `samtools`. The read mapping operation was performed with default settings as shown below:

```bash
	# Step 1: Align with Bowtie2
	bowtie2 --threads 128 --local -x GRCh38_noalt_as -U /path/to/fastq/ -S /target/path/to/sam/
	# Step 2: Convert SAM to BAM
	samtools view -@ 128 -bS /target/path/to/sam/ > /target/path/to/bam/
	# Step 3: Sort BAM
	samtools sort -@ 128 -o /target/path/to/sam/sorted/ /target/path/to/sam/
	# Step 4: Index BAM
	samtools index /target/path/to/sam/sorted/
```
#### `BIGWIG` file
An example of `BigWig` file, a [phyloP track of Zoonomia multispecies alignment , was obtained from UCSC database](http://hgdownload.cse.ucsc.edu/goldenPath/hg38/phyloP447way) and was used without modification.

#### `VCF` file
Example `VCF` files were obtained from [the gnomAD database](https://gnomad.broadinstitute.org/data#v3) and stored in the same directory. After downloading, these files were indexed by `tabix` for efficient data retrieval.

The following snippet is the command for VCF file indexing:
```bash
	tabix -p vcf /path/to/zipped/vcf/file/
```

#### `MAF` file
\label{mafprep}
Example `MAF` files of [the updated Zoonomia dataset obtained from UCSC database](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/cactus447way/maf/). As `MafIO.MafIndex()` does not support a compressed `MAF` file, all files were unzipped before indexing for fast data retrieval. For this framework, `.mafindex` files are required to have exactly the same name as their `.maf` counterparts as shown below.
```bash
	chr1.maf
	chr1.mafindex
	...
	chrY.maf
	chrY.mafindex
```

The MAF files were indexed using the code snippet below:
```python
	# build mafindex file 
	from bio import MafIO
	maf_file = '/path/to/maf/'
	target_species = 'hg38'
	maf_id = f'{target_species}.{chrom}'
	mafindex_filedir = '.'.join(str.split(maf_file, sep='.')[:-1]) #remove .maf 
	mafindex_filepath = f'{mafindex_filedir}.mafindex'
	index_maf = MafIO.MafIndex(mafindex_filepath, maf_file, maf_id) 
```  


## TE MSA conversion
Alignment conversion is mandatory for data mapping as the numerical alignment matrix is needed as a reference for gap position. The input for this step is the aligned FASTA file with specific sequence header `>ID::chromosome:start-end(strand)`. 

```python
	from ma_mapper import mapper
	alignment_filepath = '/path/to/alignment/fasta/'
	
	filtered_alignment_matrix, alignment_coordinate  = mapper.parse_and_filter(alignment_file=alignment_filepath,col_threshold = 0.10, col_content_threshold = 0.10, row_threshold = 0.10)
``` 

The outputs of this step are a filtered numerical alignment matrix `numpy.ndarray` and a metadata table `pandas.DataFrame`. By default, these outputs are already filtered internally. For pre-filtered output, users can toggle the `preprocess_out` argument to make the function returns a pre-filtered numerical alignment matrix, a pre-filtered metadata, and filters.

Filtering behavior is customizable with several arguments:
- `filter`: enable/disable filtering
- `row_threshold`: set the row threshold
- `col_threshold`: set the column threshold
- `col_content_threshold`: sets the minimum non-gap content for a column

Additionally, to examine regions beyond TE boundaries, users can specify the extension length and the path to the genome sequence for sequence extraction using `extension` and `source_fasta` arguments respectively.      

## Data extraction and Mapping
The inputs for this step are the numerical alignment matrix, a path to genome-wide data file, and a type of file. Regardless of input type, this function always outputs a numerical data matrix that has the same gapped position as the numerical alignment matrix. The specific handling and parameters for each supported file type is described below.

### `BED` file
```python
	alignment_filepath = '/path/to/alignment/fasta'
	genomewide_data_filepath = '/path/to/bed/file/'
	
	filtered_alignment_matrix, alignment_coordinate  = mapper.parse_and_filter(alignment_file=alignment_filepath,col_threshold = 0.10, col_content_threshold = 0.10, row_threshold = 0.10)
	bed_data_matrix=mapper.map_and_overlay(alignment_filepath, genomewide_data_filepath,data_format='bed',strand_overlap=True, 
	col_threshold = 0.10, col_content_threshold = 0.10, row_threshold = 0.10)
``` 
Some genome-wide dataset `BED` files may lack strand information (i.e., `.` instead of `+,-`). The `strand_overlap` argument can be used to enable/disable strand-specific processing.

### `BAM` file

```python
	alignment_filepath = '/path/to/alignment/fasta'
	genomewide_data_filepath = '/path/to/bam/file/'
	
	filtered_alignment_matrix, alignment_coordinate  = mapper.parse_and_filter(alignment_file=alignment_filepath,col_threshold = 0.10, col_content_threshold = 0.10, row_threshold = 0.10)
	data_matrix_forward=mapper.map_and_overlay(alignment_filepath, genomewide_data_filepath,data_format='read_forward', 
	col_threshold = 0.10, col_content_threshold = 0.10, row_threshold = 0.10)
	data_matrix_reverse=mapper.map_and_overlay(alignment_filepath, genomewide_data_filepath,data_format='read_reverse', 
	col_threshold = 0.10, col_content_threshold = 0.10, row_threshold = 0.10)
```

For `BAM` files, users can choose to separate `read_forward` and `read_reverse` or the combined `read_sum`. Smoothing is also supported via `normal_forward` and `normal_reverse` with parameters for for probe length `probe_length` and smoothing length `smoothing_length`. 

Even though `BAM` files are recommended to be preprocessed before using this pipeline, basic filtering using `pysam` internal flags is supported. Users can enable arguments such as:
- `drop_duplicate`: remove duplicate reads
- `drop_unmapped`: remove unmapped reads
- `drop_2`: remove secondary/supplementary reads
- `drop_multipmap`: remove multimapped reads

All filters are disabled by default.


### `BigWig` file
```python
	alignment_filepath = '/path/to/alignment/fasta'
	genomewide_data_filepath = '/path/to/bigwig/file/'
	
	filtered_alignment_matrix, alignment_coordinate  = mapper.parse_and_filter(alignment_file=alignment_filepath,col_threshold = 0.10, col_content_threshold = 0.10, row_threshold = 0.10)
	bigwig_data_matrix=mapper.map_and_overlay(alignment_filepath, genomewide_data_filepath,data_format='bigwig', 
	col_threshold = 0.10, col_content_threshold = 0.10, row_threshold = 0.10)
```  
There is no additional parameter for `BigWig`.

### `VCF` file
```python
	alignment_filepath = '/path/to/alignment/fasta'
	genomewide_data_path = '/path/to/vcf/directory/'
	
	vcf_data_matrix=mapper.map_and_overlay(alignment_filepath, genomewide_data_path,data_format='vcf',query_key='AF',vcf_format='gnomad', 
	col_threshold = 0.10, col_content_threshold = 0.10, row_threshold = 0.10)
```
Here, `vcf_format` argument has two options, one for gnomAD dataset (fully supported) and another for direct file extraction. For gnomAD extraction, users have to provide directory path to gnomAD `VCF` files. Other than alternate allele frequency value, `query_key` argument also allows the users to access other fields.

### `MAF` file
```python
	alignment_filepath = '/path/to/alignment/fasta'
	genomewide_data_path = '/path/to/maf/directory/'
	
	maf_data_matrix=mapper.map_and_overlay(alignment_filepath, genomewide_data_path,data_format='maf',separated_maf=True, count_arg='common_freq', target_species='hg38', 
	col_threshold = 0.10, col_content_threshold = 0.10, row_threshold = 0.10)
```
Similarly, `separated_maf` toggle allows either per-chromosome `MAF` files or a single file extraction. For separated `MAF` extraction, users must provide directory path. The `target_species` argument sets the reference genome, for example, `hg38` for the human genome in the updated Zoonomia dataset. 

The `count_arg` parameter control the format at nucleotide level. Avialable options are:
- `ref_freq`: reference base frequency
- `coverage`: number of species aligned (not to be confused with sequencing coverage)
- `common_freq`: most common base frequency
- `common_freq_nogap`: same as above, excluding gaps from total count
- `common_raw`: raw count of the most common base
- `base_count`: all base counts
- `base_freq`: all base frequencies
- `base_freq_nogap`: all base frequencies excluding gaps
- `raw`: returns a `pandas.DataFrame` with alignment content at each position
- `raw_genome` same as `raw` but only with genome names only, excluding chromosome labels.


## Visualization
The `plots` module was used to visualize the resulting data matrix from the previous step. The code below is the intended simplest representation of the matrix.  

```python
	from ma_mapper import plots
	from ma_mapper import mapper
	plots.plot(
	data = [bigwig_data_matrix], 
	alignment=filtered_alignment_matrix,
	show_alignment=False, 
	heatmap_color=['RdBu'],
	heatmap_mode='overlay',
	heatmap_title=['Example of plot function'],
	heatmap_title_fs = 10, 
	heatmap_xlabel = 'position (bp)',
	heatmap_xlabel_fs = 10,
	heatmap_ylabel = 'sequences',
	heatmap_ylabel_fs = 10,
	vlim = [[-1,1]], 
	colorbar=True,
	colorbar_steps=[0.1]
	)
```  
<img src="img/2athe1cmas.png" alt="THE1C alignment" style="width:40%;" />

Slicing/Cropping for specific parts of the plot is possible using plot coordinates, as shown below:
```python
	from ma_mapper import plots
	from ma_mapper import mapper
	plots.plot(
	data = [bigwig_data_matrix], 
	alignment=filtered_alignment_matrix,
	show_alignment=False, 
	heatmap_color=['RdBu'],
	heatmap_mode='overlay',
	heatmap_title=['Example of plot function'],
	heatmap_title_fs = 10, 
	heatmap_xlabel = 'position (bp)',
	heatmap_xlabel_fs = 10,
	heatmap_ylabel = 'sequences',
	heatmap_ylabel_fs = 10,
	vlim = [[-1,1]], 
	colorbar=True,
	colorbar_steps=[0.1],
	xlim=[200,250],
	)
``` 

<img src="img/2athe1cmaszoom.png" alt="THE1C alignment" style="width:40%;" />

Additionally, the metadata table which contains TE coordinates, can be expanded further with additional information from either external sources or other data matrices transfromed into row-wise data arrays and subsequently converted into plot annotation 

See the tutorials [here](https://github.com/pkncsk/ma_mapper/blob/experimental/script/documentation/tutorial_nb/00overview.md) for more complex usages such as multiple matrix overlays or row annotation.
---

## Command line

Everything above is also available without writing Python. `--help` on any
subcommand lists every option with its default.

```bash
ma-mapper sequence  --coordinates copies.bed --genome hg38.fa --output copies.fasta
ma-mapper sequence  --repeatmasker hg38.rmsk.tsv --subfamily THE1C MER11A \
                    --genome hg38.fa --output copies.fasta
ma-mapper align     --input copies.fasta --threads 8 --output copies.aligned
ma-mapper overlay   --alignment copies.aligned --data signal.bw --format bigwig \
                    --output matrix.npy
ma-mapper plot      --matrix matrix.npy --kind profile --output profile.png
ma-mapper index-maf --maf chr1.maf              # once, before any parallel run
ma-mapper info                                  # versions, and which binaries are present
```

`overlay` writes three files, and the second and third are what make the first
interpretable a month later:

| file | holds |
|---|---|
| `matrix.npy` | the copies × columns matrix |
| `matrix.rows.tsv` | the coordinates of the copies that survived filtering, in row order |
| `matrix.alignment.npy` | the filtered alignment, so you can tell a gap from a measured zero |

The filter options (`--filter`, `--col-threshold`, `--col-content-threshold`,
`--row-threshold`) mirror `mapper.create_filter`, with the same defaults.
`--legacy-second-pass` reproduces the pre-0.9.0 column arithmetic, for
rebuilding an old matrix bit-for-bit.

## Parallelism

`extract_vcf`, `extract_maf` and the MAF index build are the only parallel
parts. They use as many worker processes as the process is allowed CPUs —
`os.sched_getaffinity`, which under Slurm is the allocation rather than the
machine — and respect, in order of precedence:

| variable | applies to |
|---|---|
| `MA_MAPPER_VCF_WORKERS` | `extract_vcf.vcf_io` |
| `MA_MAPPER_MAF_WORKERS` | `extract_maf.maf_io` |
| `MA_MAPPER_WORKERS` | both |

On a shared filesystem more workers mean more I/O, not more speed: the readers
are I/O bound, and a VCF as large as gnomAD is slow to open, not slow to query.

## Logging

Importing the package configures nothing and writes nothing. Warnings and
errors reach stderr through the standard library's fallback; for progress
output, ask for it:

```python
import ma_mapper
ma_mapper.configure_logging(logfile='run.log')   # or stream=True for stderr only
```

or set `MA_MAPPER_LOG=/path/to/run.log` before importing. On the command line,
`--verbose` and `--logfile` do the same.

## Semantics worth checking before publishing a number

| thing | what it actually does |
|---|---|
| `extract_bam(offset=...)` | **defaults to 5**, shifting every read 5 bp downstream |
| `row_threshold` | non-gap *density within a copy's own span* — correlated with element age, not a length filter |
| `normalise(denominator='alignment')` | divides by the copies with a base in the column, `nan` values included; `'measured'` divides by the copies with a finite value |
| `normalise(method='perc_coverage')` | counts non-zero cells, so a genuine 0.0 reads as absent |
| MAF `count_arg='coverage'` | counts species present in the block, deletions included |
| `generate_new_id=True` | one copy per row: a fragmented element is *not* reassembled |

## Reproducing results from before 0.9.0

Three fixes in 0.9.0 change numbers: minus-strand VCF extraction, BAM read
counts above 255, and the second column pass of `create_filter`. The first two
were wrong and have no reproduction switch. For the third, pass
`legacy_second_pass=True` to `create_filter` (or `--legacy-second-pass` on the
command line); that path reproduces the old masks exactly, which was checked on
MER11A, MER11C, LFSINE_Vert, MamRep434 and L1PA2 at three threshold settings.
See the [CHANGELOG](../CHANGELOG.md).
