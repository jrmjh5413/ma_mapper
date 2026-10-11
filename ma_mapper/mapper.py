"""The core: parse an alignment, filter it, and place extracted data in its columns.

The extractor modules are imported where they are used, not here, so that
parsing and filtering an alignment needs only numpy, pandas and Biopython --
a missing pysam or cyvcf2 then only affects the format that needs it.
"""

import inspect
import os
from functools import lru_cache
from typing import Callable, Literal, Tuple

import numpy as np
import pandas as pd

from . import logger
from . import sequence_alignment

#%%
def load_alignment_file(alignment_file):
    if isinstance(alignment_file, str) == True:
        from Bio import SeqIO
        records = (r for r in SeqIO.parse(alignment_file, "fasta"))
    else:
        records = alignment_file
    return records

def extract_coordinate_from_alignment(alignment_file, save_to_file = None):
    records = load_alignment_file(alignment_file)
    meta_id_list = []; chrom_list = []; start_list = []; end_list = []; strand_list = []; score_list= []
    import re
    for record in records:
        identifier=record.name
        match = re.match(r'([^:]+)::([^:]+):([\d,-]+)\((.)\)', identifier)
        if not match:
            raise ValueError(f"Invalid format: {identifier}")
        uniq_meta_id, chrom, coord_string, strand = match.groups()
        coord_pairs = coord_string.split(',')
        for pair in coord_pairs:
            start, end = map(int, pair.split('-'))
            chrom_list.append(chrom)
            start_list.append(start)
            end_list.append(end)
            meta_id_list.append(uniq_meta_id)
            score_list.append(20)
            strand_list.append(strand)

    new_coordinate_dict = {'chrom':chrom_list,'start':start_list,'end':end_list,'name':meta_id_list,'score':score_list,'strand':strand_list}
    coordinate_table = pd.DataFrame(new_coordinate_dict)
    if save_to_file:
        if isinstance(save_to_file, str) == True:
            output_filepath = save_to_file
        else:
            import os
            output_filepath =f'{os.path.dirname(os.path.abspath(__file__))}/alignment_coordinate.txt'
        coordinate_table.to_csv(output_filepath, sep='\t', index= False)
        logger.info('save new coordinate table at', output_filepath)
    else:
        return coordinate_table

def parse_alignment(alignment_file:str,
                    save_to_file: bool =False):
    """Turn an aligned FASTA into a copies x columns matrix of small integers.

    Gap is 0, A/C/T/G are 1/2/3/4 and N is 5, so ``matrix > 0`` is the
    non-gap mask every other function keys off.

    Parameters
    ----------
    alignment_file : str or sequence of SeqRecord
        Path to an aligned FASTA, or records already in memory. A one-shot
        iterator is materialised first: before 0.9.0 it was walked twice and
        the second pass, over an exhausted generator, left the whole matrix
        silently at zero -- an all-gap alignment.
    save_to_file : bool or str
        ``True`` writes ``<alignment_file>.parsed``; a string is used as the
        path. The HDF5 dataset is named ``parsed_alignment``.

    Returns
    -------
    numpy.ndarray
        uint8, shape (copies, columns). ``None`` if the matrix was saved.
    """
    logger.info('parse alignment')
    if isinstance(alignment_file, str):
        if not os.path.isfile(alignment_file):
            logger.error("Alignment file '%s' not found.", alignment_file)
            raise FileNotFoundError(f"Alignment file '{alignment_file}' not found.")
        records = load_alignment_file(alignment_file)
        source = alignment_file
    elif isinstance(alignment_file, (list, tuple)):
        records = alignment_file
        source = '<records>'
    else:
        # A generator can only be walked once; hold it so the two passes below
        # see the same records.
        records = list(alignment_file)
        source = '<records>'

    # Pass one: count the records and check they are all the same length.
    seq_count = 0
    seq_length = None
    ragged = None
    for value in records:
        if seq_count == 0:
            seq_length = len(value)
        elif len(value) != seq_length and ragged is None:
            ragged = (getattr(value, 'name', f'record {seq_count}'), len(value))
        seq_count += 1

    if seq_count == 0:
        # An empty/malformed alignment used to fall through to the np.zeros call
        # below and surface as UnboundLocalError('seq_length'), which hid the real
        # cause. 15 high-copy families (AluJb, MIR, L2c ...) were silently skipped
        # in the 2026-08-08 subfam_matrices_000 build for exactly this reason:
        # default MAFFT is O(N^2) and had written 0-byte .fasta.aligned files.
        raise ValueError(
            "alignment '%s' contains no records (empty or malformed file)" % source)
    if ragged is not None:
        raise ValueError(
            f"alignment '{source}' is not aligned: record '{ragged[0]}' is "
            f"{ragged[1]} columns long, the first record is {seq_length}"
        )

    #create parsed array
    parsed_alignment = np.zeros((seq_count, seq_length), dtype=np.uint8)
    if isinstance(alignment_file, str):
        records = load_alignment_file(alignment_file)
    for i, value in enumerate(records):
        parsed_alignment[i, :] = np.array(str(value.seq).upper(), 'c').view(np.uint8)

    # '-' and '.' both mean gap ('.' is HMMER's insert-column gap character).
    # Values are 0-5, so uint8 is sufficient: float64 costs 8x the memory for
    # no benefit, and np.bincount() in base_count() rejects float input.
    hash_table = {45:0, 46:0, 65:1, 67:2, 84:3, 71:4, 78:5}
    # np.zeros, not np.ndarray: np.ndarray leaves cells uninitialised, so any
    # character outside hash_table used to yield garbage that made the
    # `row > 0` test in map_data() non-deterministic.
    alignment_matrix = np.zeros(parsed_alignment.shape, dtype=np.uint8)
    for k in hash_table:
        alignment_matrix[parsed_alignment == k] = hash_table[k]
    # Unmapped characters silently become gaps; say so rather than hide it.
    # Counted with a 256-bin histogram rather than a boolean mask of the whole
    # matrix, which for a 223,000-copy family would allocate another 1.5 GB.
    char_histogram = np.bincount(parsed_alignment.ravel(), minlength=256)
    char_histogram[list(hash_table)] = 0
    unmapped_total = int(char_histogram.sum())
    if unmapped_total:
        chars = sorted(chr(c) for c in np.nonzero(char_histogram)[0])
        logger.warning(
            "%d alignment characters outside {-,.,A,C,T,G,N} were treated as gaps: %s",
            unmapped_total,
            chars,
        )
    if save_to_file:
        if isinstance(save_to_file, str) == True:
            output_file = save_to_file
        else:
            if isinstance(alignment_file, str):
                output_file = f'{alignment_file}.parsed'
            else:
                output_file = f'{os.path.dirname(os.path.abspath(__file__))}/alignment.parsed'

        import h5py
        with h5py.File(output_file, 'w') as hf:
            hf.create_dataset(_H5_DATASET,  data=alignment_matrix, compression="lzf")
        logger.info('done, saving parsed alignment at: %s', output_file)
    else:
        return alignment_matrix
#%%
_H5_DATASET = 'parsed_alignment'


def import_parsed_alignment(parsed_alignment_file):
    """Read back a matrix written by ``parse_alignment(save_to_file=...)``.

    Before 0.9.0 this looked for a dataset called ``result_array`` while the
    writer created ``parsed_alignment``, so every save/load round trip raised
    ``KeyError`` and no caller could use the saved form at all.
    """
    import h5py
    with h5py.File(f'{parsed_alignment_file}', 'r') as hf:
        if _H5_DATASET in hf:
            return hf[_H5_DATASET][:]
        if 'result_array' in hf:  # files written by other tools
            return hf['result_array'][:]
        raise KeyError(
            f"{parsed_alignment_file} holds no '{_H5_DATASET}' dataset "
            f"(found: {list(hf.keys())})"
        )
#%%
def _column_occupancy(parsed_alignment, legacy_span=False):
    """Per column: how many copies are non-gap there, and how many span it.

    ``col_counts_max`` counts the copies whose aligned span covers the column,
    gap or not, so ``col_counts / col_counts_max`` is occupancy among the
    copies that reach that far. With ``legacy_span`` the span excludes each
    copy's own last non-gap column, reproducing the pre-0.9.0 arithmetic in
    which a column could be non-gap in more copies than were counted as
    spanning it.
    """
    col_counts = np.zeros(parsed_alignment.shape[1])
    col_counts_max = np.zeros(parsed_alignment.shape[1])
    for row in parsed_alignment:
        nonzeros = np.nonzero(row)[0]
        if nonzeros.size == 0:
            # An all-gap row used to raise IndexError here.
            continue
        col_counts[nonzeros] += 1
        last = nonzeros[-1] if legacy_span else nonzeros[-1] + 1
        col_counts_max[nonzeros[0]:last] += 1
    return col_counts, col_counts_max


def create_filter(parsed_alignment, col_threshold = 0.50, col_content_threshold = 0.10, row_threshold = 0.50, save_to_file = False, legacy_second_pass = False):
    """Column and row masks for an alignment matrix.

    Three thresholds, applied as columns -> rows -> columns:

    ``col_threshold``
        keep a column whose occupancy among the copies spanning it exceeds this;
    ``col_content_threshold``
        and which is non-gap in at least this fraction of **all** copies;
    ``row_threshold``
        keep a copy whose non-gap bases, over the surviving columns, are at
        least this dense **within its own span**. This is a density, not a
        length: a short intact fragment passes and a full-length gappy copy
        fails, which makes the threshold correlate with element age.

    Parameters
    ----------
    parsed_alignment : numpy.ndarray or str
        Matrix from :func:`parse_alignment`, or a path to a saved one.
    save_to_file : bool or str
        Pickle the masks instead of returning them.
    legacy_second_pass : bool
        Reproduce the pre-0.9.0 second column pass, which *added* the
        surviving rows' occupancy to the counts already accumulated from all
        rows instead of recounting, and excluded each copy's last non-gap
        column from the span. Set this to rebuild matrices bit-for-bit as they
        were before 0.9.0.

    Returns
    -------
    list
        ``[row_filter, col_filter]``, both boolean masks.
    """
    if isinstance(parsed_alignment, str) == True:
        parsed_alignment = import_parsed_alignment(parsed_alignment)

    n_rows = parsed_alignment.shape[0]
    col_counts, col_counts_max = _column_occupancy(parsed_alignment, legacy_span=legacy_second_pass)
    col_filter = (col_counts > col_counts_max * col_threshold) & (col_counts >= n_rows * col_content_threshold)

    row_ratio = np.zeros(n_rows)
    for i, row in enumerate(parsed_alignment[:, col_filter]):
        nonzeros = np.nonzero(row)[0]
        if nonzeros.size and float(nonzeros[-1] - nonzeros[0]) > 0:
            row_ratio[i] = len(nonzeros) / (float(nonzeros[-1] - nonzeros[0]))
    row_filter = row_ratio >= row_threshold

    content_rows = n_rows
    if legacy_second_pass:
        # The pre-0.9.0 arithmetic: keep adding to the first pass's counters.
        extra_counts, extra_max = _column_occupancy(parsed_alignment[row_filter, :], legacy_span=True)
        col_counts = col_counts + extra_counts
        col_counts_max = col_counts_max + extra_max
    else:
        # Recount from the surviving copies only, which is what a second pass
        # is for. The old code added the surviving rows on top of the counts
        # from all rows, so the final column filter was computed from a mixture
        # of the two populations -- roughly doubling both terms of the
        # occupancy ratio while leaving the content threshold's denominator
        # alone, and mixing populations outright whenever the row filter
        # dropped a substantial share of copies.
        col_counts, col_counts_max = _column_occupancy(parsed_alignment[row_filter, :])
        # ... and judge "non-gap in at least this fraction of the copies"
        # against the copies that are still there. Keeping the original row
        # count here makes the threshold silently stricter in proportion to how
        # many copies the row filter removed: on MamRep434 at a row threshold
        # of 0.5, which keeps 130 of 2,102 copies, it leaves no columns at all.
        content_rows = max(1, int(row_filter.sum()))
    col_filter = (col_counts > col_counts_max * col_threshold) & (col_counts >= content_rows * col_content_threshold)
    filters = [row_filter,col_filter]
    if save_to_file:
        import pickle
        if isinstance(save_to_file, str) == True:
            output_file = save_to_file 
        else:
            output_file = f'{os.path.dirname(os.path.abspath(__file__))}/alignment.filter'
        pickle.dump(filters, open(output_file, "wb" ))
    else:
        return filters
#%%
import numpy as np

def iterative_filter(parsed_alignment_np, col_threshold=0.5, row_threshold=0.5, save_to_file=False):

    if isinstance(parsed_alignment_np, str):
        parsed_alignment_matrix = import_parsed_alignment(parsed_alignment_np)
    else:
        parsed_alignment_matrix = parsed_alignment_np
    num_rows, num_cols = parsed_alignment_matrix.shape

    # initialtion
    current_row_indices = np.arange(num_rows)
    current_col_indices = np.arange(num_cols)

    while True:
        submatrix = parsed_alignment_matrix[np.ix_(current_row_indices, current_col_indices)]

        # find nonzero column
        col_nonzero_counts = np.count_nonzero(submatrix, axis=0)
        col_nonzero_ratio = col_nonzero_counts / submatrix.shape[0]

        # Filter columns by threshold
        cols_to_keep_mask = col_nonzero_ratio >= col_threshold
        new_col_indices = current_col_indices[cols_to_keep_mask]

        # find nonzero row after applying column filter
        submatrix_rows_filtered = parsed_alignment_matrix[np.ix_(current_row_indices, new_col_indices)]
        row_nonzero_counts = np.count_nonzero(submatrix_rows_filtered, axis=1)
        row_nonzero_ratio = row_nonzero_counts / submatrix_rows_filtered.shape[1]

        # Filter rows by threshold
        rows_to_keep_mask = row_nonzero_ratio >= row_threshold
        new_row_indices = current_row_indices[rows_to_keep_mask]

        # check filter change
        if (np.array_equal(new_row_indices, current_row_indices) and
            np.array_equal(new_col_indices, current_col_indices)):
            break

        # Update indices for next iteration
        current_row_indices = new_row_indices
        current_col_indices = new_col_indices

    filters = [current_row_indices, current_col_indices]

    if save_to_file:
        import pickle
        if isinstance(save_to_file, str) == True:
            output_file = save_to_file 
        else:
            output_file = f'{os.path.dirname(os.path.abspath(__file__))}/alignment.filter'
        pickle.dump(filters, open(output_file, "wb" ))
    else:
        return filters

#%%
_NUMERIC_SCALARS = (int, float, np.integer, np.floating)


def _is_nested(extracted_data):
    """Does the extracted data hold objects (dicts of base counts) or numbers?

    A numeric payload goes on a float canvas, anything else on an object
    canvas. Only the first value of each copy is inspected: before 0.9.0 this
    walked every element of every copy in pure Python, which for a
    100,000-copy family meant hundreds of millions of iterations and dominated
    the runtime of the whole pipeline.
    """
    for row in extracted_data:
        if isinstance(row, np.ndarray):
            if row.dtype != object:
                return False
            if row.size == 0:
                continue
            sample = row.flat[0]
        else:
            if len(row) == 0:
                continue
            sample = row[0]
        if isinstance(sample, _NUMERIC_SCALARS):
            return False
    return True


def map_data(extracted_data:str,
             alignment_matrix: np.ndarray,
             filter: bool|None|list = True, **kwargs)-> np.ndarray:
    """Place extracted per-base values into their alignment columns.

    ``extracted_data[i]`` must hold exactly as many values as copy ``i`` has
    non-gap bases in ``alignment_matrix``, in the same 5'->3' order; the
    extractors guarantee this when both are built from the same coordinate
    table. A mismatch raises ``ValueError`` naming the copy.

    Parameters
    ----------
    extracted_data : list of array, or str
        Output of one of the ``*_io`` functions, or a path to a pickled one.
    alignment_matrix : numpy.ndarray
        Parsed alignment, from :func:`parse_alignment`.
    filter : bool, None or list
        ``True`` builds masks with :func:`create_filter` (any extra keyword
        arguments go to it), a list is used as ``[row_filter, col_filter]``,
        and ``False``/``None`` returns the unfiltered matrix.

    Returns
    -------
    numpy.ndarray
        float32 for numeric input, object dtype for dicts of base counts.
    """
    if isinstance(extracted_data, str) == True:
        import compress_pickle
        try:
            extracted_data = compress_pickle.load(extracted_data)
        except Exception:
            # compress_pickle guesses the compression from the extension; the
            # extractors write lzma whatever the name.
            extracted_data = compress_pickle.load(extracted_data, compression="lzma")

    nested_data = _is_nested(extracted_data)
    logger.info('nested_data:%s', nested_data)
    if nested_data:
        canvas = np.full(alignment_matrix.shape, None)
    else:
        # float32, not float64: genome-wide tracks (phyloP, CADD, allele
        # frequencies, coverage) carry far fewer than 7 significant digits, and
        # this array is the single largest allocation in the pipeline.
        canvas = np.zeros(alignment_matrix.shape, dtype=np.float32)
        
    if len(extracted_data) != alignment_matrix.shape[0]:
        raise ValueError(
            f"extracted data has {len(extracted_data)} copies but the alignment has "
            f"{alignment_matrix.shape[0]} rows: both must come from the same "
            "coordinate table, in the same order"
        )
    for idx, row in enumerate(alignment_matrix):
        non_gap = row > 0
        expected = int(non_gap.sum())
        values = extracted_data[idx]
        if len(values) != expected:
            raise ValueError(
                f"copy {idx} has {expected} non-gap bases in the alignment but "
                f"{len(values)} extracted values. The usual causes are an "
                "alignment and a coordinate table built from different copy "
                "sets, or an aligner that trimmed columns (hmmalign --trim) or "
                "wrote '.' for gaps without translating them to '-'."
            )
        canvas[idx, non_gap] = values
    mapped_data = canvas
    if filter: 
        if isinstance(filter, list):
            filters = filter
        else:  
            filters=create_filter(alignment_matrix, **kwargs)   
        row_filter = filters[0]
        col_filter = filters[1]
        mapped_data=mapped_data[np.ix_(row_filter,col_filter)]

    return mapped_data
#%%
def base_count(alignment_matrix:pd.DataFrame):
    """Count each symbol per alignment column.

    Returns a (6, columns) array of counts indexed by the parse codes
    0 gap, 1 A, 2 C, 3 T, 4 G, 5 N.

    Before 0.9.0 this ran ``np.bincount(col, minlength=5)`` column by column
    through ``apply_along_axis``: the first column fixed the output height, so
    an alignment whose first column held no N but whose later columns did
    failed outright with a broadcast error. Counting the six symbols with
    vectorised comparisons also removes a Python call per column.
    """
    return np.stack(
        [np.count_nonzero(alignment_matrix == code, axis=0) for code in range(6)]
    )

_NORM_METHOD = Literal['average','perc_coverage','median']
_NORM_DENOMINATOR = Literal['alignment','measured']
def normalise(alignment_matrix: str|pd.DataFrame,
              data_matrix: pd.DataFrame,
              method:_NORM_METHOD = 'average',
              denominator:_NORM_DENOMINATOR = 'alignment',
              ):
    """Collapse a copies x columns matrix to one value per column.

    Parameters
    ----------
    alignment_matrix : numpy.ndarray
        The parsed alignment, used to count how many copies have a base in
        each column.
    data_matrix : numpy.ndarray
        The overlay, same shape.
    method : {'average', 'perc_coverage', 'median'}
        ``average`` sums the data and divides by the denominator below;
        ``perc_coverage`` reports the percentage of copies with a non-zero
        value; ``median`` ignores the alignment entirely and takes
        ``nanmedian`` down each column.
    denominator : {'alignment', 'measured'}
        What ``average`` and ``perc_coverage`` divide by. ``alignment`` (the
        default, and the only behaviour before 0.9.0) divides by the number of
        copies with a base in the column, **including copies whose data is
        nan** -- so a track that does not cover part of the family reads as a
        low value rather than as missing. ``measured`` divides by the number
        of copies that actually have a finite value there. Which one is right
        depends on whether a gap in the data means "zero" or "not measured";
        for bigWig input it means "not measured".

    Notes
    -----
    ``perc_coverage`` counts a genuine 0.0 as absent, because it uses
    ``count_nonzero``: that is only a coverage measure for data where 0 means
    "nothing here", such as peak overlap.

    A column no copy occupies yields ``nan`` rather than ``inf``.
    """
    data_matrix = np.asarray(data_matrix, dtype=float)
    finite = np.isfinite(data_matrix)

    if method == 'median':
        return np.nanmedian(data_matrix, axis=0)

    if denominator == 'measured':
        normalizer = np.count_nonzero(finite & (np.asarray(alignment_matrix) > 0), axis=0)
    elif denominator == 'alignment':
        normalizer = np.count_nonzero(alignment_matrix, axis=0)
    else:
        raise ValueError(f"denominator must be 'alignment' or 'measured', got {denominator!r}")

    if method == 'average':
        numerator = np.nansum(data_matrix, axis=0)
    elif method == 'perc_coverage':
        numerator = np.count_nonzero(np.nan_to_num(data_matrix), axis=0) * 100.0
    else:
        raise ValueError(
            f"method must be 'average', 'perc_coverage' or 'median', got {method!r}"
        )

    normalizer = np.asarray(normalizer, dtype=float)
    return np.divide(
        numerator,
        normalizer,
        out=np.full(numerator.shape, np.nan, dtype=float),
        where=normalizer > 0,
    )

def flank_sequence_io(coordinate_table: pd.DataFrame,
                    source_fasta: str,
                    coordinate_table_out:bool = False, 
                    extension_length:int = 500) -> Tuple[pd.DataFrame, pd.DataFrame]:

    front_coordinate_table = coordinate_table.reset_index(drop=True).copy()
    front_coordinate_table.end = front_coordinate_table.start.astype(int)
    # Clamp at the start of the chromosome: a copy within extension_length of
    # position 0 used to be given a negative start, which every extractor then
    # read as an invalid region. The far end cannot be clamped here because
    # the chromosome length is not known to this function.
    front_coordinate_table.start = (
        front_coordinate_table.start.astype(int) - extension_length
    ).clip(lower=0)

    back_coordinate_table = coordinate_table.reset_index(drop=True).copy()
    back_coordinate_table.start = back_coordinate_table.end.astype(int)
    back_coordinate_table.end = back_coordinate_table.end.astype(int) + extension_length
    if coordinate_table_out:
        return front_coordinate_table, back_coordinate_table
    else:
        front_seq = sequence_alignment.sequence_io(front_coordinate_table, source_fasta)
        back_seq = sequence_alignment.sequence_io(back_coordinate_table, source_fasta)
        
        front_list = []
        back_list = []

        for i, row in coordinate_table.reset_index().iterrows():
            if row.strand == '+':
                front_list.append(front_seq[i])
                back_list.append(back_seq[i])
            else:
                front_list.append(back_seq[i])
                back_list.append(front_seq[i])
        front_parsed_alignment=parse_alignment(front_list)
        back_parsed_alignment =parse_alignment(back_list)
        return front_parsed_alignment, back_parsed_alignment
#steamline functions

def parse_and_filter(alignment_file: str,
                     filter: Literal['one_pass', 'iterative', None] = 'one_pass',
                     extension_length:int | None = None,
                     preprocess_out:bool = False,
                     **kwargs)-> Tuple[np.ndarray, pd.DataFrame]:
    if filter not in ('one_pass', 'iterative', None):
        raise ValueError(f"Invalid value for `filter`: {filter!r}. Must be 'one_pass', 'iterative', or None.")
    parsed_alignment = parse_alignment(alignment_file)
    coordinate_table = extract_coordinate_from_alignment(alignment_file)
    coordinate_table['original_order'] = coordinate_table.index
    # How an alignment row maps to coordinate rows depends on how the FASTA was
    # written, and both conventions are in use:
    #
    #   one record per copy      -- sequence_io joins a copy's fragments into a
    #                               single record, whose header carries all of
    #                               them, so one alignment row can own several
    #                               coordinate rows and names are unique;
    #   one record per fragment  -- one coordinate row per record, names repeated
    #                               across the fragments of one copy.
    #
    # Positional pairing (what this function did before 0.9.0) is right for the
    # second and wrong for the first, where it raised
    # "IndexError: Boolean index has wrong length".
    n_rows_alignment = parsed_alignment.shape[0]
    copy_names = coordinate_table['name'].drop_duplicates().to_numpy()
    if len(coordinate_table) == n_rows_alignment:
        per_fragment_records = True
    elif len(copy_names) == n_rows_alignment:
        per_fragment_records = False
    else:
        raise ValueError(
            f"{n_rows_alignment} alignment rows, but their headers give "
            f"{len(coordinate_table)} coordinate rows over {len(copy_names)} distinct "
            "names: the alignment and its coordinates cannot be matched up, either "
            "one row per record or one row per name"
        )
    if filter:
        if filter == 'one_pass':
            filters=create_filter(parsed_alignment, **filter_kwargs(create_filter, kwargs))
        if filter == 'iterative':
            filters=iterative_filter(parsed_alignment, **filter_kwargs(iterative_filter, kwargs))
        row_filter = filters[0]
        col_filter = filters[1]
        parsed_alignment_filtered=parsed_alignment[np.ix_(row_filter,col_filter)]
        if per_fragment_records:
            coordinate_table_filtered=coordinate_table.iloc[row_filter,:]
        else:
            surviving_names = copy_names[row_filter]
            coordinate_table_filtered=coordinate_table[
                coordinate_table['name'].isin(surviving_names)
            ]
    else:
        filters = [np.arange(parsed_alignment.shape[0]), np.arange(parsed_alignment.shape[1])]
        parsed_alignment_filtered = parsed_alignment
        coordinate_table_filtered = coordinate_table
    #print(metadata_aligned_filtered.dtypes)
    if extension_length is not None:
        ffs_kwargs = {key: value for key, value in kwargs.items() if key in ('source_fasta')}
        front_parsed_alignment, back_parsed_alignment= flank_sequence_io(coordinate_table=coordinate_table_filtered, extension_length=extension_length, **ffs_kwargs)
        alignment_output = np.hstack((front_parsed_alignment, parsed_alignment_filtered, back_parsed_alignment))
    else:
        alignment_output = parsed_alignment_filtered
    
    if preprocess_out:
        return parsed_alignment, coordinate_table, filters
    else:
        return alignment_output, coordinate_table_filtered

def match_age_to_id_coordinate(coordinate_file: pd.DataFrame|str, 
                             age_table:str|None = None) -> pd.DataFrame:
    if isinstance(coordinate_file, str):
        if os.path.isfile(coordinate_file):
            coordinate_df = pd.read_csv(coordinate_file, sep='\t')
        else:
            logger.error('coordinate file not found')
    elif isinstance(coordinate_file, pd.DataFrame):
        coordinate_df = coordinate_file
    if isinstance(age_table,list):
        age_df_list = []
        for age_tbl in age_table:
            age_df_list.append(pd.read_csv(age_tbl, sep='\t'))
        age_df=pd.concat(age_df_list)
    elif isinstance(age_table, str):
        if os.path.isfile(age_table):
            age_df = pd.read_csv(age_table, sep='\t')
        else:
            logger.error('age table file not found')
            raise FileNotFoundError(f"age table file not found")
    elif isinstance(age_table, pd.DataFrame):
        age_df = age_table
    #assume BED format:
    bed_headers = ['chrom','start','end','name','score','strand']
    if coordinate_df.shape[1] == 6:
        coordinate_df.columns = bed_headers
    elif coordinate_df.shape[1] > 6:
        current_columns = coordinate_df.columns.tolist()
        coordinate_df.columns = bed_headers + current_columns[6:]
    merged_table=coordinate_df.merge(age_df, left_on='name', right_on='internal_id', how='left')
    merged_table=merged_table.drop(columns=['internal_id'])

    return merged_table

@lru_cache
def get_func_param_names(func: Callable):
    return set(inspect.signature(func).parameters)

def filter_kwargs(funcs: list[Callable] | Callable, kwargs: dict) -> dict:
    if not isinstance(funcs, list):
        funcs = [funcs]
    accepted_keys = set().union(*(get_func_param_names(f) for f in funcs))
    return {k: v for k, v in kwargs.items() if k in accepted_keys}

def _extract(coordinate_table, data_file, data_format, kwargs):
    """Run the extractor that matches ``data_format``.

    One dispatch for the element itself and for its flanks. Before 0.9.0 the
    two were written out separately, and the flank copy referenced
    ``extract_bam`` without importing it.
    """
    if data_format in _BAM_FORMATS:
        from . import extract_bam
        io_kwargs = filter_kwargs(extract_bam.bam_io, kwargs)
        return extract_bam.bam_io(
            coordinate_table=coordinate_table,
            bam_file=data_file,
            bam_format=data_format,
            **io_kwargs,
        )
    if data_format == 'bigwig':
        from . import extract_bigwig
        io_kwargs = filter_kwargs(extract_bigwig.bigwig_io, kwargs)
        return extract_bigwig.bigwig_io(
            coordinate_table=coordinate_table, bigwig=data_file, **io_kwargs
        )
    if data_format == 'bed':
        from . import extract_bed
        io_kwargs = filter_kwargs(extract_bed.bed_io, kwargs)
        return extract_bed.bed_io(
            coordinate_table=coordinate_table, bed=data_file, **io_kwargs
        )
    if data_format == 'maf':
        from . import extract_maf
        io_kwargs = filter_kwargs(extract_maf.maf_io, kwargs)
        return extract_maf.maf_io(
            coordinate_table=coordinate_table, maf=data_file, **io_kwargs
        )
    if data_format == 'vcf':
        from . import extract_vcf
        io_kwargs = filter_kwargs(extract_vcf.vcf_io, kwargs)
        return extract_vcf.vcf_io(
            coordinate_table=coordinate_table, vcf=data_file, **io_kwargs
        )
    raise ValueError(f"data_format must be one of {DATA_FORMATS}, got {data_format!r}")


_BAM_FORMATS = ('read_max','read_min','read_forward','read_reverse','read_sum',
                'normal_max','normal_min','normal_forward','normal_reverse')
_FORMAT = Literal['read_max','read_min','read_forward','read_reverse','read_sum',
                  'normal_max','normal_min','normal_forward','normal_reverse',
                  'bigwig','bed','maf','vcf']
DATA_FORMATS = _BAM_FORMATS + ('bigwig', 'bed', 'maf', 'vcf')


def map_and_overlay(alignment:str,
            data_file:str,
            data_format:_FORMAT,
            filter:bool = True,
            coordinate_file:str|None = None,
            extension_length: int|None = None,
            source_fasta:str = None,
            filter_mode:Literal['one_pass','iterative'] = 'one_pass',
            return_details:bool = False,
            **kwargs,
            ):
    """Alignment plus dataset in, copies x columns matrix out.

    Parameters
    ----------
    alignment : str
        Aligned FASTA of one family, with headers written by
        ``sequence_alignment.sequence_io`` (``name::chrom:start-end(strand)``),
        which is where the genomic coordinates come from.
    data_file : str
        The dataset: a BAM, bigWig, BED, MAF or VCF, matching ``data_format``.
    data_format : str
        One of ``DATA_FORMATS``. Before 0.9.0 the annotation advertised
        ``bam_max`` and friends, which the dispatch below did not accept: the
        call logged "dataformat field unspecified" and then failed with
        ``UnboundLocalError``.
    filter : bool
        Apply the filter; extra keyword arguments reach :func:`create_filter`.
    filter_mode : {'one_pass', 'iterative'}
        Which filter to use. Before 0.9.0 the iterative filter could not be
        reached through this function at all: ``parse_and_filter`` was always
        called with its own default.
    return_details : bool
        Return ``(matrix, coordinate_table, alignment_matrix)`` instead of just
        the matrix: the coordinates of the surviving copies in matrix row
        order, and the filtered alignment itself. Without them a matrix has no
        row labels and no record of which bases are gaps.
    coordinate_file : str or None
        Use these coordinates instead of the ones in the alignment headers.
        Rows are matched to alignment rows by name.
    extension_length : int or None
        Also extract this many bases of flanking sequence on each side.
        Requires ``source_fasta``.
    source_fasta : str
        Genome FASTA, needed only with ``extension_length``.
    **kwargs
        Passed to the filter and the extractor; prefix with ``md_`` for
        :func:`map_data` or ``fs_`` for the flank extraction.

    Returns
    -------
    numpy.ndarray
        The overlay matrix, filtered unless ``filter=False``.
    """
    if data_format not in DATA_FORMATS:
        raise ValueError(
            f"data_format must be one of {DATA_FORMATS}, got {data_format!r}"
        )
    combined_kwargs = filter_kwargs([parse_and_filter, create_filter], kwargs)
    fs_kwargs = {key[3:]: value for key, value in kwargs.items() if key.startswith('fs_')}   
    md_kwargs = {key[3:]: value for key, value in kwargs.items() if key.startswith('md_')}
    kwargs = {key: value for key, value in kwargs.items() if not key.startswith(('pf_','md_','fs_'))}
    effective_filter_mode = filter_mode if filter else None
    alignment_matrix, alignment_coordinate, filters  = parse_and_filter(alignment_file=alignment, filter=effective_filter_mode, preprocess_out=True, **combined_kwargs)
    #use alignment coordinate as coordinate file if there is no coordinate file input
    if coordinate_file is None:
        coordinate_table = alignment_coordinate
    elif isinstance(coordinate_file, str):
        # Before 0.9.0 this tested os.path.isfile(coordinate_table), a name
        # that does not exist yet at this point, so any coordinate_file path
        # raised NameError.
        if not os.path.isfile(coordinate_file):
            raise FileNotFoundError(f"coordinate file not found: {coordinate_file}")
        coordinate_table = pd.read_csv(coordinate_file, sep='\t', header=None)
    else:
        coordinate_table = coordinate_file
    extracted_data = _extract(coordinate_table, data_file, data_format, kwargs)
    row_filter, col_filter = filters

    if coordinate_file is None:
        alignment_matrix_sorted = alignment_matrix
        alignment_coordinate_sorted = alignment_coordinate
        row_filter_sorted = row_filter
    else: #if coordinate table is specified, alignment will be rearranged by NAME field in the coordinate table
        # The extractor returns one array per copy in the coordinate table's
        # order, so the alignment rows have to be permuted INTO that order.
        # Before 0.9.0 the permutation was built the other way round -- the
        # position of each alignment row in the coordinate table, used directly
        # as an index -- which is the inverse permutation. It happens to be
        # correct for an identity or a straight swap, and silently pairs the
        # wrong copy with the wrong data for anything else; it also walked the
        # coordinate table fragment by fragment, so a fragmented copy produced
        # more indices than the alignment has rows.
        if kwargs.get('generate_new_id'):
            raise ValueError(
                "generate_new_id cannot be combined with coordinate_file: the rows "
                "are matched to the alignment by name"
            )
        source_order = coordinate_table.iloc[:,3].unique()
        alignment_names = alignment_coordinate['name'].drop_duplicates().to_numpy()
        if len(alignment_names) != alignment_matrix.shape[0]:
            raise ValueError(
                f"the alignment has {alignment_matrix.shape[0]} rows but only "
                f"{len(alignment_names)} distinct names in its headers, so rows cannot "
                "be matched to a coordinate_file by name. Give one record per copy "
                "(as sequence_alignment.sequence_io writes), or drop coordinate_file "
                "and use the coordinates in the headers."
            )
        name_to_row = {name: index for index, name in enumerate(alignment_names)}
        missing = [name for name in source_order if name not in name_to_row]
        if missing:
            raise ValueError(
                f"{len(missing)} copies in the coordinate table are not in the "
                f"alignment, first: {missing[0]}"
            )
        unused = [name for name in alignment_names if name not in set(source_order)]
        if unused:
            raise ValueError(
                f"{len(unused)} copies in the alignment are not in the coordinate "
                f"table, first: {unused[0]}"
            )
        sorted_order = [name_to_row[name] for name in source_order]
        alignment_matrix_sorted = alignment_matrix[sorted_order]
        alignment_coordinate_sorted = (
            alignment_coordinate.set_index('name').loc[list(source_order)].reset_index()
        )
        row_filter_sorted = row_filter[sorted_order]
    if filter:
        if isinstance(filter, list):
            mapping_filters = filter
        else:
            mapping_filters = [row_filter_sorted, col_filter]
        # Filtered by name, not by position: a fragmented copy has more
        # coordinate rows than the alignment has alignment rows.
        sorted_names = alignment_coordinate_sorted['name'].drop_duplicates().to_numpy()
        surviving_names = set(sorted_names[row_filter_sorted])
        alignment_coordinate_filtered=alignment_coordinate_sorted[
            alignment_coordinate_sorted['name'].isin(surviving_names)
        ]
    else:
        mapping_filters = None
        alignment_coordinate_filtered=alignment_coordinate_sorted
    output_matrix_filtered = map_data(extracted_data=extracted_data, alignment_matrix=alignment_matrix_sorted, filter=mapping_filters, **md_kwargs)
    
    
    if extension_length is not None:
        
        front_coordinate_table, back_coordinate_table = flank_sequence_io(alignment_coordinate_filtered, coordinate_table_out=True, extension_length=extension_length,source_fasta=source_fasta, **fs_kwargs)
        front_extracted_data = _extract(front_coordinate_table, data_file, data_format, kwargs)
        back_extracted_data = _extract(back_coordinate_table, data_file, data_format, kwargs)
        logger.info(
            'flanks: %s rows, %d front, %d back',
            alignment_coordinate_sorted.shape,
            len(front_extracted_data),
            len(back_extracted_data),
        )
        output_matrix_filtered = np.hstack((front_extracted_data,output_matrix_filtered,back_extracted_data))

    if return_details:
        if mapping_filters is None:
            alignment_matrix_out = alignment_matrix_sorted
        else:
            alignment_matrix_out = alignment_matrix_sorted[
                np.ix_(mapping_filters[0], mapping_filters[1])
            ]
        return output_matrix_filtered, alignment_coordinate_filtered, alignment_matrix_out
    return output_matrix_filtered