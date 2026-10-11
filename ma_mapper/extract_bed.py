"""Paint the score of overlapping BED features onto each base of each copy.

Every base of every copy gets the sum of the scores of the features covering
it, in the element's own orientation. Bases with no feature stay 0.0.
"""

import os

import numpy as np
import pandas as pd
import pybedtools

from . import logger
from ._io_common import group_coordinates, read_coordinate_table, save_or_return

# Column layout of the loj intersect output: six columns of the query copy,
# then the overlapping feature. Read by position, because to_dataframe()'s
# BED12 auto-names only line up for a feature file of exactly six columns.
_A_CHROM, _A_START, _A_END, _A_NAME, _A_SCORE, _A_STRAND = range(6)
_B_START, _B_END, _B_SCORE = 7, 8, 10


def coordinate_table_to_bed(coordinate_table, save_to_file=False, export_as_bedtool=False):
    """Reduce a coordinate table to a six-column BED.

    Parameters
    ----------
    coordinate_table : pandas.DataFrame or str
        BED-like table. Columns are taken by position, so a frame read without
        a header works; before 0.9.0 this needed named ``meta_id`` and
        ``strand`` columns and raised ``AttributeError`` otherwise.
    save_to_file : bool or str
        ``True`` writes ``coordinate.bed`` beside the input; a string is used
        as the path.
    export_as_bedtool : bool
        Return a ``pybedtools.BedTool`` instead of a DataFrame.
    """
    coordinate_local = read_coordinate_table(coordinate_table)
    coordinate_bed = coordinate_local.iloc[:, 0:3].copy()
    coordinate_bed.columns = ["chrom", "start", "end"]
    name_column = "meta_id" if "meta_id" in coordinate_local.columns else coordinate_local.columns[3]
    coordinate_bed["name"] = coordinate_local[name_column].values
    coordinate_bed["score"] = 10
    strand_column = "strand" if "strand" in coordinate_local.columns else coordinate_local.columns[5]
    coordinate_bed["strand"] = coordinate_local[strand_column].values

    if save_to_file:
        # Before 0.9.0 this branch called str.split on the DataFrame instead of
        # the path and raised TypeError.
        if isinstance(save_to_file, str):
            output_filepath = save_to_file
        elif isinstance(coordinate_table, str):
            output_filepath = os.path.join(
                os.path.dirname(os.path.abspath(coordinate_table)), "coordinate.bed"
            )
        else:
            output_filepath = os.path.join(os.getcwd(), "coordinate.bed")
        coordinate_bed.to_csv(output_filepath, sep="\t", index=False, header=False)
        logger.info("saved coordinate bed at: %s", output_filepath)

    if export_as_bedtool:
        return pybedtools.BedTool.from_dataframe(coordinate_bed)
    return coordinate_bed


def _feature_score(raw):
    """Score of an overlapping feature as a float.

    A BED score may be an integer, a float or ``'.'`` for "no score". Before
    0.9.0 this was ``int(...)``, which truncated 0.87 to 0 and raised
    ``ValueError`` on ``'.'``; a missing score now counts as 1, i.e. "this base
    is covered".
    """
    if raw is None:
        return 1.0
    text = str(raw).strip()
    if text in (".", "", "nan"):
        return 1.0
    return float(text)


def _as_int(raw, default=-1):
    """Integer coordinate from an intersect field, tolerating bedtools' '.' filler."""
    text = str(raw).strip()
    try:
        return int(text)
    except ValueError:
        return default


def _paint_copy(group, strand):
    """Per-base feature scores for the fragments of one copy, element orientation."""
    drawings = []
    for canvas_start, fragment in group:
        canvas_end = fragment["a_end"]
        canvas = np.zeros(canvas_end - canvas_start, dtype=float)
        for b_start, b_end, b_score in fragment["hits"]:
            # cap out-of-bound positions
            intersect_start = max(b_start, canvas_start)
            intersect_end = min(b_end, canvas_end)
            if intersect_end <= intersect_start:
                continue
            if strand == "+":
                left = intersect_start - canvas_start
                right = intersect_end - canvas_start
            else:
                left = canvas_end - intersect_end
                right = canvas_end - intersect_start
            canvas[left:right] += b_score
        drawings.append(canvas)
    if strand == "-":
        drawings.reverse()
    return np.concatenate(drawings) if drawings else np.zeros(0, dtype=float)


def _group_intersect(intersect_df):
    """Index the loj intersect output by copy name.

    Returns ``{name: (strand, [(fragment_start, {...}), ...])}``. Built once for
    the whole table: before 0.9.0 the full intersect DataFrame was pickled and
    shipped to a worker process once per copy, so a family of 100,000 copies
    serialised the same table 100,000 times.
    """
    by_name = {}
    columns = intersect_df.columns
    for name, rows in intersect_df.groupby(columns[_A_NAME], sort=False):
        strand = rows.iloc[0, _A_STRAND]
        fragments = {}
        for row in rows.itertuples(index=False, name=None):
            key = _as_int(row[_A_START])
            a_end = _as_int(row[_A_END])
            fragment = fragments.get(key)
            if fragment is None:
                fragment = fragments[key] = {"a_end": a_end, "hits": []}
            elif fragment["a_end"] != a_end:
                logger.info("multiple length frags for %s at %s", name, key)
            b_start = _as_int(row[_B_START])
            b_end = _as_int(row[_B_END])
            if b_start < 0 or b_end < 0:  # loj filler: no overlap
                continue
            fragment["hits"].append((b_start, b_end, _feature_score(row[_B_SCORE])))
        by_name[name] = (strand, sorted(fragments.items()))
    return by_name


def intersect_process(intersect_df, unique_id):
    """Per-base feature scores for one copy. Kept for backwards compatibility.

    ``bed_io`` no longer calls this: it indexes the intersect table once with
    ``_group_intersect`` instead of re-scanning it per copy.
    """
    grouped = _group_intersect(intersect_df)
    if unique_id not in grouped:
        raise KeyError(f"{unique_id} is not in the intersect table")
    strand, fragments = grouped[unique_id]
    return _paint_copy(fragments, strand)


def bed_io(
    coordinate_table,
    bed,
    save_to_file=False,
    generate_new_id=False,
    strand_overlap=True,
):
    """Overlay a BED file on every copy in a coordinate table.

    Parameters
    ----------
    coordinate_table : pandas.DataFrame or str
        BED-like table, one or more rows per copy.
    bed : str, pandas.DataFrame or pybedtools.BedTool
        Features to overlay. A DataFrame is used as its first six columns;
        before 0.9.0 a DataFrame raised ``TypeError`` on a leftover second
        conversion of the same argument.
    save_to_file : bool or str
        ``True`` writes ``bed_out.p`` beside the coordinate table; a string is
        used as the output path.
    generate_new_id : bool
        Treat every row as its own copy instead of grouping by column 3.
    strand_overlap : bool
        Require the feature to be on the same strand as the copy
        (``bedtools intersect -s``).

    Returns
    -------
    list of numpy.ndarray
        One array per copy, or the output path if the result was saved.
    """
    coordinate_local = read_coordinate_table(coordinate_table)
    meta_id, chrom_list, start_list, end_list, strand_list = group_coordinates(
        coordinate_local, generate_new_id=generate_new_id
    )

    if isinstance(bed, str):
        if not os.path.isfile(bed):
            raise FileNotFoundError(f"bed file not found: {bed}")
        bed_file = pybedtools.BedTool(bed)
        logger.info("extract from bed file: %s", bed)
    elif isinstance(bed, pd.DataFrame):
        bed_file = pybedtools.BedTool.from_dataframe(bed.iloc[:, 0:6])
        logger.info("extract from dataframe")
    elif isinstance(bed, pybedtools.BedTool):
        bed_file = bed
    else:
        raise TypeError(
            f"bed must be a path, DataFrame or BedTool, got {type(bed).__name__}"
        )

    query = coordinate_local.iloc[:, 0:6].copy()
    if generate_new_id:
        query.iloc[:, _A_NAME] = [f"entry_{index}" for index in coordinate_local.index.astype(str)]
    coordinate_bed = pybedtools.BedTool.from_dataframe(query)

    intersect_bed = coordinate_bed.intersect(bed_file, loj=True, wa=True, wb=True, s=strand_overlap)
    # Positional column names: the BED12 auto-names only line up when the
    # feature file has exactly six columns.
    intersect_df = intersect_bed.to_dataframe(disable_auto_names=True, header=None)
    if intersect_df.empty or intersect_df.shape[1] <= _B_SCORE:
        raise ValueError(
            "bedtools intersect returned nothing usable: the feature file needs at "
            "least six BED columns, and the query must not be empty"
        )

    grouped = _group_intersect(intersect_df)
    missing = [name for name in meta_id if name not in grouped]
    if missing:
        raise ValueError(
            f"{len(missing)} copies are absent from the intersect output, first: {missing[0]}"
        )

    logger.info("painting %d copies", len(meta_id))
    bed_out = [_paint_copy(grouped[name][1], grouped[name][0]) for name in meta_id]

    return save_or_return(bed_out, save_to_file, coordinate_table, "bed_out.lzma", "bed_out")
