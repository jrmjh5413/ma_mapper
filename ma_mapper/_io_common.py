"""Shared input handling for the five ``*_io`` entry points.

Every extractor (BED, BAM, bigWig, VCF, MAF) starts from the same coordinate
table and ends with the same save-or-return decision. Until 0.9.0 each of them
carried its own copy of that code, and therefore its own copy of the same three
defects: ``if save_to_file == True`` silently ignored a path (``'x' == True`` is
``False`` in Python, so a caller asking to save got nothing written and no
error), a missing coordinate file was logged and then crashed with
``NameError`` further down, and the per-copy grouping ran four separate
``groupby().apply()`` passes over the table.

The coordinate table is a BED-like frame: column 0 chromosome, 1 start
(0-based), 2 end (exclusive), 3 name, 4 score, 5 strand. One copy may occupy
several rows (one per fragment of a broken-up element); rows are grouped by
name.
"""

import os

import pandas as pd

from . import logger

BED_COLUMNS = ("chrom", "start", "end", "name", "score", "strand")


def read_coordinate_table(coordinate_table):
    """Return the coordinate table as a DataFrame.

    Accepts a DataFrame (returned as-is) or a path to a headerless
    tab-separated BED-like file. A missing path raises ``FileNotFoundError``
    rather than logging an error and failing later with ``NameError``.
    """
    if isinstance(coordinate_table, str):
        if not os.path.isfile(coordinate_table):
            raise FileNotFoundError(f"coordinate_table file not found: {coordinate_table}")
        return pd.read_csv(coordinate_table, sep="\t", header=None)
    if isinstance(coordinate_table, pd.DataFrame):
        return coordinate_table
    raise TypeError(
        f"coordinate_table must be a path or a pandas DataFrame, got {type(coordinate_table).__name__}"
    )


def group_coordinates(coordinate_table, generate_new_id=False, sort_blocks=True):
    """Group a coordinate table into one entry per copy.

    Returns ``(meta_id, chrom_list, start_list, end_list, strand_list)``, where
    ``meta_id`` is in order of first appearance and ``start_list[i]`` /
    ``end_list[i]`` are the fragment coordinates of copy ``i``.

    Parameters
    ----------
    generate_new_id : bool
        Number the copies ``entry_0, entry_1, ...`` instead of using column 3.
        Each row becomes its own copy, so a fragmented element is not reassembled.
    sort_blocks : bool
        Sort the fragments of each copy by start position. The minus-strand
        handling in every extractor reverses the fragment order, which is only
        correct if the fragments arrive in ascending genomic order; an unsorted
        table would otherwise produce a silently scrambled copy.

    The input frame is not modified: before 0.9.0 these functions added a
    ``meta_id`` column to the caller's DataFrame as a side effect.
    """
    local = coordinate_table.copy()
    if local.shape[1] < 6:
        raise ValueError(
            "coordinate table needs at least 6 BED columns "
            f"(chrom, start, end, name, score, strand); got {local.shape[1]}"
        )

    if generate_new_id:
        local["meta_id"] = [f"entry_{index}" for index in local.index.astype(str)]
    else:
        local["meta_id"] = local.iloc[:, 3]

    chrom_list, start_list, end_list, strand_list, meta_id = [], [], [], [], []
    unsorted_copies = 0
    mixed_chrom = 0
    mixed_strand = 0

    # One pass. The previous four groupby().apply() calls each walked the whole
    # table, and apply() on a DataFrameGroupBy is an order of magnitude slower
    # than iterating the groups.
    for key, group in local.groupby("meta_id", sort=False):
        starts = group.iloc[:, 1].tolist()
        ends = group.iloc[:, 2].tolist()
        if len(starts) > 1 and any(b < a for a, b in zip(starts, starts[1:])):
            unsorted_copies += 1
            if sort_blocks:
                order = sorted(range(len(starts)), key=lambda i: starts[i])
                starts = [starts[i] for i in order]
                ends = [ends[i] for i in order]
        if group.iloc[:, 0].nunique() > 1:
            mixed_chrom += 1
        if group.iloc[:, 5].nunique() > 1:
            mixed_strand += 1
        meta_id.append(key)
        chrom_list.append(group.iloc[0, 0])
        start_list.append(starts)
        end_list.append(ends)
        strand_list.append(group.iloc[0, 5])

    if unsorted_copies:
        logger.warning(
            "%d copies had fragments out of ascending order; %s",
            unsorted_copies,
            "sorted them" if sort_blocks else "left them as given (sort_blocks=False)",
        )
    if mixed_chrom:
        logger.warning(
            "%d copies span more than one chromosome; only the first is used", mixed_chrom
        )
    if mixed_strand:
        logger.warning("%d copies carry more than one strand; only the first is used", mixed_strand)

    return meta_id, chrom_list, start_list, end_list, strand_list


def resolve_output_path(save_to_file, coordinate_table, default_name):
    """Work out where to write, from a path, or ``True`` plus a default name.

    ``save_to_file`` may be a path (used as given) or ``True`` (``default_name``
    next to the coordinate table, or in the current directory if the table was
    passed as a DataFrame).
    """
    if isinstance(save_to_file, str):
        return save_to_file
    if isinstance(coordinate_table, str):
        output_dir = os.path.dirname(os.path.abspath(coordinate_table))
    else:
        # Before 0.9.0 this fell back to the installed package directory, which
        # is read-only in a normal install and is not where anyone looks.
        output_dir = os.getcwd()
    return os.path.join(output_dir, default_name)


def save_or_return(payload, save_to_file, coordinate_table, default_name, label):
    """Pickle ``payload`` if asked to, otherwise return it.

    A truthy ``save_to_file`` saves; a string is treated as the destination
    path. Before 0.9.0 a string was compared against ``True`` and therefore
    silently skipped the save.
    """
    if not save_to_file:
        logger.info("done, returning %s", label)
        return payload

    import compress_pickle

    output_filepath = resolve_output_path(save_to_file, coordinate_table, default_name)
    # set_default_extension=False: compress_pickle otherwise replaces the
    # extension with the compression's own, so the old code reported
    # "saving at .../vcf_out.p" while actually writing vcf_out.lzma.
    compress_pickle.dump(
        payload, output_filepath, compression="lzma", set_default_extension=False
    )
    logger.info("done, saving %s at: %s", label, output_filepath)
    return output_filepath


def worker_count(env_var="MA_MAPPER_WORKERS", default=None):
    """Number of worker processes to use.

    Reads ``env_var``, then ``MA_MAPPER_WORKERS``, then ``default``, and
    finally the number of CPUs this process is actually allowed to use --
    ``os.sched_getaffinity`` under a Slurm allocation, not the machine's core
    count. Before 0.9.0 every extractor hardcoded 40 workers, which on a shared
    filesystem multiplies I/O rather than speeding anything up and on a 4-core
    allocation oversubscribes it tenfold.
    """
    for name in (env_var, "MA_MAPPER_WORKERS"):
        value = os.environ.get(name)
        if value:
            return max(1, int(value))
    if default is not None:
        return max(1, int(default))
    if hasattr(os, "sched_getaffinity"):
        return max(1, len(os.sched_getaffinity(0)))
    return max(1, os.cpu_count() or 1)
