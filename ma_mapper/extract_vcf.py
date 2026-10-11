"""Read a per-base value out of a VCF for each copy of a family.

For every fragment of every copy the extractor builds an array as long as the
fragment and adds the requested INFO field of each SNP to the base it sits on.
Minus-strand copies come back in the element's own 5'->3' orientation, which is
what ``mapper.map_data`` needs in order to put a value in the alignment column
of its own base.
"""

import os
from concurrent.futures import ProcessPoolExecutor
from itertools import repeat

import cyvcf2
import numpy as np

from . import logger
from ._io_common import group_coordinates, read_coordinate_table, save_or_return, worker_count

# One cyvcf2 handle per (worker process, vcf_file). The previous code opened the VCF
# inside the per-window loop, so a family with N copies performed N opens of a file that
# is 182 GB for gnomAD chr1 -- each open re-reads the tabix index and the (very large)
# gnomAD header. Multiplied by the worker pool and by concurrent jobs this saturated the
# rds-d6 filesystem (CSD3 support, 2026-08-02). Opening once per process is equivalent
# in output and reduces opens by 3-4 orders of magnitude.
_VCF_CACHE = {}


def _get_vcf(vcf_file):
    vcf = _VCF_CACHE.get(vcf_file)
    if vcf is None:
        vcf = _VCF_CACHE[vcf_file] = cyvcf2.VCF(vcf_file)
    return vcf


def _info_value(variant, query_key):
    """The INFO field as one number.

    A multiallelic record carries one value per ALT, which cyvcf2 returns as a
    tuple; adding that to a float array raised TypeError before 0.9.0. The
    values are summed, which for allele frequencies is the frequency of being
    non-reference at that base.
    """
    try:
        value = variant.INFO[query_key]
    except KeyError:
        return 0.0
    if isinstance(value, (tuple, list)):
        return float(sum(v for v in value if v is not None))
    if value is None:
        return 0.0
    return float(value)


def extract_vcf(vcf_file, chrom, start_list, end_list, strand, query_key):
    """Per-base INFO values for one copy, in element orientation.

    Parameters
    ----------
    vcf_file : str
        Path to a tabix-indexed VCF.
    chrom : str
        Chromosome name as spelled in the VCF.
    start_list, end_list : list of int
        Fragment coordinates of this copy, 0-based half-open (BED), in
        ascending order.
    strand : {'+', '-'}
        ``'-'`` reverses each fragment and the fragment order, so that index 0
        is the first base of the element rather than the leftmost genomic base.
    query_key : str
        INFO field to read, e.g. ``'AF'``. Bases with no SNP stay 0.0, as do
        SNPs whose record lacks the field.

    Returns
    -------
    numpy.ndarray
        float64, length ``sum(end - start)``.
    """
    try:
        vcf = _get_vcf(vcf_file)
    except OSError:
        logger.error("cannot open vcf %s", vcf_file)
        raise

    windows = []
    for i in range(len(start_list)):
        # zero-based to one-based conversion [BED -> VCF]
        start = start_list[i] + 1
        end = end_list[i]
        genome_position = f"{chrom}:{start}-{end}"
        window = np.zeros(end - start + 1, dtype=float)  # offset the +1 at the start
        for variant in vcf(genome_position):
            if variant.var_type != "snp":
                continue
            relative_position = variant.POS - start
            # A region query returns anything overlapping the window, so a
            # record can start before it; a negative index would otherwise
            # write silently to the far end of the array.
            if relative_position < 0 or relative_position >= window.size:
                continue
            window[relative_position] += _info_value(variant, query_key)
        # Flip ONCE, after the window is fully populated. This used to sit inside the
        # variant loop, which flipped the array once per variant: an even variant count
        # left the window unflipped entirely, an odd count left it flipped but with every
        # variant after the first written at a mirrored position. Minus-strand rows were
        # therefore anti-registered against their own alignment columns (~50% of copies).
        if strand == "-":
            window = np.flip(window)
        windows.append(window)

    if strand == "-":
        windows.reverse()
    return np.concatenate(windows).astype(float)


def vcf_io(
    coordinate_table,
    vcf,
    query_key="AF",
    vcf_format="gnomad",
    save_to_file=False,
    generate_new_id=False,
):
    """Extract a VCF INFO field for every copy in a coordinate table.

    Parameters
    ----------
    coordinate_table : pandas.DataFrame or str
        BED-like table, one or more rows per copy.
    vcf : str
        A single tabix-indexed VCF, or -- with ``vcf_format='gnomad'`` -- a
        directory holding one file per chromosome named
        ``gnomad.genomes.v3.1.1.sites.<chrom>.vcf.gz``.
    query_key : str
        INFO field to read.
    vcf_format : {'gnomad', 'single'}
        ``'gnomad'`` builds the per-chromosome filename shown above; anything
        else treats ``vcf`` as one file.
    save_to_file : bool or str
        ``True`` writes ``vcf_out.p`` beside the coordinate table; a string is
        used as the output path.
    generate_new_id : bool
        Treat every row as its own copy instead of grouping by column 3.

    Returns
    -------
    list of numpy.ndarray
        One array per copy, or the output path if the result was saved.
    """
    coordinate_local = read_coordinate_table(coordinate_table)
    meta_id, chrom_list, start_list, end_list, strand_list = group_coordinates(
        coordinate_local, generate_new_id=generate_new_id
    )

    vcf_call_list = []
    for chrom in chrom_list:
        if vcf_format == "gnomad":
            vcf_file = os.path.join(vcf, f"gnomad.genomes.v3.1.1.sites.{chrom}.vcf.gz")
        else:
            vcf_file = vcf
        vcf_call_list.append(vcf_file)

    missing = sorted({path for path in vcf_call_list if not os.path.isfile(path)})
    if missing:
        raise FileNotFoundError(
            f"{len(missing)} VCF file(s) not found, first: {missing[0]}"
            + (" (vcf_format='gnomad' expects one file per chromosome)" if vcf_format == "gnomad" else "")
        )

    workers = worker_count("MA_MAPPER_VCF_WORKERS")
    logger.info("extracting %s from %d copies with %d workers", query_key, len(meta_id), workers)
    with ProcessPoolExecutor(max_workers=workers) as executor:
        results = executor.map(
            extract_vcf,
            vcf_call_list,
            chrom_list,
            start_list,
            end_list,
            strand_list,
            repeat(query_key),
        )
        vcf_out = list(results)

    return save_or_return(vcf_out, save_to_file, coordinate_table, "vcf_out.lzma", "vcf_out")
