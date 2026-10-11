"""Helpers for the tables that surround an overlay: RepeatMasker and species ages."""

import numpy as np
import pandas as pd

from . import logger

#: Chromosomes kept by :func:`repeatmasker_prep` unless told otherwise.
HUMAN_MAIN_CHROMOSOMES = [
    "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9",
    "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17",
    "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY",
]

#: Repeat classes dropped by :func:`repeatmasker_prep` unless told otherwise:
#: simple repeats, satellites, low complexity, structural RNAs, and the
#: question-marked (uncertain) LTR and DNA classes.
DEFAULT_DISCARD_CLASS = [
    "Simple_repeat", "Satellite/telo", "Low_complexity", "snRNA", "tRNA",
    "Satellite", "srpRNA", "rRNA", "scRNA", "RNA", "Satellite/centr",
    "Satellite/acro", "LTR/ERVL?", "LTR/ERV1?", "LTR?", "LTR/Gypsy?",
    "DNA/hAT?", "RC?/Helitron?", "DNA?/hAT-Tip100?", "DNA?/PiggyBac?",
    "SINE?/tRNA",
]


def repeatmasker_prep(repeatmasker_filepath, main_chr=None, discard_class=None):
    """Load a RepeatMasker table and drop scaffolds and non-TE repeat classes.

    Parameters
    ----------
    repeatmasker_filepath : str or pandas.DataFrame
        UCSC-style ``rmsk`` table with ``genoName`` and ``repClass`` columns.
    main_chr : list or None
        Chromosomes to keep. The default is the human set in
        ``HUMAN_MAIN_CHROMOSOMES``; **pass this for any other assembly**,
        since mouse ``chr1..chr19`` survives the human list only by accident
        and ``chrM``, scaffolds and alternative haplotypes are always dropped.
        ``[]`` keeps every chromosome.
    discard_class : list or None
        Repeat classes to drop; defaults to ``DEFAULT_DISCARD_CLASS``.
        ``[]`` keeps every class.

    Returns
    -------
    pandas.DataFrame
        The filtered table.
    """
    if isinstance(repeatmasker_filepath, pd.DataFrame):
        rmskout_table = repeatmasker_filepath
    else:
        rmskout_table = pd.read_csv(repeatmasker_filepath, sep="\t", index_col=0, header=0)

    main_chr = HUMAN_MAIN_CHROMOSOMES if main_chr is None else main_chr
    discard_class = DEFAULT_DISCARD_CLASS if discard_class is None else discard_class

    keep = pd.Series(True, index=rmskout_table.index)
    if discard_class:
        keep &= ~rmskout_table["repClass"].isin(discard_class)
    if main_chr:
        keep &= rmskout_table["genoName"].isin(main_chr)

    filtered_table = rmskout_table[keep]
    logger.info(
        "repeatmasker_prep kept %d of %d rows (%d chromosomes requested)",
        len(filtered_table),
        len(rmskout_table),
        len(main_chr) if main_chr else rmskout_table["genoName"].nunique(),
    )
    return filtered_table


def divergence_table_prep(divergence_table_filepath, drop_unaligned=True):
    """Load a species divergence table and add an underscored ``meta_name``.

    Parameters
    ----------
    divergence_table_filepath : str
        Tab-separated table with ``Species`` and ``Source`` columns.
    drop_unaligned : bool
        Drop the species marked ``3. Zoonomia (not in alignment)``. Before
        0.9.0 this filter was computed and then thrown away: the function
        returned the unfiltered table, so those species were never actually
        excluded.
    """
    species = pd.read_table(divergence_table_filepath, sep="\t")
    species["meta_name"] = species["Species"].str.replace(" ", "_")
    if drop_unaligned:
        return species[species["Source"] != "3. Zoonomia (not in alignment)"]
    return species


#: Divergence times, in millions of years, of the reference points used by the
#: age tables below. The last entry is nan for elements with no age call.
AGE_CANON = [0, 6.7, 9.06, 15.76, 20.19, 29.44, 43.2, 73.8, 76, 82, 90, 96, 105, np.nan]


def age_canon():
    """The divergence times in ``AGE_CANON``."""
    return list(AGE_CANON)


def age_reference_table():
    """Divergence times with a representative taxon for each."""
    return {
        "age": list(AGE_CANON),
        "representative": [
            "human",
            "chimpanzee",
            "gorilla",
            "orangutan",
            "gibbon",
            "old-world monkey",
            "new-world monkey",
            "primates",
            "close primate relatives",
            "early primate ancestor",
            "rodents",
            "modern eutheria",
            "placental",
            "undefined",
        ],
        "representative_age": [
            "human (0)",
            "chimpanzee (6.7)",
            "gorilla (9.06)",
            "orangutan (15.76)",
            "gibbon (20.19)",
            "old-world monkey (29.44)",
            "new-world monkey (43.2)",
            "primates (73.8)",
            "close primate relatives (76)",
            "early primate ancestor (82)",
            "rodents (90)",
            "modern eutheria (96)",
            "placental mammals (105+)",
            "Undefined (NA)",
        ],
        "note": [
            "human",
            "ape",
            "ape",
            "ape",
            "gibbon",
            "old-world monkey",
            "new-world monkey",
            "primates",
            "colugos (flying Lemurs)",
            "primitive prosimian/primate",
            "rodents",
            "Boreoeutheria",
            "Atlantogenata",
            "Undefined (NA)",
        ],
    }


def age_reference_table_reorg():
    """``age_reference_table`` with the 82 and 90 My rows merged into their neighbours."""
    age_canon_reorg = [0, 6.7, 9.06, 15.76, 20.19, 29.44, 43.2, 73.8, 76, 96, 105, np.nan]
    return {
        "age": age_canon_reorg,
        "representative": [
            "human",
            "chimpanzee",
            "gorilla",
            "orangutan",
            "gibbon",
            "old-world monkey",
            "new-world monkey",
            "primates",
            "close primate relatives",
            "modern eutheria",
            "placental",
            "undefined",
        ],
        "representative_age": [
            "human (0)",
            "chimpanzee (6.7)",
            "gorilla (9.06)",
            "orangutan (15.76)",
            "gibbon (20.19)",
            "old-world monkey (29.44)",
            "new-world monkey (43.2)",
            "primates (73.8)",
            "close primate relatives (76)",
            "modern eutheria (96)",
            "placental mammals (105+)",
            "undefined (NA)",
        ],
        "note": [
            "human",
            "ape",
            "ape",
            "ape",
            "gibbon",
            "old-world monkey",
            "new-world monkey",
            "primates",
            "colugos (flying Lemurs)/primitive prosimian/primate",
            "rodents/Boreoeutheria",
            "Atlantogenata",
            "Undefinied",
        ],
    }
