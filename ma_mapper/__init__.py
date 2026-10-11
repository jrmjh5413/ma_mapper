"""ma_mapper -- map genome-wide data onto a multiple alignment of transposable elements.

The package turns (a) a multiple sequence alignment of the copies of one TE
family and (b) a genome-wide dataset (BED, BAM, bigWig, VCF or MAF) into a
copy x alignment-column matrix, so that every value sits in the alignment
column of the base it was measured on.

Quick start
-----------
>>> from ma_mapper import mapper
>>> matrix = mapper.map_and_overlay(
...     alignment='MER11A.fasta.aligned',
...     data_file='phyloP447.bigWig',
...     data_format='bigwig',
... )

Logging
-------
Importing this package does not configure logging and does not create a file:
messages of level WARNING and above reach stderr through the standard library
fallback, and anything quieter is dropped. To see progress, call
``configure_logging()``, or set ``MA_MAPPER_LOG`` to a path before import.

Before 0.9.0 the import itself attached a file handler writing ``ma_mapper.log``
to the current working directory, which failed outright on a read-only
directory and littered every working directory it ran in.
"""

import logging
import os

__version__ = "0.9.0"

# Submodules are not imported here: mapper pulls in pysam, pyBigWig, cyvcf2 and
# pybedtools, and importing them all costs seconds even when the caller needs one.
__all__ = ["configure_logging", "logger", "__version__"]

logger = logging.getLogger("ma_mapper")
# A library attaches no handler of its own. Without this, Python prints
# "No handlers could be found" style complaints or, worse, a parent
# application's configuration is silently overridden.
logger.addHandler(logging.NullHandler())


def configure_logging(level=logging.INFO, logfile=None, stream=True):
    """Attach handlers to the ``ma_mapper`` logger. Safe to call more than once.

    Parameters
    ----------
    level : int or str
        Threshold for the handlers this call attaches.
    logfile : str or None
        Path to append log records to. ``None`` (the default) writes no file.
    stream : bool
        Also write records to stderr.

    Returns
    -------
    logging.Logger
        The package logger.
    """
    logger.setLevel(level)
    formatter = logging.Formatter("%(asctime)s %(name)s %(levelname)s %(message)s")

    existing = {(type(h).__name__, getattr(h, "baseFilename", None)) for h in logger.handlers}

    if logfile is not None:
        key = ("FileHandler", os.path.abspath(logfile))
        if key not in existing:
            file_handler = logging.FileHandler(logfile)
            file_handler.setLevel(level)
            file_handler.setFormatter(formatter)
            logger.addHandler(file_handler)

    if stream and ("StreamHandler", None) not in existing:
        stream_handler = logging.StreamHandler()
        stream_handler.setLevel(level)
        stream_handler.setFormatter(logging.Formatter("%(asctime)s %(levelname)s: %(message)s"))
        logger.addHandler(stream_handler)

    return logger


# Opt-in file logging for users who relied on the pre-0.9.0 behaviour.
if os.environ.get("MA_MAPPER_LOG"):
    configure_logging(logfile=os.environ["MA_MAPPER_LOG"], stream=False)
