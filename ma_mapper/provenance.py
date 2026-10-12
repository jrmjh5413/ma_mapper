"""Record how a matrix was made, next to the matrix.

A `.npy` on disk says nothing about which alignment it came from, which track,
which filter thresholds, or which version of this package placed the values.
Six months later that is the difference between a result and a file. Every
``ma-mapper overlay`` run writes a JSON record beside its output, and
:func:`write_run_record` does the same from Python.

Large inputs are not hashed: a checksum of a 182 GB gnomAD VCF would cost more
than the analysis. Files above ``MAX_HASH_BYTES`` are identified by size and
modification time instead, and the record says which of the two it used.
"""

import hashlib
import json
import os
import platform
import socket
import sys
import time

from . import __version__, logger

#: Inputs larger than this are identified by size and mtime rather than hashed.
MAX_HASH_BYTES = 256 * 1024 * 1024


def _sha256(path, chunk=1 << 20):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(chunk), b""):
            digest.update(block)
    return digest.hexdigest()


def describe_file(path, max_hash_bytes=MAX_HASH_BYTES):
    """Identify one file: path, size, mtime, and a checksum if it is small enough."""
    if path is None:
        return None
    record = {"path": os.path.abspath(path)}
    try:
        stat = os.stat(path)
    except OSError as exc:
        record["error"] = str(exc)
        return record
    record["bytes"] = stat.st_size
    record["mtime_utc"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime(stat.st_mtime))
    if stat.st_size <= max_hash_bytes:
        try:
            record["sha256"] = _sha256(path)
        except OSError as exc:  # unreadable, or vanished between stat and read
            record["error"] = str(exc)
    else:
        record["identified_by"] = "size and mtime (too large to hash)"
    return record


def _package_versions():
    import importlib

    versions = {"python": sys.version.split()[0], "ma_mapper": __version__}
    for name in ("numpy", "pandas", "Bio", "pysam", "pyBigWig", "pybedtools", "cyvcf2", "h5py"):
        try:
            versions[name] = getattr(importlib.import_module(name), "__version__", "unknown")
        except Exception:
            pass
    return versions


def run_record(inputs, parameters, outputs=None, command=None, extra=None):
    """Build the record. Returns a plain dict, ready to serialise.

    Parameters
    ----------
    inputs : dict
        ``{role: path}``, e.g. ``{'alignment': ..., 'data': ...}``. ``None``
        values are kept as ``None`` so an unused role is explicit rather than
        absent.
    parameters : dict
        Everything that would have to be repeated to get the same answer:
        the data format, the filter mode and its thresholds, and any extractor
        options that are not the defaults.
    outputs : dict, optional
        ``{role: {...}}``, typically the path, shape and dtype of each file.
    command : list or str, optional
        The command line, if there was one.
    extra : dict, optional
        Anything else worth keeping.
    """
    record = {
        "ma_mapper_version": __version__,
        "created_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "host": socket.gethostname(),
        "platform": platform.platform(),
        "working_directory": os.getcwd(),
        "command": command if isinstance(command, (str, type(None))) else " ".join(command),
        "inputs": {role: describe_file(path) for role, path in inputs.items()},
        "parameters": parameters,
        "outputs": outputs or {},
        "versions": _package_versions(),
    }
    slurm = {
        key: value for key, value in os.environ.items() if key.startswith("SLURM_JOB")
    }
    if slurm:
        record["slurm"] = slurm
    if extra:
        record.update(extra)
    return record


def write_run_record(path, inputs, parameters, outputs=None, command=None, extra=None):
    """Write :func:`run_record` to ``path`` as JSON. Returns the path."""
    record = run_record(inputs, parameters, outputs=outputs, command=command, extra=extra)
    with open(path, "w") as handle:
        json.dump(record, handle, indent=2, sort_keys=True, default=str)
        handle.write("\n")
    logger.info("wrote run record %s", path)
    return path


def read_run_record(path):
    """Read a record back."""
    with open(path) as handle:
        return json.load(handle)
