"""One small synthetic family, and a dataset in every supported format.

Everything is built from scratch at test time, so the repository carries no
binary fixtures and the expected answers are arithmetic rather than recorded
output.

The toy family has four copies of a short element on one chromosome:

======  ====================  ======  ========================================
copy    coordinates (BED)     strand  shape
======  ====================  ======  ========================================
copyA   chrT 100-120          +       one fragment, 20 bases
copyB   chrT 200-220          -       one fragment, 20 bases
copyC   chrT 300-310,320-330  +       two fragments, 10 bases each
copyD   chrT 50-65            -       one fragment, 15 bases
======  ====================  ======  ========================================

Every signal fixture is built so that the value at genomic position ``p`` is
``p`` itself. An extractor that returns the wrong orientation, the wrong
fragment order or the wrong offset therefore cannot accidentally look right:
the expected arrays are ``[100..119]``, ``[219..200]``, ``[300..309, 320..329]``
and ``[64..50]``.
"""

import os
import shutil

import numpy as np
import pandas as pd
import pytest


@pytest.fixture(autouse=True, scope="session")
def _cap_worker_processes():
    """Two worker processes for the whole suite, not one per available CPU.

    The VCF and MAF extractors size their pool from the CPUs the process may
    use. On a login node or a shared runner that is a lot of processes for
    fixtures holding four copies, and it makes the timings noisy.
    """
    os.environ.setdefault("MA_MAPPER_WORKERS", "2")
    yield


def pytest_collection_modifyitems(config, items):
    """Skip the tests that need an external binary when it is not installed."""
    for marker, binary in (("requires_bedtools", "bedtools"), ("requires_mafft", "mafft")):
        if shutil.which(binary) is None:
            skip = pytest.mark.skip(reason=f"{binary} is not on PATH")
            for item in items:
                if marker in item.keywords:
                    item.add_marker(skip)

CHROM = "chrT"
CHROM_LENGTH = 400

COPIES = {
    "copyA": {"fragments": [(100, 120)], "strand": "+"},
    "copyB": {"fragments": [(200, 220)], "strand": "-"},
    "copyC": {"fragments": [(300, 310), (320, 330)], "strand": "+"},
    "copyD": {"fragments": [(50, 65)], "strand": "-"},
}

#: What a correct extractor returns for each copy, for a track whose value at
#: position p is p.
EXPECTED_POSITIONS = {
    "copyA": list(range(100, 120)),
    "copyB": list(range(219, 199, -1)),
    "copyC": list(range(300, 310)) + list(range(320, 330)),
    "copyD": list(range(64, 49, -1)),
}

# 21 alignment columns. Non-gap counts are 20, 20, 20 and 15, matching the
# fragment lengths above, which is the invariant map_data checks.
ALIGNMENT_ROWS = [
    ("copyA::chrT:100-120(+)", "ACGTACGTAC-GTACGTACGT"),
    ("copyB::chrT:200-220(-)", "ACGTACGTA-CGTACGTACGT"),
    ("copyC::chrT:300-310,320-330(+)", "ACGTACGTACG-TACGTACGT"),
    ("copyD::chrT:50-65(-)", "ACGTACG---TAC---GTACG"),
]


@pytest.fixture
def coordinate_table():
    """The copies as a six-column BED-like frame, one row per fragment."""
    rows = []
    for name, spec in COPIES.items():
        for start, end in spec["fragments"]:
            rows.append([CHROM, start, end, name, 10, spec["strand"]])
    return pd.DataFrame(rows)


@pytest.fixture
def genome_fasta(tmp_path):
    """A 400 bp chromosome with a repeating sequence, plus its .fai index."""
    import pysam

    sequence = ("ACGT" * (CHROM_LENGTH // 4))[:CHROM_LENGTH]
    path = tmp_path / "genome.fa"
    with open(path, "w") as handle:
        handle.write(f">{CHROM}\n")
        for i in range(0, len(sequence), 60):
            handle.write(sequence[i : i + 60] + "\n")
    pysam.faidx(str(path))
    return str(path)


@pytest.fixture
def alignment_fasta(tmp_path):
    """The toy alignment, with the coordinates in the headers."""
    path = tmp_path / "family.fasta.aligned"
    with open(path, "w") as handle:
        for name, sequence in ALIGNMENT_ROWS:
            handle.write(f">{name}\n{sequence}\n")
    return str(path)


@pytest.fixture
def position_bigwig(tmp_path):
    """A bigWig whose value at position p is p, with 250-260 left unmeasured."""
    import pyBigWig

    path = tmp_path / "positions.bw"
    bw = pyBigWig.open(str(path), "w")
    bw.addHeader([(CHROM, CHROM_LENGTH)])
    starts = [p for p in range(CHROM_LENGTH) if not 250 <= p < 260]
    bw.addEntries(
        [CHROM] * len(starts),
        starts,
        ends=[p + 1 for p in starts],
        values=[float(p) for p in starts],
    )
    bw.close()
    return str(path)


@pytest.fixture
def position_vcf(tmp_path):
    """SNPs at known positions, including three inside one minus-strand copy.

    copyB spans chrT 200-220 on the minus strand and contains SNPs at
    one-based 201, 203 and 207 -- three records in one window, which is what
    the pre-0.9.0 per-variant flip mangled. One record is multiallelic.
    """
    import pysam

    lines = [
        "##fileformat=VCFv4.2",
        f"##contig=<ID={CHROM},length={CHROM_LENGTH}>",
        '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
        f"{CHROM}\t101\t.\tA\tC\t.\tPASS\tAF=0.5",  # copyA, first base
        f"{CHROM}\t201\t.\tA\tC\t.\tPASS\tAF=0.1",
        f"{CHROM}\t203\t.\tA\tG\t.\tPASS\tAF=0.2",
        f"{CHROM}\t207\t.\tA\tT,G\t.\tPASS\tAF=0.25,0.15",  # multiallelic
        f"{CHROM}\t305\t.\tA\tG\t.\tPASS\tAF=0.3",
    ]
    plain = tmp_path / "variants.vcf"
    plain.write_text("\n".join(lines) + "\n")
    path = str(tmp_path / "variants.vcf.gz")
    pysam.tabix_compress(str(plain), path, force=True)
    pysam.tabix_index(path, preset="vcf", force=True)
    return path


@pytest.fixture
def feature_bed(tmp_path):
    """Features over copyA and copyB, with an integer, a float and a '.' score."""
    path = tmp_path / "features.bed"
    path.write_text(
        "\n".join(
            [
                f"{CHROM}\t100\t105\tfeat_int\t3\t+",
                f"{CHROM}\t110\t112\tfeat_float\t0.75\t+",
                f"{CHROM}\t205\t210\tfeat_dot\t.\t-",
                f"{CHROM}\t380\t390\tfeat_elsewhere\t9\t+",
            ]
        )
        + "\n"
    )
    return str(path)


@pytest.fixture
def read_bam(tmp_path):
    """A BAM with 300 forward reads starting at 105, 2 at 110, 1 reverse read.

    300 is deliberately above 255: the pre-0.9.0 uint8 profile wrapped round.
    The reverse read spans 105-115, so its 5' end (114) is inside copyA.
    """
    import pysam

    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": CHROM, "LN": CHROM_LENGTH}]}
    path = tmp_path / "reads.bam"
    with pysam.AlignmentFile(str(path), "wb", header=header) as out:
        def add(name, start, reverse):
            segment = pysam.AlignedSegment()
            segment.query_name = name
            segment.query_sequence = "ACGTACGTAC"
            segment.flag = 16 if reverse else 0
            segment.reference_id = 0
            segment.reference_start = start
            segment.mapping_quality = 60
            segment.cigarstring = "10M"
            segment.query_qualities = pysam.qualitystring_to_array("I" * 10)
            out.write(segment)

        for i in range(300):
            add(f"fwd{i}", 105, False)
        add("rev0", 105, True)  # 5' end at 114, inside copyA
        for i in range(2):
            add(f"fwd_b{i}", 110, False)
    pysam.index(str(path))
    return str(path)


@pytest.fixture
def three_species_maf(tmp_path):
    """A MAF over copyA with an insertion in one species and a duplicate row.

    ``Homo_sapiens.chrT`` is the reference. ``Mus_musculus`` carries one extra
    base, so the reference row holds a gap column that has to be dropped for
    the result to be 20 values long. ``Canis_lupus`` appears twice, under two
    chromosome names, which is what the species-level dedup is for.
    """
    ref = "ACGTACGTACGTACGTACGT"  # chrT:100-120
    lines = [
        "##maf version=1 scoring=none",
        "",
        "a score=0.0",
        f"s Homo_sapiens.{CHROM} 100 20 + {CHROM_LENGTH} {ref[:10]}-{ref[10:]}",
        f"s Mus_musculus.chr9 500 21 + 1000 {ref[:10]}T{ref[10:]}",
        f"s Canis_lupus.chr1 700 20 + 1000 {ref[:10]}-{ref[10:]}",
        f"s Canis_lupus.chr7 900 20 + 1000 {ref[:10]}-{ref[10:12]}AA{ref[14:]}",
        "",
    ]
    path = tmp_path / "toy.maf"
    path.write_text("\n".join(lines))
    return str(path)


def expected_array(copy_name):
    """The correct per-base result for ``copy_name`` on a position-valued track."""
    return np.array(EXPECTED_POSITIONS[copy_name], dtype=float)
