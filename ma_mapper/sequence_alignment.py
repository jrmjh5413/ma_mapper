#%%
"""Pull genomic sequence for a set of copies and align it.

The alignment headers written here (``name::chrom:start-end(strand)``) are what
``mapper.extract_coordinate_from_alignment`` reads back, so the sequence and
its coordinates can never drift apart.
"""

import os
import subprocess

import pandas as pd

from . import logger
from ._io_common import read_coordinate_table

#%%
def extract_coord_from_repeatmasker_table(
        subfamily: str|list,
        repeatmasker_table: str|pd.DataFrame,
        internal_id_table: str|list|None=None,
        save_to_file:bool = False,
        ):
    
    if isinstance(subfamily, str):
        subfamilies = [subfamily]
    else:
        subfamilies = subfamily
    coord_df = []
    for idx,subfamily in enumerate(subfamilies):
        if isinstance(repeatmasker_table, str):
            repeatmasker_table = pd.read_csv(repeatmasker_table, sep='\t', index_col=0, header=0)

        if internal_id_table is not None:
            if isinstance(internal_id_table, list):
                internal_id_tbl = internal_id_table[idx]
            else:
                internal_id_tbl = internal_id_table
            if isinstance(internal_id_tbl, str):
                if os.path.isfile(internal_id_tbl):
                    internal_id_df = pd.read_csv(internal_id_tbl, sep='\t')
                else:
                    raise FileNotFoundError(f"Internal ID table {internal_id_tbl} not found")

            elif isinstance(internal_id_tbl, pd.DataFrame):
                internal_id_df = internal_id_tbl
            # Restrict to this subfamily first. Before 0.9.0 the merge used the
            # whole RepeatMasker table regardless of `subfamily`, so a list of N
            # subfamilies returned the same full table N times over.
            subfam_rows = repeatmasker_table[repeatmasker_table.repName == subfamily]
            subfam_w_internal_id=pd.merge(subfam_rows, internal_id_df, left_index = True, right_on = 'rmsk_index')
            subfam_coords = subfam_w_internal_id[['genoName','genoStart','genoEnd','internal_id']].copy()
            subfam_coords['score'] = 10
            subfam_coords['strand'] = subfam_w_internal_id.strand
        else:
            subfam_table=repeatmasker_table[repeatmasker_table.repName == subfamily].copy().reset_index()
            subfam_coords = subfam_table[['genoName','genoStart','genoEnd']].copy()
            subfam_coords['internal_id'] =subfam_table.index.map(lambda x: f'{subfamily}_{x}')
            subfam_coords['score'] = 10
            subfam_coords['strand'] = subfam_table.strand
        coord_df.append(subfam_coords)

    coords = pd.concat(coord_df).reset_index(drop=True)
    if save_to_file:
        if isinstance(save_to_file, str):
            output_filepath = save_to_file  
        else:
            output_filepath = f'{os.getcwd()}/{subfamily}.txt'
        coords.to_csv(output_filepath, sep='\t', index= False, header=False)
    return coords
#%%
class _IndexedGenome:
    """Random access to a genome FASTA through its .fai index.

    Reading hg38 with ``SeqIO.to_dict`` holds the whole assembly in memory as
    Python strings -- several gigabytes to fetch a few hundred base pairs,
    which is what every call did before 0.9.0.
    """

    def __init__(self, path):
        import pysam

        self._handle = pysam.FastaFile(path)
        self.references = set(self._handle.references)

    def fetch(self, chrom, start, end):
        return self._handle.fetch(chrom, start, end)

    def close(self):
        self._handle.close()

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()


class _InMemoryGenome:
    """Fallback for a FASTA that cannot be indexed (no .fai and no write access)."""

    def __init__(self, path):
        from Bio import SeqIO

        with open(path) as handle:
            self._records = SeqIO.to_dict(SeqIO.parse(handle, 'fasta'))
        self.references = set(self._records)

    def fetch(self, chrom, start, end):
        return str(self._records[chrom][start:end].seq)

    def close(self):
        self._records.clear()

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()


def open_genome(source_fasta):
    """Open a genome FASTA for random access, indexed if possible."""
    try:
        return _IndexedGenome(source_fasta)
    except Exception as exc:  # no .fai and no write access, or a non-FASTA file
        logger.warning(
            'could not use the FASTA index for %s (%s); reading the whole file '
            'into memory instead', source_fasta, exc
        )
        return _InMemoryGenome(source_fasta)


def extract_sequence(chrom, start, end, strand, records):
    """Genomic sequence for one interval, reverse complemented on the minus strand.

    ``records`` may be an object returned by :func:`open_genome` or, as before,
    a dict of Biopython records keyed by chromosome.
    """
    if hasattr(records, 'fetch'):
        seq_string = records.fetch(chrom, int(start), int(end))
    else:
        seq_string = str(records[chrom][int(start):int(end)].seq)
    if strand == '-':
        from Bio.Seq import Seq

        seq_string = str(Seq(seq_string).reverse_complement())
    return seq_string

def sequence_io(coordinate_table,
                source_fasta, 
                save_to_file=False, 
                generate_new_id = False, 
                ):
    """Fetch the genomic sequence of every copy as FASTA records.

    One record per copy, fragments concatenated in the element's own
    orientation, with the header
    ``name::chrom:start-end,start-end(strand)`` that
    ``mapper.extract_coordinate_from_alignment`` parses back.

    Parameters
    ----------
    coordinate_table : pandas.DataFrame or str
        BED-like table, one or more rows per copy.
    source_fasta : str
        Genome FASTA. Read through its .fai index when there is one.
    save_to_file : bool or str
        ``True`` writes ``sequence.fasta`` beside the coordinate table; a
        string is used as the path. Before 0.9.0 the ``True`` case built the
        path by appending ``.fasta`` to the table's *parent directory*.
    generate_new_id : bool
        Treat every row as its own copy instead of grouping by column 3.

    Returns
    -------
    list of SeqRecord
        Or the output path, if saved.
    """
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

    coordinate_local = read_coordinate_table(coordinate_table).copy()
    if generate_new_id:
        coordinate_local['meta_id'] = [f'entry_{index}' for index in coordinate_local.index.astype(str)]
    else:
        coordinate_local['meta_id'] = coordinate_local.iloc[:, 3]

    if hasattr(source_fasta, 'fetch'):
        genome, close_genome = source_fasta, False
    else:
        if not os.path.isfile(source_fasta):
            raise FileNotFoundError(f'source_fasta not found: {source_fasta}')
        genome, close_genome = open_genome(source_fasta), True

    seq_records = []
    try:
        missing_chroms = set()
        # Grouped in one pass. Before 0.9.0 the whole table was re-scanned for
        # every copy, which is quadratic: 100,000 copies meant 100,000 full
        # passes over 100,000 rows.
        for uniq_meta_id, group in coordinate_local.groupby('meta_id', sort=False):
            group = group.sort_values(group.columns[1])
            chrom = group.iloc[0, 0]
            strand = group.iloc[0, 5]
            if getattr(genome, 'references', None) and chrom not in genome.references:
                missing_chroms.add(chrom)
                continue
            seq_strings, start_list, end_list = [], [], []
            for row in group.itertuples(index=False, name=None):
                start, end = int(row[1]), int(row[2])
                seq_strings.append(extract_sequence(chrom, start, end, strand, genome))
                start_list.append(str(start))
                end_list.append(str(end))
            if strand == '-':
                seq_strings.reverse()
            coord_string = ",".join(f"{a}-{b}" for a, b in zip(start_list, end_list))
            seqname = f"{uniq_meta_id}::{chrom}:{coord_string}({strand})"
            seq_records.append(
                SeqRecord(Seq(''.join(seq_strings)), id=seqname, name='', description='')
            )
        if missing_chroms:
            logger.warning(
                'skipped %d copies on chromosomes absent from the FASTA: %s',
                len(missing_chroms), sorted(missing_chroms)[:5]
            )
    finally:
        if close_genome:
            genome.close()

    if save_to_file:
        if isinstance(save_to_file, str):
            output_filepath = save_to_file
        elif isinstance(coordinate_table, str):
            output_filepath = os.path.join(
                os.path.dirname(os.path.abspath(coordinate_table)), 'sequence.fasta'
            )
        else:
            output_filepath = os.path.join(os.getcwd(), 'sequence.fasta')
        with open(output_filepath, "w") as output_handle:
            SeqIO.write(seq_records, output_handle, "fasta")
        logger.info('wrote %d sequences to %s', len(seq_records), output_filepath)
        return output_filepath
    return seq_records


def parse_alignment_to_fasta(parsed_alignment,name=None, save_to_file=False):
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio import SeqIO

    # Define the mapping. 0 is a GAP, not an N: parse_alignment codes gap 0 and
    # N 5, so the old table wrote every gap out as an N and raised KeyError on
    # the first genuine N it met.
    base_mapping = {0: '-', 1: 'A', 2: 'C', 3: 'T', 4: 'G', 5: 'N'}

    if name is None:
        name = 'entry'
    # Convert each row to DNA sequence
    dna_sequences = []
    for idx, row in enumerate(parsed_alignment):
        seq_name = f'{name}_{idx}'
        dna_sequence = ''.join(base_mapping.get(int(val), 'N') for val in row)
        dna_sequences.append(SeqRecord(Seq(dna_sequence),seq_name , '', ''))

    # Write the sequences to a FASTA file
    if save_to_file:
        if isinstance(save_to_file, str):
            output_filepath = save_to_file  
        else:
            output_filepath = f'{os.getcwd()}/parsed_to_seq.fasta'
        with open(output_filepath, "w") as output_handle:
            SeqIO.write(dna_sequences, output_handle, "fasta")
    return dna_sequences

def mafft_align(input_filepath, nthread = None, nthreadtb = None, nthreadit =None, output_filepath = None, mafft_arg = 
                  #--quiet --memsave 
                  '--treeout --reorder'):
    """Align a FASTA with MAFFT's default algorithm (FFT-NS-2).

    Parameters
    ----------
    input_filepath : str
        Unaligned FASTA.
    nthread, nthreadtb, nthreadit : int or None
        MAFFT's ``--thread``, ``--threadtb`` and ``--threadit``. Before 0.9.0
        the last two were filled in with the value of ``nthread``, so setting
        them alone produced ``--threadb None``.
    output_filepath : str or None
        Defaults to ``<input_filepath>.aligned``.
    mafft_arg : str
        Extra flags, as one string.

    Notes
    -----
    Memory grows as about n^1.19 and runs out above roughly 65,000 copies; see
    ``docs/aligning_large_families.md`` for the profile-alignment route.

    The command is run without a shell, so paths containing spaces or shell
    metacharacters are passed through safely.
    """
    command = ['mafft']
    if nthread is not None:
        command += ['--thread', str(nthread)]
    if nthreadtb is not None:
        command += ['--threadtb', str(nthreadtb)]
    if nthreadit is not None:
        command += ['--threadit', str(nthreadit)]
    if mafft_arg:
        command += mafft_arg.split()
    command.append(input_filepath)
    if output_filepath is None:
        output_filepath = f'{input_filepath}.aligned' 
    logger.info('running: %s', ' '.join(command))
    try:
        completed = subprocess.run(command, capture_output=True, check=True)
    except FileNotFoundError:
        raise FileNotFoundError(
            'mafft was not found on PATH; install it (conda install -c bioconda mafft)'
        ) from None
    except subprocess.CalledProcessError as e:
        logger.error('MAFFT failed (exit %s): %s', e.returncode, e.stderr.decode(errors='replace')[-2000:])
        raise
    with open(output_filepath, 'wb') as file:
        file.write(completed.stdout)
    return output_filepath
#%%
