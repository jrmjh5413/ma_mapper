#%%
import numpy as np
import pandas as pd
import os
import sys
from concurrent.futures import ProcessPoolExecutor
from itertools import repeat 
from Bio.AlignIO import MafIO
from . import logger
if sys.version_info >= (3, 8, 0):
    from typing import Literal, Tuple, List
else:
    from typing_extensions import Literal, Tuple, List

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import Counter, defaultdict
import time
from bisect import bisect_left, bisect_right
def get_spliced_mod(self, starts, ends, strand=1):
    # Dictionary for IUPAC ambiguity codes for 2-base combinations
    iupac_code = {
        frozenset(['A', 'G']): 'R', frozenset(['C', 'T']): 'Y',
        frozenset(['G', 'C']): 'S', frozenset(['A', 'T']): 'W',
        frozenset(['G', 'T']): 'K', frozenset(['A', 'C']): 'M',
        frozenset(['C', 'G', 'T']): 'B', frozenset(['A', 'G', 'T']): 'D',
        frozenset(['A', 'C', 'T']): 'H', frozenset(['A', 'C', 'G']): 'V',
        frozenset(['A', 'C', 'G', 'T']): 'N'
    }
    from typing import List, Tuple
    def trim_overlapping_fragments(starts: List[int], ends: List[int]) -> Tuple[List[int], List[int]]:
        """
        Trim overlapping TE fragments by adjusting the start positions to avoid base duplication.

        This is used to **prevent per-base duplication** when processing overlapping TE fragments,
        especially in multi-species alignments. It ensures that each base in the reference genome
        is only used once, even if it is part of multiple fragment calls.

        This function is necessary for **multi-species sequence alignments** where 
            overlapping TE fragments (e.g. from RepeatMasker) can introduce duplicated sequence 
            content if not properly merged. Without this step, fragments with physical overlap 
            in one species can be duplicated in the alignment, causing:

                - False inflation of sequence length
                - Artificially high TE conservation across species
                - Misleading phylogenetic signals
                - Incorrect ancestral reconstruction

            In contrast, if working in **single-genome**, retaining overlaps 
            may preserve internal deletions or tandem repeat events and should be handled differently.

        Parameters:
        -----------
        starts : List[int]
            List of start positions.
        ends : List[int]
            List of end positions (must be same length as `starts`).

        Returns:
        --------
        Tuple[List[int], List[int]]
            Adjusted (trimmed) start and end positions with no overlapping bases.
        """
        if not starts or not ends or len(starts) != len(ends):
            raise ValueError("Starts and ends must be non-empty and of equal length.")

        trimmed_starts = [starts[0]]
        trimmed_ends = [ends[0]]

        for i in range(1, len(starts)):
            prev_end = trimmed_ends[-1]
            curr_start = starts[i]
            curr_end = ends[i]

            # Trim current start to be >= previous end
            adjusted_start = max(curr_start, prev_end)
            if adjusted_start >= curr_end:
                # Skip completely overlapped region
                continue
            trimmed_starts.append(adjusted_start)
            trimmed_ends.append(curr_end)

        return trimmed_starts, trimmed_ends

    def convert_to_iupac(sequence):
        unique_bases = frozenset(sequence)
        if len(unique_bases) == 1:
            return sequence[0].upper()  
        return iupac_code.get(unique_bases, 'N')  # Default to 'N' for any unhandled cases
    
    def process_sequence_localized(sequence):
        sequence = sequence.upper()
        filtered_sequence = [base for base in sequence if base != '-']

        #base_counts = Counter(sequence)
        #most_common_bases = base_counts.most_common()
        #max_count = most_common_bases[0][1]
        #consensus_bases = [base for base, count in most_common_bases if count == max_count]
        new_base = convert_to_iupac(filtered_sequence)
        return new_base

    if strand not in (1, -1): 
        raise ValueError("Strand must be 1 or -1, got %s" % strand)
    
    starts, ends = trim_overlapping_fragments(starts, ends)
    fetched = list(self.search(starts, ends))
    #return fetched
    expected_letters = sum(end - start for start, end in zip(starts, ends))
    if len(fetched) == 0:
        return pd.DataFrame({'seqid': [self._target_seqname], 'seq': [Seq("N" * expected_letters)]})

    all_seqnames = {sequence.id for multiseq in fetched for sequence in multiseq}
    split_by_position = {seq_name: {} for seq_name in all_seqnames}
    
    exon_positions = []
    for start, end in zip(starts, ends):
        
        exon_positions.extend(range(start, end))
    exon_positions = sorted(set(exon_positions))  # Ensure unique & sorted
    
    total_rec_length = 0
    ref_first_strand = None
    split_by_position = defaultdict(lambda: defaultdict(str))
    for multiseq in fetched:
        if len(multiseq) == 1 and multiseq[0].id == self._target_seqname:
            # Fast-path for human-only alignment block
            rec = multiseq[0]
            rec_start = rec.annotations["start"]
            gapped_seq = rec.seq
            if ref_first_strand is None:
                ref_first_strand = rec.annotations["strand"]
            if "-" not in gapped_seq:
                #print('debug-humanonly: superfast route')
                # Fully ungapped, directly assign bases by offset
                for exonstart, exonend in zip(starts, ends):
                    # Clamp the range to only what's within the alignment block
                    effective_start = max(exonstart, rec_start)
                    effective_end = min(exonend, rec_start + len(gapped_seq))
           
                    if effective_start >= effective_end:
                        # Exon lies completely outside this block — skip
                        continue

                    start_idx = effective_start - rec_start
                    end_idx = effective_end - rec_start
                    for offset, base in zip(range(effective_start, effective_end), gapped_seq[start_idx:end_idx]):
                        split_by_position[rec.id][offset] = base
            else:
                #print('debug-humanonly: fast route')
                # Gapped version: build position map efficiently
                valid_indexes = []
                positions = []
                realpos = rec_start
                for i, base in enumerate(gapped_seq):
                    if base != "-":
                        valid_indexes.append((i, realpos))
                        positions.append(realpos)
                        realpos += 1

                for exonstart, exonend in zip(starts, ends):
                    i_start = bisect_left(positions, exonstart)
                    i_end = bisect_right(positions, exonend - 1)
                 
                    for i in range(i_start, i_end):
                        aln_idx, real_pos = valid_indexes[i]
                    
                        split_by_position[rec.id][real_pos] += gapped_seq[aln_idx]

            continue  # Skip rest of loop for this block
        # Find the human reference sequence
        ref = next((r for r in multiseq if r.id == self._target_seqname), None)
        if ref is None:
            raise ValueError(
                "Did not find %s in alignment bundle" % (self._target_seqname,)
            )
        
        try:
            if ref_first_strand is None:
                ref_first_strand = ref.annotations["strand"]
                if ref_first_strand not in (1, -1):
                    raise ValueError("Strand must be 1 or -1")
            elif ref_first_strand != ref.annotations["strand"]:
                raise ValueError(
                    "Encountered strand='%s' on target seqname, expected '%s'"
                    % (strand, ref_first_strand)
                )
        except KeyError:
            raise ValueError(
                "No strand information for target seqname (%s)" % self._target_seqname
            ) from None

        ref_start = ref.annotations["start"]
        ref_seq = ref.seq

        # Gapped path: build ref_pos_map
        #print("debug-humanref: fast route")
        ref_pos_map = []
        realpos = ref_start
        for i, base in enumerate(ref_seq):
            if base != "-":
                ref_pos_map.append((i, realpos))
                realpos += 1
        
        block_positions  = [pos for _, pos in ref_pos_map]
     
        # Now extract aligned fragments for all sequences using the same alignment coordinates
        for exonstart, exonend in zip(starts, ends):
            i_start = bisect_left(block_positions, exonstart)
            i_end = bisect_right(block_positions, exonend - 1)
   
            if i_start == i_end:
                continue  # This exon not represented in the alignment block

            aln_indexes = [ref_pos_map[i][0] for i in range(i_start, i_end)]
            ref_coords = block_positions[i_start:i_end]  # genomic coords covered in this exon by this block
      
            for rec in multiseq:
                rec_id = rec.id
                gapped_seq = rec.seq
                frag = "".join(gapped_seq[i] for i in aln_indexes)

                for pos, base in zip(ref_coords, frag):
                    split_by_position[rec_id][pos] += base
          
            # Count how many times each seqrec.id appears (detect duplicates)
            id_counts = Counter(seqrec.id for seqrec in multiseq)
            duplicate_ids = {seq_id for seq_id, count in id_counts.items() if count > 1}

            # Process duplicates only
            for seq_id in duplicate_ids:
                for pos, seq_str in split_by_position[seq_id].items():
                    # If sequence length > 1, means concatenation happened (duplication)
                    if len(seq_str) > 1:
                        new_seq = process_sequence_localized(seq_str)
                        split_by_position[seq_id][pos] = new_seq        
    realpos_to_len = {
        pos: len(gapped_fragment)
        for pos, gapped_fragment in split_by_position[self._target_seqname].items()
        if len(gapped_fragment) > 1
    }

    seqid_list = []
    seq_list = []
    
    for seqid in all_seqnames:
        seq_split = split_by_position[seqid]
        seq_splice = []
        filler_char = "N" if seqid == self._target_seqname else "-"
        append = seq_splice.append

        for exonstart, exonend in zip(starts, ends):
            for real_pos in range(exonstart, exonend):
                if real_pos in seq_split:
                    append(seq_split[real_pos])
                elif real_pos in realpos_to_len:
                    append(filler_char * realpos_to_len[real_pos])
                else:
                    append(filler_char)

        seqid_list.append(seqid)
        seq_list.append(Seq("".join(seq_splice))) 
    target_index = seqid_list.index(self._target_seqname)
    if len(seq_list[target_index].replace("-", "")) != expected_letters:
        raise ValueError(
            "coordinates %s , %s Returning %s letters for target seqname (%s), expected %s"
            % (
                starts,
                ends,
                len(seq_list[target_index].replace("-", "")),
                self._target_seqname,
                expected_letters,
            )
        )

    ref_subseq_len = len(seq_list[target_index])
    for seqid, seq in zip(seqid_list, seq_list):
        if len(seq) != ref_subseq_len:
            raise ValueError(
                "Returning length %s for %s, expected %s"
                % (len(seq), seqid, ref_subseq_len)
            )
  
    # Create a DataFrame
    df = pd.DataFrame({
        'seqid': seqid_list,
        'seq': [seq.reverse_complement() if strand != ref_first_strand else seq for seq in seq_list]
    })
    return df

MafIO.MafIndex.get_spliced = get_spliced_mod
#%%
_AGEARG = Literal[None,'calibrate']
_COUNTARG = Literal['ref_freq','coverage','common_freq', 'common_raw','common_freq_nogap','a','t','c','g','base_count','base_freq','base_freq_nogap','raw']
def extract_maf(name:str,
                maf_file:str, 
                chrom:str, 
                start:int, 
                end:int, 
                strand:str,
                target_species:str = 'Homo_sapiens', 
                count_arg:_COUNTARG = 'common_freq',
                species_list:list=None,
                e_value_df:pd.DataFrame = None, 
                internal_id_df:pd.DataFrame = None,
                ):
    iupac_codes = {
        'A': {'A'},
        'C': {'C'},
        'G': {'G'},
        'T': {'T'},
        'R': {'A', 'G'},
        'Y': {'C', 'T'},
        'S': {'G', 'C'},
        'W': {'A', 'T'},
        'K': {'G', 'T'},
        'M': {'A', 'C'},
        'B': {'C', 'G', 'T'},
        'D': {'A', 'G', 'T'},
        'H': {'A', 'C', 'T'},
        'V': {'A', 'C', 'G'},
        'N': {'A', 'C', 'G', 'T'},
        '-': {'-'}
    }

    def count_bases_with_ambiguity(sequences):
        base_counts = Counter()
        for sequence in sequences:
            for base in sequence:
                possible_bases = iupac_codes.get(base.upper(), {base.upper()})
                # Distribute counts among all possible bases
                for possible_base in possible_bases:
                    base_counts[possible_base] += 1 / len(possible_bases)
        
        return dict(base_counts)
    def remove_hg38_gap_columns(df, target_species):
        # Find hg38 sequence
        target_row = df[df["seqid"].str.startswith(target_species)]
        if target_row.empty:
            return df  # no hg38, return unchanged
        
        target_row = target_row.iloc[0]["seq"]
        keep_indices = [i for i, base in enumerate(target_row) if base != '-']
        
        def filter_seq(seq):
            return Seq(''.join(seq[i] for i in keep_indices))

        df = df.copy()
        df["seq"] = df["seq"].apply(filter_seq)
        return df
    
    maf_id = f'{target_species}.{chrom}'
    mafindex_filedir = '.'.join(str.split(maf_file, sep='.')[:-1])
    mafindex_filepath = f'{mafindex_filedir}.mafindex'
    index_maf = MafIO.MafIndex(mafindex_filepath, maf_file, maf_id) 
    n_strand = -1 if strand == '-' else 1
    results =index_maf.get_spliced(start,end,n_strand)
    collector = remove_hg38_gap_columns(results, target_species)
    if species_list is not None:
        pattern = '|'.join(species_list)  
        collector = results[results['seqid'].str.contains(pattern)]
    else:
        if internal_id_df is not None:
            _internal_id = internal_id_df[internal_id_df['meta_id'] == internal_id]['internal_id'].values[0]
        else:
            _internal_id = name
        if e_value_df is not None:
            e_value_df['seqid'] = e_value_df.species + '.' + e_value_df.chr_code
            e_value_internal_id=e_value_df[e_value_df['internal_id'] == _internal_id]
        
            if e_value_internal_id.shape[0] <= 1:
                collector = results[results.seqid.str.contains(target_species)]
            else:
                collector = results[results.seqid.isin(e_value_internal_id.seqid.values)]

    # Deduplicate by species name: the same organism can appear with different chromosome
    # identifiers across MAF blocks (e.g. Pan_troglodytes.chr1 and Pan_troglodytes.chr7).
    # Without deduplication, that organism contributes multiple rows, inflating coverage
    # and biasing all other per-position statistics.
    species_names = collector['seqid'].str.split('.').str[0]
    if species_names.duplicated().any():
        iupac_rev = {frozenset(v): k for k, v in iupac_codes.items() if k != '-'}
        merged_rows = []
        for species, group in collector.groupby(species_names, sort=False):
            if len(group) == 1:
                merged_rows.append({'seqid': group.iloc[0]['seqid'], 'seq': group.iloc[0]['seq']})
            else:
                seqs = [str(s) for s in group['seq'].tolist()]
                merged = []
                for bases in zip(*seqs):
                    non_gap = [b.upper() for b in bases if b != '-']
                    if not non_gap:
                        merged.append('-')
                    else:
                        unique = frozenset(non_gap)
                        if len(unique) == 1:
                            merged.append(next(iter(unique)))
                        else:
                            merged.append(iupac_rev.get(unique, 'N'))
                merged_rows.append({'seqid': group.iloc[0]['seqid'], 'seq': Seq(''.join(merged))})
        collector = pd.DataFrame(merged_rows)

    if count_arg in ['raw','raw_genome']:
        seqids = collector['seqid'].copy()
        if count_arg == 'raw_genome':
            seqids = seqids.str.split('.').str[0]
        sequence_length = len(collector.iloc[0]['seq'])
        result = []
        for i in range(sequence_length):
            bases = collector['seq'].apply(lambda x: x[i])
            entry = '|'.join(f'{s}:{b}' for s, b in zip(seqids, bases))
            result.append(entry)
        return result

    try:
        ref_alleles = np.char.upper(collector[collector.seqid.str.contains(target_species)]['seq'].to_list())[0]
    except IndexError:
        print(f'IndexError:{name}\t{maf_file}\t{maf_id}\t{chrom}\t{start}\t{end}\t{strand}')
       
    
    array_transposed=np.array(collector['seq'].to_list()).transpose()
    output_array=[]
    for idx, pos_array in enumerate(array_transposed):
        frequencies =count_bases_with_ambiguity(np.char.upper(pos_array))
        ref_allele = str(ref_alleles[idx])
        total = sum(frequencies.values())
        if count_arg == 'ref_freq':
            alt_count = total - frequencies.get(ref_allele, 0)
            alt_freq=alt_count/total
            if total == 1:
                alt_freq = np.nan
            output_array.append(alt_freq)
        elif count_arg == 'coverage':
            output_array.append(total)
        elif count_arg in ['common_freq', 'common_raw']:
            common_allele = max(frequencies, key=frequencies.get)
            common_count = frequencies.get(common_allele)
            if count_arg == 'common_freq':
                common_freq = common_count/total
                output_array.append(common_freq)
            else:
                output_array.append(common_count)
        elif count_arg == 'common_freq_nogap':
            frequencies.pop('-', None)
            try:
                common_allele = max(frequencies, key=frequencies.get)
            except ValueError:
                print(f'ValueError:{name}\t{maf_file}\t{maf_id}\t{chrom}\t{start}\t{end}\t{strand}')
            common_count = frequencies.get(common_allele)
            common_freq = common_count/total
            output_array.append(common_freq)
        elif count_arg == 'base_count':
            output_array.append(frequencies)
        elif count_arg in ['base_freq','base_freq_nogap']:
            if count_arg == 'base_freq_nogap':
                frequencies.pop('-', None)
            
            count_freq = {base: count / total for base, count in frequencies.items()}
            output_array.append(count_freq)
    if isinstance(output_array, str):
        print('debug',name,output_array)
    return output_array
#%%
def get_maf_filepath(maf_dir, chrom):
    files = os.listdir(maf_dir)
    maf_filename = f"{chrom}.maf" 
    maf_files = [f for f in files if f.endswith(maf_filename)]

    # Determine the appropriate file to use
    maf_filepath = None
    maf_match = next((f for f in maf_files if f.endswith('.maf') and not f.endswith('.maf.gz')), None)
    gz_match = next((f for f in maf_files if f.endswith('.maf.gz')), None)
    if maf_match:
        maf_filepath = f"{maf_dir}/{maf_match}"
    elif gz_match:
        maf_filepath = f"{maf_dir}/{gz_match}"
    else:
        raise FileNotFoundError(f"No .maf or .maf.gz file found for chromosome {chrom} in {maf_dir}")

    return maf_filepath
#%%
def maf_io(coordinate_table: pd.DataFrame|str, 
           maf:str,
           separated_maf:bool = False, 
           target_species:str = 'Homo_sapiens', 
           count_arg = 'common_freq', 
           save_to_file:bool = False, 
           generate_new_id:bool = False, 
           species_list:list=None,
           e_value_table: str|pd.DataFrame|None=None, 
           internal_id_table: str|pd.DataFrame|None=None, **kwargs):
    if isinstance(coordinate_table, str):
        if os.path.isfile(coordinate_table):
            coordinate_local = pd.read_csv(coordinate_table, sep='\t', header=None)
        else:
            print('coordinate_table file not found')
    else:
        coordinate_local = coordinate_table
    

    print(f'extract from maf target: {maf}')
    if generate_new_id == True:
        meta_id = [f'entry_{index}' for index in coordinate_local.index.astype(str)]
        coordinate_local['meta_id'] = meta_id
    else:
        coordinate_local['meta_id'] = coordinate_local.iloc[:,3]
        meta_id = coordinate_local.meta_id.unique()

    grouped = coordinate_local.groupby('meta_id', sort=False)
    chrom_list = grouped.apply(lambda x: x.iloc[:,0].unique()[0]).tolist()
    start_list = grouped.apply(lambda x: x.iloc[:,1].tolist()).tolist()
    end_list = grouped.apply(lambda x: x.iloc[:,2].tolist()).tolist()
    strand_list = grouped.apply(lambda x: x.iloc[:,5].unique()[0]).tolist()
    maf_call_list = []
    for chrom in chrom_list:
        if separated_maf == True:
            maf_file = get_maf_filepath(maf, chrom)
        else:
            maf_file = maf
        maf_call_list.append(maf_file)
    e_value_df = None
    if isinstance(e_value_table, str):
        e_value_df = pd.read_csv(e_value_table, sep='\t')
    else:
        e_value_df = e_value_table
    internal_id_df = None
    if isinstance(internal_id_table, str):
        internal_id_df = pd.read_csv(internal_id_table, sep='\t')
    else:
        internal_id_df = internal_id_table
    
    with ProcessPoolExecutor(max_workers=40) as executor:
        results = executor.map(extract_maf, meta_id, maf_call_list, chrom_list, start_list, end_list, strand_list, repeat(target_species), repeat(count_arg), repeat(e_value_df), repeat(internal_id_df), repeat(species_list))

    maf_out = []
    for result in results:
        maf_out.append(result)
    if save_to_file == True:
        if isinstance(save_to_file, str):
            output_filepath = save_to_file
        else:
            if isinstance(coordinate_table, str):
                output_dir = '/'.join(str.split(coordinate_table, sep ='/')[:-1])
            else: 
                output_dir = os.path.dirname(os.path.abspath(__file__))
            output_filepath = f'{output_dir}/maf_output.p'
        import compress_pickle
        compress_pickle.dump(maf_out, output_filepath, compression="lzma")
        print('done, saving maf_out at: ', output_filepath)
    else:
        print('done, returning maf_out as object')
        return maf_out
#%%
