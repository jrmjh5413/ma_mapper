#%%
import pandas as pd

from ma_mapper import extract_vcf
from ma_mapper.extract_bigwig import bigwig_io
from . import sequence_alignment
import numpy as np
import sys
import os
from . import logger
import inspect
from functools import lru_cache
#from typing import Literal <- python3.8+
if sys.version_info >= (3, 8, 0):
    from typing import Literal, Tuple, List, Callable
else:
    from typing_extensions import Literal, Tuple, List

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
    logger.info('parse alignment')
    if isinstance(alignment_file, str):
        if (os.path.isfile(alignment_file) == True):
            records = load_alignment_file(alignment_file)
        else:
            logger.error(f"Alignment file '{alignment_file}' not found.")
            raise FileNotFoundError(f"Alignment file '{alignment_file}' not found.")
    else:
        records = alignment_file
    #counting sequence records
    seq_count = 0
    for i, value in enumerate(records):
        if seq_count == 0:
            seq_length = len(value)
        seq_count = seq_count + 1
    if seq_count == 0:
        # An empty/malformed alignment used to fall through to the np.zeros call
        # below and surface as UnboundLocalError('seq_length'), which hid the real
        # cause. 15 high-copy families (AluJb, MIR, L2c ...) were silently skipped
        # in the 2026-08-08 subfam_matrices_000 build for exactly this reason:
        # default MAFFT is O(N^2) and had written 0-byte .fasta.aligned files.
        src = alignment_file if isinstance(alignment_file, str) else '<records>'
        raise ValueError(
            "alignment '%s' contains no records (empty or malformed file)" % src)
    else:
        seq_id_list = list()
        for idx, content in enumerate(records):
            seq_id_list.append(content.name)
    #create parsed array
    parsed_alignment = np.zeros((seq_count, seq_length), dtype=np.uint8)
    if isinstance(alignment_file, str):
        if (os.path.isfile(alignment_file) == True):
            records = load_alignment_file(alignment_file)
        else:
            logger.error('alignment file not found')
    else:
        records = alignment_file
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
    # Unmapped characters now silently become gaps; say so rather than hide it.
    known = np.zeros(256, dtype=bool)
    known[list(hash_table)] = True
    unmapped = ~known[parsed_alignment]
    if unmapped.any():
        chars = sorted({chr(c) for c in np.unique(parsed_alignment[unmapped])})
        logger.warning(
            f"{int(unmapped.sum())} alignment characters outside "
            f"{{-,.,A,C,T,G,N}} were treated as gaps: {chars}"
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
            hf.create_dataset("parsed_alignment",  data=alignment_matrix, compression="lzf")
        logger.info('done, saving parsed alignment at: ',output_file)
    else: 
        return alignment_matrix
#%%
def import_parsed_alignment(parsed_alignment_file):
    import h5py
    with h5py.File(f'{parsed_alignment_file}', 'r') as hf:
        parsed_alignment = hf['result_array'][:]  
    return parsed_alignment
#%%
def create_filter(parsed_alignment, col_threshold = 0.50, col_content_threshold = 0.10, row_threshold = 0.50, save_to_file = False):
    if isinstance(parsed_alignment, str) == True:
        parsed_alignment = import_parsed_alignment(parsed_alignment)
    col_counts = np.zeros(parsed_alignment.shape[1])
    col_counts_max = np.zeros(parsed_alignment.shape[1])  
    for i, row in enumerate(parsed_alignment):
        nonzeros = np.nonzero(row)[0]
        col_counts[nonzeros] += 1
        col_counts_max[nonzeros[0]:nonzeros[-1]] += 1
    col_filter = (col_counts > col_counts_max * col_threshold) & (col_counts >= parsed_alignment.shape[0] * col_content_threshold)
    row_ratio = np.zeros(parsed_alignment.shape[0])
    for i, row in enumerate(parsed_alignment[:, col_filter]):
        nonzeros = np.nonzero(row)[0]
        if nonzeros.size and float(nonzeros[-1] - nonzeros[0]) > 0:
            row_ratio[i] = len(nonzeros) / (float(nonzeros[-1] - nonzeros[0]))
            #print(row_ratio[i])
    row_filter = row_ratio >= row_threshold
    for i, row in  enumerate(parsed_alignment[row_filter, :]):
        nonzeros = np.nonzero(row)[0]

        col_counts[nonzeros] += 1
        col_counts_max[nonzeros[0]:nonzeros[-1]] += 1
    col_filter = (col_counts > col_counts_max * col_threshold) & (col_counts >= parsed_alignment.shape[0] * col_content_threshold)
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
def map_data(extracted_data:str,
             alignment_matrix: np.ndarray, 
             filter: bool|None|list = True, **kwargs)-> np.ndarray:
    if isinstance(extracted_data, str) == True:
        import compress_pickle
        extracted_data = compress_pickle.load(extracted_data)
    def get_unique_element_types(matrix):
        unique_types = set()
        for row in matrix:
            for element in row:
                unique_types.add(type(element))
        return unique_types
    

    unique_types = get_unique_element_types(extracted_data)
    logger.info(unique_types)
    if any(t in unique_types for t in (int,float,np.uint8,np.uint16, np.float32, np.float64)) :
        nested_data = False
    else:
        nested_data = True
    logger.info(f'nested_data:{nested_data}')
    if nested_data:
        canvas = np.full(alignment_matrix.shape, None)
    else:
        # float32, not float64: genome-wide tracks (phyloP, CADD, allele
        # frequencies, coverage) carry far fewer than 7 significant digits, and
        # this array is the single largest allocation in the pipeline.
        canvas = np.zeros(alignment_matrix.shape, dtype=np.float32)
        
    for idx, row in enumerate(alignment_matrix):
        #print(f"data_file[{idx}]: {data_file[idx]}")
        #print(f"canvas[{idx}]: {canvas[idx, row>0]}")
        canvas[idx, row>0] = extracted_data[idx]
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

    base_count = np.apply_along_axis(lambda col: np.bincount(col, minlength=5), axis=0, arr=alignment_matrix)
    return base_count

#FIXED: wrong nonzero
_NORM_METHOD = Literal['average','perc_coverage','median']
def normalise(alignment_matrix: str|pd.DataFrame,
              data_matrix: pd.DataFrame, 
              method:_NORM_METHOD = 'average', 
              ):
    if method == 'average':
        normalizer = np.count_nonzero(alignment_matrix, axis=0)
        data_sum = np.nansum(data_matrix, axis=0)
        normalised_data = data_sum/normalizer 
    if method == 'perc_coverage':
        normalizer = np.count_nonzero(alignment_matrix, axis=0)
        data_count = np.count_nonzero(data_matrix, axis=0)
        normalised_data = data_count/normalizer * 100
    if method == 'median':
        normalised_data = np.nanmedian(data_matrix, axis = 0)
    return normalised_data

def flank_sequence_io(coordinate_table: pd.DataFrame,
                    source_fasta: str,
                    coordinate_table_out:bool = False, 
                    extension_length:int = 500) -> Tuple[pd.DataFrame, pd.DataFrame]:

    front_coordinate_table = coordinate_table.reset_index(drop=True).copy()
    front_coordinate_table.end = front_coordinate_table.start.astype(int)
    front_coordinate_table.start = front_coordinate_table.start.astype(int) - extension_length
    
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
    if filter:
        if filter == 'one_pass':
            filters=create_filter(parsed_alignment, **filter_kwargs(create_filter, kwargs))
        if filter == 'iterative':
            filters=iterative_filter(parsed_alignment, **filter_kwargs(iterative_filter, kwargs))
        row_filter = filters[0]
        col_filter = filters[1]
        parsed_alignment_filtered=parsed_alignment[np.ix_(row_filter,col_filter)]
        coordinate_table_filtered=coordinate_table.iloc[row_filter,:]
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

_FORMAT = Literal['bam_max','bam_min','bam_forward','bam_reverse','bigwig','bed']
def map_and_overlay(alignment:str,
            data_file:str,
            data_format:_FORMAT,
            filter:bool = True,
            coordinate_file:str|None = None,
            extension_length: int|None = None,
            source_fasta:str = None,
            **kwargs,
            ):
    combined_kwargs = filter_kwargs([parse_and_filter, create_filter], kwargs)
    fs_kwargs = {key[3:]: value for key, value in kwargs.items() if key.startswith('fs_')}   
    md_kwargs = {key[3:]: value for key, value in kwargs.items() if key.startswith('md_')}
    kwargs = {key: value for key, value in kwargs.items() if not key.startswith(('pf_','md_','fs_'))}
    alignment_matrix, alignment_coordinate, filters  = parse_and_filter(alignment_file=alignment, preprocess_out=True, **combined_kwargs)
    #use alignment coordinate as coordinate file if there is no coordinate file input
    if coordinate_file is None:
        coordinate_table = alignment_coordinate
    elif isinstance(coordinate_file, str):
        if (os.path.isfile(coordinate_table) == True):
            coordinate_table = pd.read_csv(coordinate_file, sep='\t', header=None)
        else:
            logger.error('coordinate file not found')
    else:
        coordinate_table = coordinate_file   
    if data_format in ['read_max','read_min','read_forward','read_reverse', 'normal_max','normal_min','normal_forward','normal_reverse','read_sum']:
        from . import extract_bam
        bam_io_kwargs = filter_kwargs(extract_bam.bam_io, kwargs)
        extracted_data=extract_bam.bam_io(coordinate_table= coordinate_table, bam_file=data_file,bam_format=data_format,**bam_io_kwargs)
    elif data_format in ['bigwig']:
        from . import extract_bigwig
        bigwig_io_kwargs = filter_kwargs(extract_bigwig.bigwig_io, kwargs)
        extracted_data=extract_bigwig.bigwig_io(coordinate_table=coordinate_table, bigwig=data_file,**bigwig_io_kwargs)
    elif data_format in ['bed']:
        from . import extract_bed
        bed_io_kwargs = filter_kwargs(extract_bed.bed_io, kwargs)
        extracted_data=extract_bed.bed_io(coordinate_table=coordinate_table, bed=data_file,**bed_io_kwargs)
    elif data_format in ['maf']:
        from . import extract_maf
        maf_io_kwargs = filter_kwargs(extract_maf.maf_io, kwargs)
        extracted_data=extract_maf.maf_io(coordinate_table=coordinate_table, maf=data_file,**maf_io_kwargs)
    elif data_format in ['vcf']:
        vcf_io_kwargs = filter_kwargs(extract_vcf.vcf_io, kwargs)
        extracted_data=extract_vcf.vcf_io(coordinate_table=coordinate_table, vcf=data_file, **vcf_io_kwargs)
    else:
        logger.error('dataformat field unspecified')
    row_filter, col_filter = filters

    if coordinate_file is None:
        alignment_matrix_sorted = alignment_matrix
        alignment_coordinate_sorted = alignment_coordinate
        row_filter_sorted = row_filter
    else: #if coordinate table is specified, alignment will be rearranged by NAME field in the coordinate table
        source_order = coordinate_table.iloc[:,3].unique()
        sorted_order = []
        for _, row in alignment_coordinate.iterrows():
            sorted_order.append(np.where(source_order == row['name'])[0][0])
        alignment_matrix_sorted = alignment_matrix[sorted_order]
        alignment_coordinate_sorted = alignment_coordinate.iloc[sorted_order]
        row_filter_sorted = row_filter[sorted_order]
    if filter:
        if isinstance(filter, list):
            mapping_filters = filter
        else:
            mapping_filters = [row_filter_sorted, col_filter]
        alignment_coordinate_filtered=alignment_coordinate_sorted.iloc[row_filter_sorted]
    else:
        mapping_filters = None
        alignment_coordinate_filtered=alignment_coordinate_sorted
    output_matrix_filtered = map_data(extracted_data=extracted_data, alignment_matrix=alignment_matrix_sorted, filter=mapping_filters, **md_kwargs)
    
    
    if extension_length is not None:
        
        front_coordinate_table, back_coordinate_table = flank_sequence_io(alignment_coordinate_filtered, coordinate_table_out=True, extension_length=extension_length,source_fasta=source_fasta, **fs_kwargs)
        if data_format in ['read_max','read_min','read_forward','read_reverse', 'normal_max','normal_min','normal_forward','normal_reverse','read_sum']:
            bam_io_kwargs = filter_kwargs(extract_bam.bam_io, kwargs)
            front_extracted_data=extract_bam.bam_io(coordinate_table= front_coordinate_table, bam_file=data_file,bam_format=data_format, **bam_io_kwargs)
            back_extracted_data=extract_bam.bam_io(coordinate_table= back_coordinate_table, bam_file=data_file,bam_format=data_format, **bam_io_kwargs)
        elif data_format in ['bigwig']:
            bigwig_io_kwargs = filter_kwargs(extract_bigwig.bigwig_io, kwargs)
            front_extracted_data=extract_bigwig.bigwig_io(coordinate_table=front_coordinate_table, bigwig=data_file,**bigwig_io_kwargs)
            back_extracted_data=extract_bigwig.bigwig_io(coordinate_table=back_coordinate_table, bigwig=data_file,**bigwig_io_kwargs)
        elif data_format in ['bed']:
            bed_io_kwargs = filter_kwargs(extract_bed.bed_io, kwargs)
            front_extracted_data=extract_bed.bed_io(coordinate_table=front_coordinate_table, bed=data_file,**bed_io_kwargs)
            back_extracted_data=extract_bed.bed_io(coordinate_table=back_coordinate_table, bed=data_file,**bed_io_kwargs)
        elif data_format in ['maf']:
            maf_io_kwargs = filter_kwargs(extract_maf.maf_io, kwargs)
            front_extracted_data=extract_maf.maf_io(coordinate_table=front_coordinate_table, maf=data_file,**maf_io_kwargs)
            back_extracted_data=extract_maf.maf_io(coordinate_table=back_coordinate_table, maf=data_file,**maf_io_kwargs)
        elif data_format in ['vcf']:
            vcf_io_kwargs = filter_kwargs(extract_vcf.vcf_io, kwargs)
            front_extracted_data=extract_vcf.vcf_io(coordinate_table=front_coordinate_table, vcf=data_file,**vcf_io_kwargs)
            back_extracted_data=extract_vcf.vcf_io(coordinate_table=back_coordinate_table, vcf=data_file,**vcf_io_kwargs)
  
        print(alignment_coordinate_sorted.shape, len(front_extracted_data), len(back_extracted_data))
        output_matrix_filtered = np.hstack((front_extracted_data,output_matrix_filtered,back_extracted_data))


    return output_matrix_filtered