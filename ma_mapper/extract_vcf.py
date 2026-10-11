import numpy as np
import pandas as pd
import cyvcf2
import os
from concurrent.futures import ProcessPoolExecutor
from itertools import repeat
from . import logger
#%%
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


def extract_vcf(vcf_file, chrom, start_list, end_list, strand, query_key):
    windows = []
    try:
        vcf = _get_vcf(vcf_file)
    except OSError:
        logger.error("cannot open vcf ", vcf_file)
        raise
    for i in range(len(start_list)):
        # zero-based to one-based conversion [BED -> VCF]
        start = start_list[i] +1
        end = end_list[i]
        genome_position = f'{chrom}:{start}-{end}'
        window = np.zeros(end-start +1, dtype=float) #offset the +1 at the start
        for variant in vcf(genome_position):
            if variant.var_type == 'snp':
                relative_position = variant.POS-start
                try:
                    window[relative_position] +=  variant.INFO[query_key]
                except KeyError:
                    window[relative_position] +=  0.0
        # Flip ONCE, after the window is fully populated. This used to sit inside the
        # variant loop, which flipped the array once per variant: an even variant count
        # left the window unflipped entirely, an odd count left it flipped but with every
        # variant after the first written at a mirrored position. Minus-strand rows were
        # therefore anti-registered against their own alignment columns (~50% of copies).
        if strand == '-':
            window = np.flip(window)
        windows.append(window)
    if strand == '-':
        windows.reverse()
    window_out = np.concatenate(windows)
    return window_out.astype(float)
#%%
def vcf_io(coordinate_table, vcf, query_key = 'AF', vcf_format = 'gnomad', save_to_file = False, generate_new_id = False):
    if isinstance(coordinate_table, str):
        if (os.path.isfile(coordinate_table) == True):
            coordinate_local = pd.read_csv(coordinate_table, sep='\t', header=None)
        else:
            logger.error('coordinate_table file not found')
    else:
        coordinate_local = coordinate_table

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
    vcf_call_list = []
    for chrom in chrom_list:
        if vcf_format == 'gnomad':
            vcf_file = f'{vcf}/gnomad.genomes.v3.1.1.sites.{chrom}.vcf.gz'
        else:
            vcf_file = vcf
        vcf_call_list.append(vcf_file)
    # 40 was hardcoded; on a shared filesystem that is an I/O multiplier, not a speedup.
    # MA_MAPPER_VCF_WORKERS lets a caller (or a Slurm script) match it to allocated cores.
    _workers = int(os.environ.get('MA_MAPPER_VCF_WORKERS', '40'))
    with ProcessPoolExecutor(max_workers=_workers) as executor:
        results  = executor.map(extract_vcf, vcf_call_list, chrom_list, start_list,end_list,strand_list, repeat(query_key))
    vcf_out = []
    for result in results:
        vcf_out.append(result)
    if save_to_file == True:
        if isinstance(save_to_file, str):
            output_filepath = save_to_file
        else:
            if isinstance(coordinate_table, str):
                output_dir = '/'.join(str.split(coordinate_table, sep ='/')[:-1])
            else:
                output_dir = os.path.dirname(os.path.abspath(__file__))
            output_filepath = f'{output_dir}/vcf_out.p'
        import compress_pickle        
        compress_pickle.dump(vcf_out, output_filepath, compression="lzma")
        logger.info('done, saving vcf_out at: ', output_filepath)
    else:
        logger.info('done, returning vcf_out as object')
        return vcf_out
#%%