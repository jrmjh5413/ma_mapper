# ma_mapper.extract_maf
a Bio.AlignIO.MafIO wrapper for MAF file data extraction
## list of functions
- ma_mapper.extract_maf.maf_io()
- ma_mapper.extract_maf.extract_maf()
- ma_mapper.extract_maf.get_spliced_mod()

## ma_mapper.extract_maf.maf_io()

`ma_mapper.extract_maf.maf_io(coordinate_table, maf, separated_maf=False, target_species='Homo_sepiens, count_arg='common',save_to_file=False, generate_new_id=False, species_list=None, e_value_table=None, internal_id_table=None)`

Extract data from the input MAF file using the coodinates from the coordinate table. Noted that MAF file is actually multiple alingment file formats, the output of this package is expressed as base counts or frequencies.

### Parameters
**coordinate_table : str,pd.DataFrame, default: None**  
The coordinate table for data extraction target. Accept both filepath or pandas dataframe. The table has to be in BED file format.

**maf : str**  
The path to the MAF file or the directory of MAF files

**separated_maf : bool, default: False**  
Specify the type of MAF input. For some datasets, alignments are contained in a single MAF file, however, for a big dataset such as Zoonomia, separating alignment by chromosomes can improve speed.

**target_species : str, default: 'Homo_sapiens'**  
The species id used as a reference seqeunce for the alignment.

**count_arg : str, default: 'common'**  
The mode of base counting assigned to extract_maf().

- 'ref_freq': count reference base frequencies
- 'coverage': count total bases
- 'common_raw','common_freq','common_freq_nogap': count common base (raw count, freqeuncies, and frequencies without gaps)
- 'a','t','c','g': count ATCG
- 'base_count','base_freq','base_freq_nogap': count all bases (raw count, freqeuncies, and frequencies without gaps)
- 'raw': extract bases without counting 

**save_to_file : bool,str, default: False**  
Save the output as a pickled file with lzma compression. The filepath can be specified, using save_to_file = True would save the list to the same directory as the coordinate table or current working directory.

**generate_new_id : bool, default: False**  
Generate new id for using index of the coordinate table. This parameter dictate the fragment joining behavior. By default, if the coordinate table have more than one coordinates that share the same 'name' or id, the function would join the results from these fragments.

**species_list : list[str], default: None**  
The list of species in the alignments to be counted or extracted. The names need to be the same as the names in the alignments

**e_value_table : str, pandas.DataFrame, default: None**  
The path to the table of E-value calculated by teatime package or its pandas.DataFrame, works as the species list to filter for species with significant alignment to the reference.

**internal_id_table: str, pandas.DaraFrame, default: None**  
The path to the table of internal ID or its pandas.DataFrame, works as the reference linking the coordinates of TEs to internal ID used by teatime outputs.

#### Returns
**maf_out: list[numpy.ndarray]**  
The list of extracted data from MAF file

## ma_mapper.extract_maf.extract_maf()

`ma_mapper.extract_maf.extract_maf(name, maf_file, chrom, start, end, strand, target_species='Homo_sapiens', count_arg='common', species_list=None, e_value_df=None, internal_id_df=None)`

Extract data from the input MAF file from the specific genomic location. This is the wrapper for Bio.AlignIO.MafIO.

### Parameters
**name : str**  
The name of the coordinate (from BED file)

**maf_file: str**  
The path to the MAF file or the directory of MAF files.

**chrom : str**  
The chromosome of interest

**start : list[int]**  
The list of start locations of TE fragments on the chromosome

**end : list[int]**  
The list of end locations of TE fragments on the chromosome

**strand : str**  
The DNA strand of TE fragments, must be either '+' or '-'

**target_species : str, default: 'Homo_sapiens'**  
The species id used as a reference seqeunce for the alignment.

**count_arg : str, default: 'common'**  
The mode of base counting assigned to extract_maf().

- 'ref_freq': count reference base frequencies
- 'coverage': count total bases
- 'common_raw','common_freq','common_freq_nogap': count common base (raw count, freqeuncies, and frequencies without gaps)
- 'a','t','c','g': count ATCG
'base_count','base_freq','base_freq_nogap': count all bases (raw count, freqeuncies, and frequencies without gaps)
- 'raw': extract bases without counting 

**species_list : list[str], default: None**  
The list of species in the alignments to be counted or extracted. The names need to be the same as the names in the alignments

**e_value_df : pandas.DataFrame, default: None**  
The table of E-value calculated by teatime package as pandas.DataFrame, works as the species list to filter for species with significant alignment to the reference.

**internal_id_df: pandas.DaraFrame, default: None**  
The table of internal ID as pandas.DataFrame, works as the reference linking the coordinates of TEs to internal ID used by teatime outputs.

### Returns
**output_array: numpy.ndarray**  
The extracted base counts/frequencies from the MAF file.

## ma_mapper.extract_maf.get_spliced_mod()

`ma_mapper.extract_maf.get_spliced_mod(starts, ends, strand=1)`

The modified version of MafIO.MafIndex.get_spliced(). This version was made to solves the issue when entries with duplicated names are found in the alignment, collapsing them into ambiguous base with IUPAC instead of giving an error. Additionally, the output was modified to be pandas.DataFrame to fit the extract_maf pipeline.

### Parameters
**starts : list[int]**  
The list of start locations of TE fragments on the chromosome

**ends : list[int]**  
The list of end locations of TE fragments on the chromosome

**strand : str**  
The DNA strand of TE fragments, must be either 1 or -1

#### Returns
**df: pandas.DataFrame**  
The dataframe containing multiple alignments extracted from the MAF file. 