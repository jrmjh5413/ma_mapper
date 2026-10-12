# ma_mapper.sequence_alignment
a wrapper for MAFFT, for seqeunce extraction and multiple alignment 
## list of functions
- ma_mapper.sequence_alignment.seqeunce_io()
- ma_mapper.sequence_alignment.mafft_align()
- ma_mapper.sequence_alignment.extract_coord_from_repeatmasker_table()
- ma_mapper.sequence_alignment.extract_sequence()
- ma_mapper.sequence_alignment.parse_array_to_sequence()

## ma_mapper.sequence_alignment.seqeunce_io()

`ma_mapper.sequence_alignment.seqeunce_io(coordinate_table,source_fasta, save_to_file=False, generate_new_id = False)`

Extract DNA seqeunces from the genome sequence using the coordinates from the coordinate table. This is the wrapper for Bio.SeqIO

### Parameters
**coordinate_table : str,pd.DataFrame, default: None**  
The coordinate table for sequence extraction target. Accept both filepath or pandas dataframe. The table has to be in BED file format.

**source_fasta : str, default: None**  
The filepath to genome sequence for sequence extraction

**save_to_file : bool,str, default: False**  
Save the output sequence as a FASTA file. The filepath can be specified, using save_to_file = True would save the list to the same directory as the coordinate table or current working directory.

**generate_new_id : bool, default: False**  
Generate new id for using index of the coordinate table. This parameter dictate the fragment joining behavior. By default, if the coordinate table have more than one coordinates that share the same 'name' or id, the function would join the results from these fragments.

### Returns
**seq_records: Bio.SeqRecord.SeqRecord**  
Bio.SeqRecord.SeqRecord object containing extracted DNA sequences and their IDs from the coordinate table.

## ma_mapper.sequence_alignment.mafft_align()

`ma_mapper.sequence_alignment.mafft_align(input_filepath, nthread = None, nthreadtb = None, nthreadit =None, output_filepath = None, mafft_arg = '--treeout --reorder')`

The wrapper for MAFFT multiple alignment program. 

### Parameters
**input_filepath : str**  
The path to the input DNA sequences in FASTA format

**nthread, nthreadtb, nthreadit : int, default: None**  
The arguments for parallelization implementation of MAFFT.

**output_filepath : str, default: None**  
The path of the output file, if None the function will use the same directory as the input.

**mafft_arg : str, default: '--treeout --reorder'**  
The additional MAFFT commands/argument for fine-tuning.

## ma_mapper.sequence_alignment.extract_coord_from_repeatmasker_table()

`ma_mapper.sequence_alignment.extract_coord_from_repeatmasker_table(subfamily, repeatmasker_table, internal_id_table=None, save_to_file=False)`

Extract TE coordinates from repeatmasker table using a subfamily name as a filter, then convert and export the table into a BED file. 

### Parameters
**subfamily : str, list[str]**  
The name of the subfamily of interest. It is possible to extract more than one subfamily at once.

**repeatmasker_table : str, pandas.DataFrame**  
The path to the output table from RepeatMasker or its pandas.DataFrame.

**internal_id_table : str, list[str], pandas.DataFrame, list[pandas.DataFrame]**  
The path to the internal id table or its pandas.DataFrame from teatime package. If this option is used, the TEs will be assigned with names that can be use internally with teatime workflow.

**save_to_file : bool, default: False**  
Save the output table as a .csv file. The file path can be specified, using save_to_file = True would save the table to the current working directory 

### Returns
**coords: pandas.DataFrame**  
The pandas.DataFrame for the coordinate table of TE subfamily of interest.

## ma_mapper.sequence_alignment.extract_sequence()

`ma_mapper.sequence_alignment.extract_sequence(chrom, start, end,strand, records)`

Extract the DNA sequence string on the specific genomic location from Bio.SeqRecord.SeqRecord 

### Parameters
**chrom : str**  
The chromosome of interest

**start : int**  
The start location of TE fragments on the chromosome

**end : list**  
The end location of TE fragments on the chromosome

**strand : str**  
The DNA strand of TE fragments, must be either '+' or '-'

**records : Bio.SeqRecord.SeqRecord**  
The Bio.SeqRecord.SeqRecord object of the genome of interest

### Returns
**seq_string: str**  
The string of DNA sequeunce extracted 

## ma_mapper.sequence_alignment.parse_array_to_sequence()

`ma_mapper.sequence_alignment.parse_alignment_to_fasta(parsed_alignment,name=None, save_to_file=False)`

Convert the alignment matrix (parsed alignment) from ma_mapper.mapper into the list of Bio.SeqRecord.SeqRecord, mainly for exporting the matrix as a FASTA file

## Parameters
**parsed_alignment : numpy.ndarray**  
The parsed alignment matrix from ma_mapper.mapper

**name : str, default: None**  
The name of DNA seqeunces, it will be attached by running numbers

**save_to_file : bool,str, default: False**  
Save the output sequence as a FASTA file. The filepath can be specified, using save_to_file = True would save the list to the same directory as the coordinate table or current working directory.

### Returns
**dna_seqeunces: list[Bio.SeqRecord.SeqRecord]**
The list of Bio.SeqRecord.SeqRecord of DNA sequences converted from the parsed alignment.