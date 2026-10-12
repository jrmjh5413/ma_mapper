# ma_mapper.mapper
This module contains functions related to alignment parsing and data mapping.
## list of functions
- ma_mapper.mapper.map_and_overlay()
- ma_mapper.mapper.parse_and_filter()
- ma_mapper.mapper.map_data()  
- ma_mapper.mapper.extract_coordinate_from_alignment()
- ma_mapper.mapper.parse_alignment()
- ma_mapper.mapper.create_filter()
- ma_mapper.mapper.normalise()
- ma_mapper.mapper.load_alignment_file()
- ma_mapper.mapper.base_count()
- ma_mapper.mapper.import_parsed_alignment()
- ma_mapper.mapper.flank_sequence_io()

## ma_mapper.mapper.map_and_overlay()

`ma_mapper.mapper.map_and_overlay(alignment, data_file, data_format, filter=True,coordinate_file=None, extension_length=None, source_fasta=None, **kwargs)`

the main API of ma_mapper package, streamline genome wide data extraction and mapping/overlaying on TE multiple alignment

### Parameters:
**alignment : str**  
The path to the alignment FASTA file

**data_file : str**  
The path to genome-wide data file

**data_format : str**  
data/output format of data to be extracted  

BAM/SAM:  
'read_max','read_min','read_forward','read_reverse', 'normal_max',
'normal_min','normal_forward','normal_reverse','read_sum'

- 'read' for raw read counts and normal for read smooted by normal curve 
- 'read_forward','read_reverse','normal_forward','normal_reverse': extract forward or reverse reads only
- 'read_max','read_min', 'normal_max','normal_min': find minimum/maximum reads between forward and reverse reads
- 'read_sum' : summation of both forward and reverse reads

BIGWIG:  
'bigwig', BED: 'bed', MAF: 'maf', VCF: 'vcf'

**filter : bool,list[numpy.ndarray, numpy.ndarray], default: True**  
apply filter on the output matrix. By default, the module will use internal filter calculated from alignment file processing but custom filters can also be used by specifying the list of filters.

**coordinate_file : str,pd.DataFrame, default: None**  
The coordinate table for data extraction target. Accept both filepath or pandas dataframe. The table has to be in BED file format.

**extension_length : int, default: None**  
The extension length for data extraction, in case when observation beyond alignment is needed.

**source_fasta : str, default: None**  
The filepath to human genome sequence for additional sequnce extraction when extension_length is used.

### Returns:
**output_matrix_filtered: numpy.ndarray**  
The matrix of extracted genome-wide data with the same gap structure as TE multiple sequence alignment 

## ma_mapper.mapper.parse_and_filter()

`ma_mapper.mapper.parse_and_filter(alignment_file, filter=True, extension_length=None, preprocess_out=False,**kwargs)`

The wrapper for both parse_alignment() and create_filter(), parse alignment from FASTA format into a numerical matrix and extract coordinate table from alignment headers, then create filter from using alignment content, with additional function to apply filters on the matrix. 

### Parameters:

**alignment_file : str**  
The path to the alignment FASTA file

**filter : bool,list[numpy.ndarray, numpy.ndarray], default: True**  
create filters from parsed alignment content and apply the filter on the matrix. 

**extension_length : int, default: None**  
The extension length for seqeunce extraction, in case when observation beyond alignment is needed. the extracted sequence will be attached to both side of the alignment without filtering 

**preprocess_out : bool, default: False**  
export alignment and coordinate table before filtering, as well as the filters instead of the filtered ones

### Returns:

**alignment_output: numpy.ndarray**  
The numerical matrix of parsed alignment

**coordinate_table: pandas.DataFrame**  
The table of TE coordinates extracted from alignment file

**filters: list[numpy.ndarray, numpy.ndarray]**  
The list of row and column filers. Both are boolean numpy.1Darray to indicate row/column that need to be filtered out, only return when preprocess_out=True

## ma_mapper.mapper.map_data()

`ma_mapper.mapper.map_data(extracted_data, alignment_matrix, filter=True, nested_data=None,**kwargs)`

map/overlay extracted data onto the numerical alignment matrix by creating a blank container matrix then transferring extracted data from a list into the container using alignment matrix as a reference for gap positions

### Parameters

**extracted_data : str,numpy.ndarray**  
extracted data output from extraction modules in ma_mapper, could also be the filepath to pickled file of extracted data

**alignment_matrix : numpy.ndarray**  
The numerical matrix converted from multiple alignment in FASTA format

**filter : bool,list[numpy.ndarray, numpy.ndarray], default: True**  
create filters from parsed alignment content and apply the filter on the matrix.

### Returns:
**mapped_data_matrix: numpy.ndarray**  
The output matrix of extracted data with the same gap structure as the reference alignments. The matrix is ready to overlay with other data/alignment in visualization or downstream analyses.

biopython.SeqIOSeqRecord of multiple alignement from biopython.SeqIO.parse()

## ma_mapper.mapper.extract_coordinate_from_alignment()

`ma_mapper.mapper.extract_coordinate_from_alignment(alignment_file, save_to_file = False)`

Extract TE coordinates from headers of multiple aingnment FASTA. The headers must be generated with the following format: 

`>sequence_id::chrom:start-stop(strand)`

### Parameters

**aligment_file : str, biopython.SeqIO.SeqRecord**  
The path to the alignment FASTA file, or biopython.SeqIOSeqRecord of the alignment

**save_to_file : bool,str, default: False**  
Save the output table as a .csv file. The file path can be specified, using save_to_file = True would save the table to the current working directory

### Returns

**coordinate_table: pandas.DataFrame**  
The table in BED file format containing coordinates of seqeunces in the multiple alignment.

## ma_mapper.mapper.parse_alignment()

`ma_mapper.mapper.parse_alignment(alignment_file, save_to_file)`

Parse multiple alignment FASTA file into matrix and then convert the matrix into numerical using a hash table. 

### Parameters
**aligment_file : str, biopython.SeqIO.SeqRecord**  
The path to the alignment FASTA file, or biopython.SeqIOSeqRecord of the alignment

**save_to_file : bool,str, default: False**  
Save the output matrix as a pickled file. The file path can be specified, using save_to_file = True would save the table to the current working directory

### Returns

**alignment_matrix: numpy.ndarray**  
The numerical matrix of the multiple alignment. Ready for visualization by matplotlib.

## ma_mapper.mapper.create_filter()

`ma_mapper.mapper.create_filter(parsed_alignment, col_threshold = 0.50, col_content_threshold = 0.10, row_threshold = 0.50, save_to_file=False)`

Create column and row filter for numerical alignment, mainly for visualization as it is common for multiple alignment to have brigde structure where a seqeunce has a long line of unaligned portion that separated aligned blocks from each other. 

The function basically counts numbers of rows or columns against their total numbers.

### Parameters

**parse_alignment : numpy.ndarray**  
The numerical matrix of the multiple alignment

**col_threshold : float, default = 0.5**  
The threshold for fraction of non-gap contents of a column that appear in a possible maximum non-gap range (calculated by excusion of the gaps at the front and the back).

**col_content_threshold : float, default = 0.1**  
The threshold for fraction of non-gap contents of a column against total row.

**col_content_threshold : float, default = 0.1**  
The threshold for fraction of non-gap contents of a row against total column.

**save_to_file : bool,str, default: False**  
Save the output filters as a pickled file. The filepath can be specified, using save_to_file = True would save the filters to the current working directory

### Returns
**filters: list[numpy.ndarray, numpy.ndarray]**  
The list of numpy.1Darray of row filter and column filter for the numerical matrix of the multiple alignment.

## ma_mapper.mapper.normalise()

`ma_mapper.mapper.normalise(alignment_matrix, data_matrix, method = 'average')`

Normalize columns in multiple alignment with values adjusted to gaps in multiple aignment using numpy functions

### Parameters

**alignment_matrix : numpy.ndarray**  
The numerical matrix of the multiple alignment

**data_matrix : numpy.ndarray**  
The matrix of extracted genome-wide data with the same gap structure as TE multiple sequence alignment 

**method : str ('average','perc_coverage', 'median'), default: 'average'**  
The normalization method. 
-'average': averaging with nonzero alignment rows (excluding NA value)  
-'perc_coverage': find percentage of rows with nonzero data value compared to nonzero alignment rows 
-'median': find median value (excluding NA value)  

### Returns
**normalized_data: np.1Darray**  
The array of normalized data. The array is ready to be used as aggregated data matrix in visualization and downstream analyses

## ma_mapper.mapper.load_alignment_file()

`ma_mapper.mapper.load_alignment_file(alignment_file)`

The wrapper for biopython.SeqIO for FASTA file handling

### Parameters

**aligment_file: str, biopython.SeqIO.SeqRecord**  
The path to the alignment FASTA file, or biopython.SeqIOSeqRecord of the alignment

### Returns

**records: biopython.SeqIO.SeqRecord**  
biopython.SeqIO.SeqRecord containing alignment ready to be converted into a numerical matrix

## ma_mapper.mapper.base_count()

`ma_mapper.mapper.base_count(alignment_matrix)`

Count unique bases in numerical alignment matrix by column

### Parameters

**alignment_matrix : numpy.ndarray**  
The numerical matrix of the multiple alignment

### Returns
**base_count: numpy.ndarray**  
The base count from the multiple alignment, by column.

## ma_mapper.mapper.import_parsed_alignment()

`ma_mapper.mapper.import_parsed_alignment(parsed_alignment_file)`

Handle the alignment matrix that is stored in h5 format

### Parameters
**parsed_alignment_file : str**  
The path to the alignment matrix in h5 format

### Returns
**parse_alignment: numpy.ndarray**  
The numerical matrix of the multiple alignment

## ma_mapper.mapper.flank_sequence_io()

`ma_mapper.mapper.flank_sequence_io(coordinate_table, source_fasta, coordinate_table_out=False, extension_length=500)`

Extend the alignment matrix on both sides (front and back) by extracting flanking sequences and converting them into numerical matrices. The coordianate tables for flanking sequences can instead be exported in case they are needed for data extraction.

### Parameters
**coordinate_file : pd.DataFrame**  
The coordinate table for data extraction target. The table has to be in BED file format.

**source_fasta : str, default: None**  
The filepath to human genome sequence for sequnce extraction

**coordinate_table_out : bool, default: False**  
Export coordinate tables of the flanking sequences instead of alignment matrix

**extension_length : int, default: 500**  
The extension length for data extraction

### Returns
coordinate_table_out = False

**front_parsed_alignment, back_parsed_alignment: numpy.ndarray**  
The alignment matrices of the front and back flanking regions of the original matrix.

coordinate_table_out = True

**front_coordinate_table, back_coordinate_table: pandas.DataFrame**  
The tables in BED file format containing coordinates of flanking regions of the multiple alignment.