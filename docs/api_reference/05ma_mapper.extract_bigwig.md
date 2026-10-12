# ma_mapper.extract_bigwig
a pyBigWig wrapper for BIGWIG/WIG file data extraction
## list of functions
- ma_mapper.extract_bigwig.bigwig_io()
- ma_mapper.extract_bigwig.extract_bigwig()

## ma_mapper.extract_bigwig.bigwig_io()

`ma_mapper.extract_bigwig.bigwig_io(coordinate_table, bigwig, save_to_file = False, generate_new_id=False)`

Extract data from the input BIGWIG file using the coodinates from the coordinate table. 

### Parameters
**coordinate_table : str,pd.DataFrame, default: None**  
The coordinate table of data extraction targets. Accept both filepath or pandas dataframe. The table has to be in BED file format.

**bigwig : str**  
The path to the BIGWIG file.

**save_to_file : bool,str, default: False**  
Save the output as a pickled file with lzma compression. The filepath can be specified, using save_to_file = True would save the output to the same directory as the coordinate table or current working directory.

**generate_new_id : bool, default: False**  
Generate new id for using index of the coordinate table. This parameter dictate the fragment joining behavior. By default, if the coordinate table have more than one coordinates that share the same 'name' or id, the function would join the results from these fragments.

### Returns
**extracted_data: list[numpy.ndarray]**  
The list of extracted data from BIGWIG file

## ma_mapper.extract_bam.extract_bigwig()

`ma_mapper.extract_bigwig.extract_bigwig(bigwig_file, chrom ,start_list, end_list, strand)`

Extract data signal the input BIGWIG file from the specified genomic location. This is the wrapper for pyBigWig package.

### Parameters
**bam_file : str, pyBigWig.bigWigFile**  
The path to the BIGWIG file or pyBigWig.bigWigFile.

**chrom : str**  
The chromosome of interest

**start_list : list[int]**  
The list of start locations of TE fragments on the chromosome

**end_list : list[int]**  
The list of end locations of TE fragments on the chromosome

**strand : str**  
The DNA strand of TE fragments, must be either '+' or '-'

### Returns
**bigwig_out: numpy.ndarray**  
The array of extracted data from BIGWIG file


