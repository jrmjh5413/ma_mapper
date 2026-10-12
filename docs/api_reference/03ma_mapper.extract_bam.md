# ma_mapper.extract_bam
a pysam wrapper for BAM/SAM file data extraction
## list of functions
- ma_mapper.extract_bam.bam_io()
- ma_mapper.extract_bam.extract_bam()  
- ma_mapper.extract_bam.normal_array()

## ma_mapper.extract_bam.bam_io()

`ma_mapper.extract_bam.bam_io(coordinate_table, bam_file, bam_format,save_to_file=False, generate_new_id=False)`

Extract mapped reads from the input BAM file using the coodinates from the coordinate table. The user can choose between raw read counts and normalized reads. 

### Parameters
**coordinate_table : str,pd.DataFrame, default: None**  
The coordinate table for data extraction target. Accept both filepath or pandas dataframe. The table has to be in BED file format.

**bam_file : str**  
The path to the BAM file

**bam_format : str**  
output format of BAM data to be extracted.
- 'read' for raw read counts 
- 'normal' for Read coverage smoothed using a normal (Gaussian) curve
'read_forward','read_reverse','normal_forward','normal_reverse': extract forward or reverse reads only
- 'read_max','read_min', 'normal_max','normal_min': find minimum/maximum reads between forward and reverse reads
- 'read_sum' : sum of both forward and reverse reads

**save_to_file : bool,str, default: False**  
Save the output as a pickled file with lzma compression. The filepath can be specified, using save_to_file = True would save the list to the same directory as the coordinate table or current working directory.

**generate_new_id : bool, default: False**  
Generate new id for using index of the coordinate table. This parameter dictate the fragment joining behavior. By default, if the coordinate table have more than one coordinates that share the same 'name' or id, the function would join the results from these fragments.

### Returns
**extracted_data: list[numpy.ndarray]**  
The list of extracted data from BAM file

## ma_mapper.extract_bam.extract_bam()

`ma_mapper.extract_bam.extract_bam(bam_file, chrom ,start_list, end_list, strand, bam_format, offset=5, probe_length=100, smoothing_length=100)`

Extract mapped reads from the input BAM file from the specified genomic location. The user can choose between raw read counts and normalized reads. This is the wrapper for pysam package

### Parameters
**bam_file : str**  
The path to the BAM file

**chrom : str**  
The chromosome of interest

**start_list : list[int]**  
The list of start locations of TE fragments on the chromosome

**end_list : list[int]**  
The list of end locations of TE fragments on the chromosome

**strand : str**  
The DNA strand of TE fragments, must be either '+' or '-'

**bam_format : str**  
output format of BAM data to be extracted.
- 'read' for raw read counts 
- 'normal' for Read coverage smoothed using a normal (Gaussian) curve
- 'read_forward','read_reverse','normal_forward','normal_reverse': extract forward or reverse reads only
- 'read_max','read_min', 'normal_max','normal_min': find minimum/maximum reads between forward and reverse reads
- 'read_sum' : sum of both forward and reverse reads

**offset : int, default: 5**  
The positional offset for read mapping

**probe_length: int, default: 100**  
The region around each mapped read where read counts are distributed

**smoothing_length: int, default: 100**  
The spread of the Gaussian smoothing kernel applied to the read coverage

### Returns
**np_results: numpy.ndarray**  
The array of extracted data from BAM file

## ma_mapper.extract_bam.normal_array()

`ma_mapper.extract_bam.normal_array(width=1, sigma=1, odd=0)`

Create array of normal distribution

### Parameters
**width : int, default: 1**  
The width ot the normal array

**sigma : int, default: 1**  
The standard diviation of the normal distribution

**odd : int, default: 0**  
Control whether the output normal array has odd or even elements

### Returns
**normal_array: numpy.ndarray**  
The array of the normal distribution 
