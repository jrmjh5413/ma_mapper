# ma_mapper.extract_bed
a pybedtools wrapper for BED file data extraction
## list of functions
- ma_mapper.extract_bed.bed_io()
- ma_mapper.extract_bed.intersect_process()
- ma_mapper.coordinate_table_to_bed()

## ma_mapper.extract_bed.bed_io()

`ma_mapper.extract_bed.bed_io(coordinate_table, bed, save_to_file = False, generate_new_id=False, strand_overlap=True)`

Extract scores from the input BED file from the overlaps between the coodinates from the coordinate and from the BED file. This is the wrapper for pybedtools package.

### Parameters
**coordinate_table : str,pd.DataFrame, default: None**  
The coordinate table of data extraction targets. Accept both filepath or pandas dataframe. The table has to be in BED file format.

**bed : str, pandas.DataFrame, pybedtools.BedTool**  
The path to the BED file or the BED table as pandas.DataFrame object.

**save_to_file : bool,str, default: False**  
Save the output as a pickled file with lzma compression. The filepath can be specified, using save_to_file = True would save the output to the same directory as the coordinate table or current working directory.

**generate_new_id : bool, default: False**  
Generate new id for using index of the coordinate table. This parameter dictate the fragment joining behavior. By default, if the coordinate table have more than one coordinates that share the same 'name' or id, the function would join the results from these fragments.

**strand_overlap : bool, default: True**  
Include DNA strands of TE fragments as overlap condition. Need to be disabled in the case where DNA strand is not specified (. instead of +/-)

### Returns
**extracted_data: list[numpy.ndarray]**  
The list of extracted data from BED file

## ma_mapper.extract_bed.intersect_process()

`ma_mapper.extract_bed.intersect_process(intersect_df, unique_id)`

Create the array of scores for the overlapping ranges between two BED files

### Parameters
**intersect_df : pandas.DataFrame**  
The pandas.DataFrame of the overlapping ranges between two BED files, which is the output from pybedtools.BedTool.intersect()

**unique_id : str**  
The ID of the TE fragment of interest

### Returns
**bed_output_array: np.ndarray**  
The numpy.ndarray of accumulated score from the input BED file that overlaps with the TE of interest.

## ma_mapper.coordinate_table_to_bed()

`ma_mapper.coordinate_table_to_bed(coordinate_table, save_to_file=False, export_as_bedtool=False)`

Convert the coordinate table that might not have score column to fit the BED format into a proper BED table

### Parameters
**coordinate_table : str, pandas.DataFrame**  
The coordinate table of data extraction targets. Accept both filepath or pandas dataframe. 

**save_to_file : bool, default: False**  
Save the output table in TSV format. The filepath can be specified, using save_to_file = True would save the table to the same directory as the coordinate table or current working directory.

**export_as_bedtool : bool, default: False**  
Export the output as pybedtools.BedTool for further bedtools operations.

### Returns
**coordinate_bed: pandas.DataFrame, pybedtools.BedTool**  
The coordinate table in BED format.

