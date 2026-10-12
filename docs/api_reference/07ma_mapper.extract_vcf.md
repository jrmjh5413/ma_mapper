# ma_mapper.extract_vcf
a cyvcf2 wrapper for VCF file data extraction
## list of functions
- ma_mapper.extract_vcf.vcf_io()
- ma_mapper.extract_vcf.extract_vcf()

## ma_mapper.extract_vcf.vcf_io()

`ma_mapper.extract_vcf.vcf_io(coordinate_table, vcf, query_key = 'AF', vcf_format = 'gnomad', save_to_file = False, generate_new_id = False)`

Extract data from the input VCF file using the coodinates from the coordinate table. As of now, the module supports gnomAD dataset.

### Parameters
**coordinate_table : str,pd.DataFrame, default: None**  
The coordinate table for data extraction target. Accept both filepath or pandas dataframe. The table has to be in BED file format.

**vcf : str**  
the directory of VCF files

**query_key : str, default: 'AF'**  
The keyword to access specific data field in VCF files

**vcf_format : str, default: 'gnomad'**  
The format or the name for VCF databases. As VCF file format offers a certain degree of freedom in data fields, different databases will likely store data in different ways.

**save_to_file : bool,str, default: False**  
Save the output as a pickled file with lzma compression. The filepath can be specified, using save_to_file = True would save the list to the same directory as the coordinate table or current working directory.

**generate_new_id : bool, default: False**  
Generate new id for using index of the coordinate table. This parameter dictate the fragment joining behavior. By default, if the coordinate table have more than one coordinates that share the same 'name' or id, the function would join the results from these fragments.

### Returns
**vcf_out: list[numpy.ndarray]**  
The list of extracted data from VCF file

## ma_mapper.extract_vcf.extract_vcf()

`ma_mapper.extract_vcf.extract_vcf(vcf_file, chrom, start_list, end_list, strand, query_key)`

Extract data from the input VCF file from the specific genomic location. This is the wrapper for cyvcf2.

### Parameters
**vcf_file: str**  
The path to the VCF file.

**chrom : str**  
The chromosome of interest

**start_list : list[int]**  
The list of start locations of TE fragments on the chromosome

**end_list : list[int]**  
The list of end locations of TE fragments on the chromosome

**strand : str**  
The DNA strand of TE fragments, must be either '+' or '-'

**query_key : str, default: 'AF'**  
The keyword to access specific data field in the VCF file

#### Returns
**window_out: np.ndarray**  
The extracted data from the VCF file.