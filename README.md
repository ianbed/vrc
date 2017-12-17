Documentation for beddows-vrc.pl.

Author: Ian Beddows 12/17/2017.
Contact: beddowsi4@gmail.com.
Tested 17 Dec 2017 on MAC OS X El Capitan v10.11.2.

The program has 2 dependencies: JSON.pm and Getopt/Long.pm
        
Description:
This program takes a VCF file, extracts some of its information (see below) and returns it into a .txt file. The program also gets variant annotations from the ExAC database (exac.broadinstitute.org) if these are available.
        
To run this program on the command line type: perl beddows-vrc.pl -i path/to/vcf/file.vcf

The -o option can also be used to specify the outfile (defaults to results.txt)
        
Requires a VCFv4.2 file.
Requires the following IDs in column 8 of the VCF: TYPE, DP, AF, AO, and RO.
        
The returned txt file has 1 row per variant & 8 columns:

1) chr 

2) pos 

3) variant_type    

4) variant_effect # from the ExAC database

5) total_read_depth

6) n_reads_supporting_variant  

7) perc_reads_supporting_variant # relative to those supporting the reference

8) variant_allele_freq     
