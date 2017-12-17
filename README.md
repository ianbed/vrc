documentation for beddows-vrc.pl
a program written for the vrc coding challenge.

author: Ian Beddows 12/17/2017
contact: beddowsi4@gmail.com
tested 17 Dec 2017 on MAC OS X El Capitan v10.11.2

dependencies:
        JSON.pm
        Getopt/Long.pm
        
description:
        This program takes a VCF file, extracts some of its 
        information (see below) and returns it into a .txt file.
        
        To run this program on the command line type:
        "perl beddows-vrc.pl -i path/to/vcf/file.vcf"
        
        The -o option can also be used to specify the outfile (defaults to results.txt)
        
        Requires a VCFv4.2 file.
        Requires the following IDs in column 8 of the VCF: TYPE, DP, AF, AO, and RO.
        Requires JSON.pm and Getopt::Long.pm.
        
        The returned txt file has 1 row per variant & 8 columns:
            chr 
            pos 
            variant_type    
            variant_effect # from the ExAC database
            total_read_depth
            n_reads_supporting_variant  
            perc_reads_supporting_variant # relative to those supporting the reference
            variant_allele_freq     
