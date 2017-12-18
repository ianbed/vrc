#! /usr/bin/perl
use strict;
use warnings;

=head1
Documentation for beddows-vrc.pl.

Author: Ian Beddows 12/17/2017.

Contact: beddowsi4@gmail.com.

Tested 17 Dec 2017 on MAC OS X El Capitan v10.11.2.

The program has two perl module dependencies:
 
1) JSON.pm

2) Getopt/Long.pm
        
This program takes a VCF file, extracts some of its information (see below) and returns it into a .txt file.
Two pieces of information from the ExAC database (exac.broadinstitute.org) are also reported: 
    
1) The most deleterious variant effect.

2) Whether or not the variant is found in the ExAC database.

To run this program on the command line type: perl beddows-vrc.pl -i path/to/vcf/file.vcf

The -o option can also be used to specify the outfile (defaults to results.txt)

Requires a VCFv4.2 file.

Requires the following IDs in column 8 of the VCF: TYPE, DP, AF, AO, and RO.

The results.txt file has 1 row per variant and 9 columns:

    1) chr 
    
    2) pos 
    
    3) variant_type    
    
    4) variant_effect # from the ExAC database
    
    5) total_read_depth
    
    6) n_reads_supporting_variant  
    
    7) perc_reads_supporting_variant # relative to those supporting the reference
    
    8) variant_allele_freq 
    
    9) any_covered # from the ExAC database (true if this variant is found in individuals in the database)
     
=cut
use Getopt::Long;
use JSON qw( decode_json ); 

my $usage = <<EOF;
To run this program on the command line type: perl beddows-vrc.pl -i path/to/vcf/file.vcf
REQUIRED OPTIONS:
-infile = path to the VCFv4.2 file with the data to be reworked
OPTIONAL:
-outfile = user may specify the outfile (defualt = results.txt)
-h|help = print usage
EOF

my($help,$infile,$outfile);
GetOptions(
	'infile=s' => \$infile,	# string
	'outfile=s' => \$outfile,	# string
	'h|help' => \$help	# flag for help
);
if(defined($help)){die print "HELP REQUESTED:\n",$usage;}
if(!defined($infile)){die print "\nNeed to define a VCF file using -i /path/to/vcf/file.vcf\n\n$usage\n";}
if(!defined($outfile)){$outfile = 'results.txt'} # default
if($outfile!~/\.txt$/){$outfile = $outfile . '.txt'} # add .txt ending to outfile if needed
open(my $out,'>',$outfile) or die print "\"$outfile\" was not opened ... exiting\n";
print $out join("\t",'chr','pos','variant_type','variant_effect','total_read_depth','n_reads_supporting_variant','perc_reads_supporting_variant','variant_allele_freq','any_covered'),"\n"; # print header

my $c=0; # total positions  
my $total_variants=0; # total number of alternate alleles
open(my $in,'<',$infile) or die print "The infile \"$infile\" was not found\n";
my $total_positions = `cat $infile|grep -v '^#' | wc -l`; chomp($total_positions); # get the total number of positions in the file
while(<$in>){ # open the file
    chomp;
    if(1..1){ # first line
        if($_!~/##fileformat=VCFv4.2/){ # must be the first line
            die print "The infile \"$infile\" is not a valid VCFv4.2 file\n";
        }
    }elsif($_=~/^##/){
        ##reference=/primary/projects/bbc/references/human/annotation/hg19/special/GATK/hg19_broad_bundle/ucsc.hg19.fasta
        ##exac data from Gencode v19 transcript set   
    }elsif($_=~/^#/){
        # header line standard
    }else{
        # now dealing with the individual variants
        $c++;
        if($c % 100==0){ print "$c of $total_positions positions done\n"; } # progress reporting...
        
        my @line = split('\s+',$_); # first split the line into an array
        my $chr = shift @line; # and get the individual fields as strings
        my $pos = shift @line;
        my $id = shift @line;
        my $ref = shift @line;
        my $_alt = shift @line;
        my $qual = shift @line;
        my $filter = shift @line;
        my $_info = shift @line; 
        my $format = shift @line; 
        # now @line has only the individual samples remaining
        
        my @alts = split(',',$_alt); # for dealing with multiple alternate alleles, split them into an array
        my @info = split(';',$_info); # to access the INFO, first make an array with all key/value pairs
        my %tmp = map { my ( $id, $value ) = split "\="; $id => $value } @info; # then turn this array into a temporary hash (dictionary) with ID=value
        my %info=();
        foreach my $id (keys %tmp){ # after this, the ID will point to an array where array[i] is the ID information for alternate allele i
            @{$info{$id}} = split(',',$tmp{$id});
        }
        # now check to make sure all required INFO IDs are there (not really necessary, but may prevent issues if the investigator tries the script on a different file)
        foreach my $id ('TYPE','DP','RO','AO','AF'){
            if(!exists $info{$id}){die print "required INFO ID $id not in $infile ... exiting\n";}
        }
        
        for(my $i=0; $i<=$#alts; $i++){ # go through all alternate alleles & get the requested information
            $total_variants++;
            my $alt = $alts[$i]; # assign the first alternate allele            
            my $json = `curl -s http://exac.hms.harvard.edu/rest/variant/ordered_csqs/$chr-$pos-$ref-$alt`; # query ExAC using the api to get the most deleterious variant effect
            my $any_covered = `curl -s http://exac.hms.harvard.edu/rest/variant/any_covered/$chr-$pos-$ref-$alt`; # query ExAC again
            my $variant_effect='null'; # if no info in ExAC, then effect the is not available

            if($json ne 'null'){ # returns a null if the variant is not in the db and for a few that are in the db
                my $decoded = decode_json($json); # decodes the json into an array reference using the JSON package
                $variant_effect = $$decoded[1] or $variant_effect = 'null'; # the first array element is empty & the second is the most deleterious (https://useast.ensembl.org/info/genome/variation/predicted_data.html#consequences) 
                # many positions (e.g. chr4-107279450-C-T and chr4-100048506-GCAT-ACAT) are in the db, but have no effects given, therefore the 'or $variant_effect eq null'
            }
            
            
            
            my $perc_supp_alt; # get the percentage of reads supporting the variant versus those supporting the reference
            if($info{'RO'}[0]!=0){ # if there are NOT 0 reads supporting the reference
                $perc_supp_alt = $info{'AO'}[$i] / ($info{'RO'}[0]+$info{'AO'}[$i]) * 100; # do not divide by DP because not all reads need to support ref or this alt allele
            }elsif($info{'AO'}[$i]>0){ # if 0 reads support the ref and one or more support the alt, then we will say all reads support the alt
                $perc_supp_alt = 100;
            }else{
                die print "0 reads for the reference AND 0 reads for the alt allele\n"; # bigger problems if this happens
            }
            $perc_supp_alt = sprintf('%.2f',$perc_supp_alt); # round answer to two decimal points
            
            # we now have all of the requested information and can print to the outfile
            print $out join("\t",
                $chr,
                $pos,
                #$ref, # could include these if you wanted
                #$alt,
                $info{'TYPE'}[$i],
                $variant_effect,
                $info{'DP'}[0], # total read depth at the locus
                $info{'AO'}[$i],
                $perc_supp_alt,
                $info{'AF'}[0], # allele frequency
                $any_covered
            ),"\n";         
        }
    }
}
print("\nFinished $0\nFinal data available in: $outfile\nData includes $c positions and $total_variants total alternate alleles\n");

#=======================================================================
#( Subroutines                  )
# ------------------------------- 
#  o
#   o   \_\_    _/_/
#    o      \__/
#           (oo)_______
#           (__)       )/
#               ||----w |
#               ||     ||
