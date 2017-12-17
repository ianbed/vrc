#! /usr/bin/perl
use strict;
use warnings;

=head1

Author: Ian Beddows 12/17/2017.

Contact: beddowsi4@gmail.com.

Tested 17 Dec 2017 on MAC OS X El Capitan v10.11.2.

Dependencies:
        JSON.pm
        Getopt/Long.pm
        
This program takes a VCF file, extracts some of its information (see below) and returns it into a .txt file. The program also gets variant annotations from the ExAC database (exac.broadinstitute.org) if these are available.

To run this program on the command line type: perl beddows-vrc.pl -i path/to/vcf/file.vcf

The -o option can also be used to specify the outfile (defaults to results.txt)

Requires a VCFv4.2 file.

Requires the following IDs in column 8 of the VCF: TYPE, DP, AF, AO, and RO.

The returned txt file (-o option  or results.txt by default) has 1 row per variant & 8 columns:
    1) chr 
    2) pos 
    3) variant_type    
    4) variant_effect # from the ExAC database
    5) total_read_depth
    6) n_reads_supporting_variant  
    7) perc_reads_supporting_variant # relative to those supporting the reference
    8) variant_allele_freq            

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
print $out join("\t",'chr','pos','variant_type','variant_effect','total_read_depth','n_reads_supporting_variant','perc_reads_supporting_variant','variant_allele_freq'),"\n"; # print header

my $c=0; # total positions  
my $total_variants=0; # total number of alternate alleles
open(my $in,'<',$infile) or die print "The infile \"$infile\" was not found\n";
my $total = `cat $infile|grep -v '^#' | wc -l`; chomp($total); # get the total number of variants in the file
while(<$in>){ # open the file
    chomp;
    if(1..1){ # first line
        if($_!~/##fileformat=VCFv4.2/){ # must be the first line
            die print "The infile \"$infile\" is not a valid VCF4.2 file\n";
        }
    }elsif($_=~/^##/){
        ##reference=/primary/projects/bbc/references/human/annotation/hg19/special/GATK/hg19_broad_bundle/ucsc.hg19.fasta
        ##exac data from Gencode v19 transcript set   
    }elsif($_=~/^#/){
        # header line standard
    }else{
        # now dealing with the individual variants
        $c++;
        if($c % 100==0){ print "$c of $total positions done\n"; } # progress reporting...
        
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
        my %tmp = map { my ( $id, $value ) = split "\="; $id => $value } @info; # then turn this into a temporary hash (dictionary) with ID=value
        my %info=();
        foreach my $id (keys %tmp){ # after this, the ID will point to an array where array[i] is the ID information for alternate allele i
            @{$info{$id}} = split(',',$tmp{$id});
        }
        # now check to make sure all required INFO IDs are there (not really necessary, but may prevent issues if the investigator tries the script on an 'incomplete' VCF)
        foreach my $id ('TYPE','DP','RO','AO','AF'){
            if(!exists $info{$id}){die print "required INFO ID $id not in $infile ... exiting\n";}
        }
        
        for(my $i=0; $i<=$#alts; $i++){ # go through all alternate alleles & get the requested information
            $total_variants++;
            my $alt = $alts[$i]; # assign the first alternate allele            
            my $json = `curl -s http://exac.hms.harvard.edu/rest/variant/consequences/$chr-$pos-$ref-$alt`; # query ExAC using the api to get the variant effect
            my $variant_effect='NA'; # if no info in ExAC, then effect is not available
            if($json ne 'null' and $json ne '{}' ){ # some returns are empty like http://exac.hms.harvard.edu/rest/variant/consequences/chr4-107279450-C-T 
                my $decoded = decode_json($json); # decode the json into a hash reference using the JSON package
                my @effect_options = keys %$decoded; # the json hash keys are the different effect options (and now they are in an array)
                $variant_effect = get_most_deleterious(\@effect_options); # pass effect options array to a subroutine to get the most deleterious
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
                $info{'AF'}[0] # allele frequency
            ),"\n";         
        }
    }
}
print("\nFinished $0\nFinal data available in: $outfile\nData includes $c positions and $total_variants total alternate alleles\n");
#=======================================================================
#( Subroutines                  )
# ------------------------------------ 
#  o
#   o   \_\_    _/_/
#    o      \__/
#           (oo)_______
#           (__)       )/
#               ||----w |
#               ||     ||

sub get_most_deleterious {
    my $arr_ref = shift @_;
    my @choices = @$arr_ref;
    if($#choices==0){ # if there is only one choice, then return it (faster than the next option)
        return($choices[0]);
    }else{ # find the most deleterious of the effects
        # The ExAC consequences are from Variant Effect Predictor version 81 
        my @instances = ( # ordered array from most to least deleterious 
        # according to https://useast.ensembl.org/info/genome/variation/predicted_data.html#consequences
            'splice_acceptor_variant',
            'splice_donor_variant',
            'stop_gained',
            'frameshift_variant',
            'stop_lost',
            'start_lost',
            'transcript_amplification',
            'inframe_insertion',
            'inframe_deletion',
            'missense_variant',
            'protein_altering_variant',
            'splice_region_variant',
            'incomplete_terminal_codon_variant',
            'stop_retained_variant',
            'synonymous_variant',
            'coding_sequence_variant',
            'mature_miRNA_variant',
            '5_prime_UTR_variant',
            '3_prime_UTR_variant',
            'non_coding_transcript_exon_variant',
            'intron_variant',
            'NMD_transcript_variant',
            'non_coding_transcript_variant',
            'upstream_gene_variant',
            'downstream_gene_variant',
            'TFBS_ablation',
            'TFBS_amplification',
            'TF_binding_site_variant',
            'regulatory_region_ablation',
            'regulatory_region_amplification',
            'feature_elongation',
            'regulatory_region_variant',
            'feature_truncation',
            'intergenic_variant',
            'initiator_codon_variant' # this one was not in the documentation, but appears for e.g. chr8-104312432-A-G
        );
        my $i = 0;
        my %L = map { $_ => $i++ } @instances; # make a dictionary, L, ordered from the @instances array
        my @sorted = sort { $L{$a} <=> $L{$b} } @choices; # sort the choices based on their values in the L dictionary, e.g. $L{'splice_acceptor_variant'} = 0 and $L{'splice_donor_variant'} = 1 ... etc.
        return($sorted[0]); # return the most deleterious of the choices
    }
}
