#!/usr/bin/env perl

#  
#  Copyright 2017 INRA-CIRAD
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, see <http://www.gnu.org/licenses/> or 
#  write to the Free Software Foundation, Inc., 
#  51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#  
#  

=pod

=head1 NAME

Multiple_GATK.pl : launch GATK on several samples in order to generate a fitered vcf

=head1 DESCRIPTION

This little tool will able to launch all the different steps of GATK to generatd a vcf that will be processed by vcftools

=head1 SYNOPSIS / USAGE

quality_begin_end_middle	-l input_file_list \
[-r reference_genome] \
[-p ploidy] \


=cut

use strict;
use warnings;
use utf8;
use Pod::Usage;
use Getopt::Long;

sub help
{ 
  my $msg = shift;

  print "-----------------------------------------------------------------------------\n"
  ."-- ERROR MESSAGE --\n"."$msg\n".
  "-----------------------------------------------------------------------------\n"
  if( defined  $msg );
  pod2usage(0);
}

help "Not enough arguments" unless (@ARGV);

=pod

=head1 OPTIONS

=head2 required parameters

=head3 At least one of this two parameter must be specified

=over 1

=item B<[-l]> ([input_file_list]): 

A list of the files that must be analysed. one line must contain only one file name and nothing more. Files must be in Bam or Sam format

=item B<[-b]> ([input_bed])

A bed file containing the targeted regions.

=back

=head2 optional arguments

=over 2

=item B<[-o]> ([output_file]) [./results_depth.csv]

The output primary csv file containing the results of the analyses will be stored here. It will contain the name of each individual sample, followed by informations about the number of sequences having a depth present in different intervals

=item B<[-g]> ([output_plot]) [False]

Save of plot showing the differences between the different analyzed parts of the sequences

=item B<[-n]> ([interval_size]) [10]

The size of an interval to analyse. 5 intervals are analysed, and can be resumed by : [0;n] ; [n;2n] ; ]2n;3n] ; ]3n;4n] ; ]4n;infinity[

=item B<[-i]> ([number_of_interval]) [5]

The number of interval to analyse. Must be an integer >= 2 

=item B<[-s]> ([sorting_col]) [1]

The column on which the output graph will be sorted. Must be between 1 and the number of intervals.

=back

=cut

my $man=0;
my $help=0;
#my $inputDepth="depth_by_sample.csv";
my $input;
my $outputCSV="outputDepthBySample.csv";
my $outputPlot='';
my $bedFile;
my $outputNamedCSV;
my $setSize=10; #Must be set by the user.
my $inputList;
my $sortingCol=1;
my $nbrInt=5;

GetOptions(
  "help|?"      => \$help,
  "man"         => \$man,
  "input_list|l=s"   => \$inputList,
  "interval_size|n=i" => \$setSize,
  "number_of_interval|i=i" => \$nbrInt,
  "bed|b=s"   => \$bedFile,
  "output|o=s"  => \$outputCSV,
  "outputPlot|g!" => \$outputPlot,
  "sorting_col|s=i" => \$sortingCol
) or pod2usage(2);

pod2usage(-verbose => 2) if ($man);
pod2usage(0) if $help;

pod2usage(
  -message => "Mandatory argument is missing or is not a file",
  -verbose => 1,
) if ((!defined $inputList || !-f $inputList) or (!defined $bedFile || !-f $bedFile));

pod2usage(
  -message => "The number of interval must be an integer >= 2 and it's not.",
  -verbose => 1,
) if ($nbrInt < 2 || int($nbrInt) != $nbrInt);

pod2usage(
  -message => "The size of the intervals must be > 0",
  -verbose => 1,
) if ($setSize <= 0);

pod2usage(
  -message => "The sorting column can't be > to the number of intervals",
  -verbose => 1,
) if ($sortingCol > $nbrInt);


