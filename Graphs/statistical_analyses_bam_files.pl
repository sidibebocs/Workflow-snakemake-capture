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
#
#FIXME : regarder la présence ou non des scripts.
#      : onglet dépendances.

=pod

=head1 NAME

Statistical_analysis_bam_file - analyse a set of bam files in order to get some stats about mapping

=head1 DESCRIPTION

Giving raw reads, bam files, modified bam files and a bed, output some stats and a graph about the quality of the different steps

=head1 SYNOPSIS / USAGE

quality_begin_end_middle	-r raw_read_directory \
[-rb raw_bam_directory] \
[-pb proceeded_bam_directory] \
[-b bed_file] \
[-g output_plot] 

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

=over 2

=item B<[-r]> ([raw_read_directory]): 

The directory containing unmodified read files in fastq format. 

=item B<[-rb]> ([raw_bam_directory])

The directory containing unmodified bam files given by the mapping step

=back

=head2 optional arguments

=over 3

=item B<[-pb]> ([proceed_bam_directory])

Directory of the ready for further analysis bams, only required if some treatment have been applied on the raw bams.

=item B<[-b]> ([bed_file])

Path to a bed file containing regions of interest. Will add on the graph an "on target" information containing the amount of reads that have been mapped in the intervals specified in the bed

=item B<[-g/-nog]> ([output_plot]) [False]

If true, will output a plot (in pdf) showing the differences between the different steps of the analysis. [default = false]

=back

=cut

my $man=0;
my $help=0;
my $raw_read_directory;
my $raw_bam_directory;
my $proceed_bam_directory;
my $bed;
my $outputPlot;
my $verbose;

my $scriptPath=$0;
$scriptPath =~ s/[\w\.]*$//i; 
if ($scriptPath eq ""){
  $scriptPath=".";
}

GetOptions(
  "help|?"      => \$help,
  "man"         => \$man,
  "raw_read_directory|r=s"   => \$raw_read_directory,
  "raw_bam_directory|rb=s"  => \$raw_bam_directory,
  "proceed_bam_directory|pb=s"  => \$proceed_bam_directory,
  "bed_file|b=s"  => \$bed,
  "outputPlot|g!" => \$outputPlot,
  "verbose|v!" => \$verbose
) or pod2usage(2);

pod2usage(-verbose => 2) if ($man);
pod2usage(0) if $help;

pod2usage(
  -message => "Mandatory argument is missing",
  -verbose => 1,
) if !defined ($raw_read_directory && $raw_bam_directory);

pod2usage(
  -message => "The read directory or the bam directory are probably not directories",
  -verbose => 1,
) if (!-d $raw_read_directory || !-d $raw_bam_directory);

pod2usage(
  -message => "The proceed bam directory is probably not a directory",
  -verbose => 1,
) if  (defined $proceed_bam_directory && !-d $proceed_bam_directory);

pod2usage(
  -message => "The bed file is not a file",
  -verbose => 1,
) if  (defined $bed && !-f $bed);

pod2usage(
  -message => "One of the following file is missing in the script directory : autoGraph.R ; graphe_aire_sum_auto.R ; statistical_analyses_bam_files.sh",
  -verbose => 1,
) if  (!-f "$scriptPath/autoGraph.R" || !-f "$scriptPath/graphe_aire_sum_auto.R" || !-f "$scriptPath/statistical_analyses_bam_files.sh");

my $ligneCmd='';


$ligneCmd="bash $scriptPath/statistical_analyses_bam_files.sh -a $raw_read_directory -b $raw_bam_directory";


if (defined $proceed_bam_directory){
  $ligneCmd.=" -c $proceed_bam_directory";
}

if (defined $bed && -e $bed){
  $ligneCmd.=" -d $bed ";
}

if (defined $outputPlot){
  $ligneCmd.=" -e";
}

if ($verbose){
  print "command line built : $ligneCmd\n\n";
}

system "$ligneCmd";

=head1 DEPENDENCIES

Samtools

=head1 INCOMPATIBILITIES

None

=head1 AUTHORS

=over2

=item Alexandre SORIANO (INRA), soriano.alexandre30@gmail.com

=back

=head1 VERSION

1.0.0

=head1 DATE

17.O5.2017

=cut
