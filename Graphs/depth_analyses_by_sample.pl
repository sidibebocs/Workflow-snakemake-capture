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

depth_analyses_bam_files.pl - analyse the depth of some BAM file in order to see the mean quality

=head1 DESCRIPTION

Analyse some provided bam files and outplut informations about the average coverage for the different individuals analized.

=head1 SYNOPSIS / USAGE

quality_begin_end_middle	-l input_file_list \
[-b input_bed_file] \
[-o output_file] \
[-n interval_size] \
[-i number_of_interval] \
[-g output_plot] \
[-s sorting_col] \

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

open (my $fileList, '<', $inputList);
open (my $output, '>', $outputCSV);

my $line;
my $compteur=0;
my $intervalLength=0;
my $numInterval=0;
my $scriptPath=$0;
my $grapheName="depth_by_sample_n".$setSize."_s".$sortingCol."_i".$nbrInt.".";
my $inputDepth="samtools_depth.csv";

my $fourSetSize=$setSize*($nbrInt-1);
my @tabLine;
my @tabMean;
my @tabResults;
my @intervalDepth;
my @tabName;

while (<$fileList>){
  $line=$_;
  chomp($line);
  $tabName[$.-1]=$line;
}

##########################################
# Samtools depth

if (!-T $inputDepth){ #If the file does not exist, use samtools -aa on her.
  print "samtools depth on $inputList and $bedFile...\n\n";
  system "samtools depth -aa -b $bedFile -f $inputList > $inputDepth";
}
else{
  print "The file $inputDepth already exists or is not a appropriate file. If you are not certain that this is the required file, please delete it or move it.\n\n";
}

open (my $inputDepthFile, '<', $inputDepth);

print "Getting intervals from samtools depth output...\n\n"; #Create a bed using the file generated before.
get_interval_from_samtools_depth($inputDepth, "interval.bed");

print "Storing Bed in Array...\n\n"; #Store the different intervals of this Bed in an array.
my @bedArray = store_bed("interval.bed");

##########################################
# Analyses

print "Processing...\n";

while (<$inputDepthFile>){
  if ($.%250000==0){
    print STDERR "$. positions...\n";
  }
  $line=$_;
  chomp($line);
  @tabLine = split /\s/, $line;
  if ($compteur==0){ #Mean that we are entering a new interval
    for (my $i=0;$i<scalar(@tabLine-2);$i++){
      $tabMean[$i]=0;
    }
    $compteur = $bedArray[$numInterval][2] - $bedArray[$numInterval][1];
    $numInterval++;
    $intervalLength=$compteur;
    if ($. == 1){ #Array initialisation, depend of the number of line stored in tabline.
      for (my $i=0;$i<scalar(@tabLine)-2;$i++){
        for (my $j=0;$j<$nbrInt;$j++){
          $tabResults[$i][$j]=0;
        }
      }
    }
    next;
  }
  else{#Here we will fill the array tabMean, which will wontaint the mean for each individual for one interval
    for (my $i=2;$i<scalar(@tabLine);$i++){
      $tabMean[$i-2]+=$tabLine[$i];
    }
    $compteur--;
    if ($compteur==0){ #Mean that we reached the end of an interval, so we need to store the calculated values.
      for (my $i=0;$i<scalar(@tabMean);$i++){
        $tabMean[$i]=$tabMean[$i]/$intervalLength;
        if ($tabMean[$i] < $fourSetSize){
          #print "poney\n";
          $tabResults[$i][int($tabMean[$i]/$setSize)]++;
        }
        else{
          #print "petit cheval\n";
          $tabResults[$i][$nbrInt-1]++;
        }
      }
    }
  }
}

##########################################
# Output CSV Part

print $output "Name\t[".$fourSetSize."-inf[";
for (my $i=$nbrInt-1;$i>0;$i--){
  print $output "\t[".$setSize*$i."-".($setSize*($i-1))."[";
}
print $output "\n";


for (my $i=0;$i<scalar(@tabName);$i++){ #Output results in a csv file
  print $output "$tabName[$i]\t";
  for (my $j=$nbrInt-1;$j>=0;$j--){
    if ($j != 0){
      print $output "$tabResults[$i][$j]\t";
    }
    else{
      print $output "$tabResults[$i][$j]\n";
    }
  }
}

##########################################
# Ouput Graphe


if ($outputPlot){
  $scriptPath =~ s/[\w\.]*$//i; 
  if ($scriptPath eq ""){
    $scriptPath=".";
  }
  my $RscriptPath = $scriptPath."autoGraph_graphe_aire.R";
  my $RfunctionGraphePath = $scriptPath."graphe_aire.R";
  my $Xtitle="Sample";
  my $Ytitle="Targeted_region";
  my $grapheTitle="Intervals_of_coverage_for_targeted_regions_depending_on_sample";

  print "\nCreating plot...\n\n";

  #system "Rscript $RscriptPath $outputCSV $sortingCol $grapheName $RfunctionGraphePath";
  system "Rscript $RscriptPath $outputCSV $sortingCol $grapheName $Xtitle $Ytitle $grapheTitle $RfunctionGraphePath";

  print "The plot has been saved as \'$grapheName\'\n\n";
}

close $inputDepth;

#####################
##### Functions #####
#####################

sub store_bed { #Store a bed file in an array : 0 = chr ; 1 = begin ; 2 = end
  my ($bedFile) = @_;
  my @genomic_intervals;
  open (my $bed, '<', $bedFile);
  while (<$bed>){
    chomp($_);
    ($genomic_intervals[$.-1][0], $genomic_intervals[$.-1][1], $genomic_intervals[$.-1][2]) = split /\t/, $_;
  }
  return @genomic_intervals;
}

#FIXME
#We have +1 at the end of the interval and +1 at the begining
sub get_interval_from_samtools_depth { #Output a bed file using the output of samtools depth as an input.
  my ($fileInput, $outputFile) = @_;
  my $line;
  my $numero=0;
  my $compteur=0;
  my $storeChr=0;
  my $storeBeg;
  my $prevNum=0;
  my $storeEnd; 

  open (my $depthFile, '<', $fileInput);
  open (my $outputBed, '>', $outputFile);

  while(<$depthFile>){
    $line = $_;
    chomp($line);
    ($numero) = (split /\t/, $line)[1, 1];
    if ($compteur == 0){
      ($storeChr, $storeBeg) =(split /\t/, $line)[0, 1];
      if ($. != 1){
        $storeBeg--;
      }
      $prevNum=$numero;
      $prevNum++;
      $compteur++;
      next;
    }
    if ($numero == $prevNum){ #Same interval
      $prevNum++;
      $compteur++;
    }
    else{
      $storeEnd=$prevNum-1; #Use to be -2, don't know why.
      $compteur=0;
      print $outputBed "$storeChr\t$storeBeg\t$storeEnd\n";
    }
  }
  print $outputBed "$storeChr\t$storeBeg\t$numero\n";
}

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

02.O6.2017

=cut
