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


#TODO
#Able the user to set the different intervals.
#Corect the storeBed sub (which is now pretty usable but fail to recover 1pb long intervals (which are useless and unrealistic)

=pod

=head1 NAME

quality_begin_end_middle.pl - analyse the depth of some BAM file begining, middle and end.

=head1 DESCRIPTION

Analyse some provided bam files in order to see the differences between the depth of the beginning of targeted sequences, their middle and their end.

=head1 SYNOPSIS / USAGE

quality_begin_end_middle	-b input_bed_file \
[-f input_bam_list] \
[-o output_file] \
[-n bases_to_get] \
[-l interval_size] \
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

=over 1

=item B<[-b]> ([input_bed]): 

A bed file containing the regions of interest.

=item B<[-f]> ([input_bam_list]):

A file containing the path to every bam files to analyse. only one bam per line and nothing more.

=back

=head2 optional arguments

=over 11

=item B<[-o]> ([output_file]) [./results.csv]

The output primary csv file containing the results of the analyses will be stored here.

=item B<[-n]> ([bases_to_get]) [50]

The number of base to analyse.

=item B<[-l]> ([interval_size]) [500]

Size of the intervals that will be analized

=item B<[-g]> ([output_plot]) [False]

Save of plot showing the differences between the different analized parts of the sequences

=back

=cut

my $man=0;
my $help=0;
my $inputDepth;
my $inputBed;
my $inputFiles;
my $toGet=50;
my $outputCSV="results.csv";
my $outputPlot='';
my $intervalSize=500;
my $verbose;

my $scriptPath=$0;
$scriptPath =~ s/[\w\.]*$//i; 
if ($scriptPath eq ""){
  $scriptPath=".";
}

GetOptions(
  "help|?"      => \$help,
  "man"         => \$man,
  "input_bed|b=s"   => \$inputBed,
  "input_bam_list|f=s"   => \$inputFiles,
  "output|o=s"  => \$outputCSV,
  "bases_to_get|n=i" => \$toGet,
  "interval_size|l=i" => \$intervalSize,
  "outputPlot|g!" => \$outputPlot,
  "verbose|v!" => \$outputPlot
) or pod2usage(2);

pod2usage(-verbose => 2) if ($man);
pod2usage(0) if $help;

pod2usage(
  -message => "Mandatory argument input_bed or input_bam_list is missing or is not a file",
  -verbose => 1,
) if ((!defined $inputBed || !-f $inputBed) || (!defined $inputFiles || !-f $inputFiles));

pod2usage(
  -message => "The number of base to lood for must be an integer",
  -verbose => 1,
) if (! $toGet  =~ /^[0-9,.E]+$/ );

pod2usage(
  -message => "Both autoGraph_begMidEnd.R and fonctionWIP.R must be present in the main script directory",
  -verbose => 1,
) if (!-f "$scriptPath/autoGraph_begMidEnd.R" || !-f "$scriptPath/fonctionWIP.R" );

pod2usage(
  -message => "The size of the interval to analyse must be two time greater than the number of bases to get",
  -verbose => 1,
) if ($intervalSize < 2*$toGet);

if(defined($outputCSV) && -e $outputCSV && !(-f $outputCSV))
{
  warn "$outputCSV already exists, but it's not a file";
  exit(0);
}


my @sum_begin;
my @sum_end;
my @sum_center;
my @sum_before;
my @sum_after;

my $line;
my $doubleToGet=$toGet*2;
my $tripleToGet=$toGet*3;
my @tabLine;
my $compteur=0;
my $intervalLength=0;
my $numInterval=-1; #Interval count
my $lineInterval;
my $sum=0;
my $compteurSup=0;
my $positionBeg=0;
my $positionEnd=0;
my $positionMid=0;
my $positionBefore=0;
my $positionAfter=0;
my $middle;
my $middleLess=0;
my $middleMore=0;
my @count;

#FIXME : samtools depth reconstruit les intervalles sur le bed étendu, donc les données sont faussées. Prendre ça en compte, les intervalles font 2 toget de moins. On a aussi un soucis dans la mesure où on obtient de nouveaux chevauchements en ajoutant des plages à analyser (on perd 1000 intervalles...)
#Utiliser le bed : pour CHAQUE position du depth, on regarde si cette position est dans notre intervale. int[1][0,1,2], si c'est pas le cas on passe à l'int suivant.
##########################################
# Preliminary analysis
#
# here we build a bed file using the data of samtools depth. The different steps are indicated in the different print

#We first expand a bed file, next we use samtools depth on it and next we rebuild a bed file from the output of samtools depth in order to merge overlapping intervals.

#TODO : ajouter des conditions pour pas que l'analyse tourne si les fichiers ont déjà été créés. Ajouter un élément distinctif dans le nom avec les paramètres n et l.
print "Expanding extremities of the bed file by $toGet...\n\n";
open (my $myBed, '<', $inputBed);
open (my $outputBed, '>', "modified_bed.bed");

#Building an extanded bed. We expand the extremities of each intervals from $toGet
while (<$myBed>){
  $line=$_;
  chomp($line);
  @tabLine = split /\s/, $line;
  $tabLine[1]=$tabLine[1]-$toGet;
  $tabLine[2]=$tabLine[2]+$toGet;
  print $outputBed "$tabLine[0]\t$tabLine[1]\t$tabLine[2]\n";
}

print "Running samtools depth on the specified sequences...\n\n";
system "samtools depth -aa -b modified_bed.bed -f $inputFiles > depth_data.csv";

print "Getting intervals from samtools depth output...\n\n";
get_interval_from_samtools_depth("depth_data.csv", "interval.bed");

print "Storing Bed in Array...\n\n";
my @bed_array = store_bed("interval.bed");


open (my $depthFile, '<', "depth_data.csv");
open (my $output, '>', $outputCSV);


for (my $i=0;$i<=5;$i++){ #Array Initialisation
  $count[$i]=0;
  for (my $j=0;$j<=$toGet;$j++){
    $sum_begin[$i][$j]=0;
    $sum_end[$i][$j]=0;
    $sum_center[$i][$j]=0; 
    $sum_before[$i][$j]=0; 
    $sum_after[$i][$j]=0; 
  }
}

#Moche : faire un passage sur le bed étendu uniquement pour chopper les extrémités qui vont en utilisant les données du bed normal. Ensuite on fait un passage classique sur le bed de base.

##########################################
# Processing the file 

print "Processing...\n";

while (<$depthFile>){ #"Compteur" will be decremented every time, it's our position on the sequence. "position" is the relative position in the interval of interest.
  if ($.%500000 == 0){ #Only to show the progression of the analyses.
    print "$. positions...\n";
  }
  $line = $_;
  chomp($line);
  @tabLine = split /\s/, $line;
  if ($compteur == 0){ #We go here when we have a new interval, so we need to reset some variable.
    positionBeg=0;
    $positionEnd=0;
    $positionMid=0;
    $positionBefore=0;
    $positionAfter=0;
    $numInterval++;
    $compteur = $bed_array[$numInterval][2]-$bed_array[$numInterval][1]; #These two values are the extremity of our sequence.
    $middle = $compteur/2; #FIXME : Soucis si c'est pas un entier ?
    $intervalLength=$compteur; 
    if ($intervalLength-$doubleToGet < $intervalSize*5){
      $count[int(($intervalLength-$doubleToGet)/$intervalSize)]++;
    }
    else{
      $count[5]++;
    }
    if ($compteur > $doubleToGet){ #Tester en prenant uniquement les intervalles > xxx
      #print "$numInterval - $intervalLength\n";
      $compteurSup = $compteur - $toGet;
      $middleLess=$middle-($toGet/2);
      $middleMore=$middle+($toGet/2);
    }
    next; #Because we dont wan't compteur--
  }
  if ($intervalLength > $doubleToGet){ #If the sequence is smaller than the bases we try to analyse, we will get some problems FIXME the if statement above is useless
    if ($compteur <= $toGet){ #after#
      for (my $i=2;$i<scalar(@tabLine);$i++){ #$sum will contain the sum of all the quality score for a gived position
        $sum+=$tabLine[$i];
      }
      $sum=$sum/(scalar(@tabLine)-2);
      if ($intervalLength-$doubleToGet < $intervalSize*5){
        $sum_after[int(($intervalLength-$doubleToGet)/$intervalSize)][$positionAfter]+=$sum;
      }
      else {
        $sum_after[5][$positionAfter]+=$sum;
      }
      $positionAfter++;
      $sum=0;
    }
    elsif ($compteur > $toGet && $compteur <= $doubleToGet){ #End#
      for (my $i=2;$i<scalar(@tabLine);$i++){ #$sum will contain the sum of all the quality score for a gived position
        $sum+=$tabLine[$i];
      }
      $sum=$sum/(scalar(@tabLine)-2);
      if ($intervalLength-$doubleToGet < $intervalSize*5){
        $sum_end[int(($intervalLength-$doubleToGet)/$intervalSize)][$positionEnd]+=$sum;
      }
      else {
        $sum_end[5][$positionEnd]+=$sum;
      }
      $positionEnd++;
      $sum=0;
    }
    elsif ($compteur < $compteurSup && $compteur >= $intervalLength-$doubleToGet){ #Begin#
      for (my $i=2;$i<scalar(@tabLine);$i++){
        $sum+=$tabLine[$i];
      }
      $sum=$sum/(scalar(@tabLine)-2);
      if ($intervalLength-$doubleToGet < $intervalSize*5){
        $sum_begin[int(($intervalLength-$doubleToGet)/$intervalSize)][$positionBeg]+=$sum;
      }
      else{
        $sum_begin[5][$positionBeg]+=$sum;
      }
      $positionBeg++;
      $sum=0;
    }
    elsif ($compteur > $compteurSup){ #before#
      for (my $i=2;$i<scalar(@tabLine);$i++){
        $sum+=$tabLine[$i];
      }
      $sum=$sum/(scalar(@tabLine)-2);
      if ($intervalLength-$doubleToGet < $intervalSize*5){
        $sum_before[int(($intervalLength-$doubleToGet)/$intervalSize)][$positionBefore]+=$sum;
      }
      else{
        $sum_before[5][$positionBefore]+=$sum;
      }
      $positionBefore++;
      $sum=0;
    }
    elsif ($compteur > $middleLess && $compteur <= $middleMore) { # If we are in a specific interval located in the middle of the sequence
      for (my $i=2;$i<scalar(@tabLine);$i++){
        $sum+=$tabLine[$i];
      }
      $sum=$sum/(scalar(@tabLine)-2);
      if ($intervalLength-$doubleToGet < $intervalSize*5){
        $sum_center[int(($intervalLength-$doubleToGet)/$intervalSize)][$positionMid]+=$sum;
      }
      else{
        $sum_center[5][$positionMid]+=$sum;
      }
      $positionMid++;
      $sum=0;
    }
  }
  $compteur--; #A position has been analysed
}

##########################################
# Calculating the mean for each position


for (my $i=0;$i<=5;$i++){ #Mean for each position.
  if ($count[$i]!=0)
  {
    print "$i : $count[$i]\n";
    for (my $j=0;$j<$toGet;$j++){
      $sum_center[$i][$j]=$sum_center[$i][$j]/$count[$i];
      $sum_begin[$i][$j]=$sum_begin[$i][$j]/$count[$i];
      $sum_end[$i][$j]=$sum_end[$i][$j]/$count[$i];
      $sum_before[$i][$j]=$sum_before[$i][$j]/$count[$i];
      $sum_after[$i][$j]=$sum_after[$i][$j]/$count[$i];
    }
  }
}

##########################################
# Output the results

for (my $i=0;$i<5;$i++){
  print $output "int ".$intervalSize*$i."-".$intervalSize*($i+1)."\t";
}
print $output "int ".$intervalSize*5 ."-inf\n ";

for (my $j=0;$j<$toGet;$j++){
  printf $output ("%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n", $sum_before[0][$j], $sum_before[1][$j], $sum_before[2][$j], $sum_before[3][$j], $sum_before[4][$j], $sum_before[5][$j]);
}

for (my $i=0;$i<10;$i++){
  print $output "NA\tNA\tNA\tNA\tNA\tNA\n";
}

for (my $j=0;$j<$toGet;$j++){
  printf $output ("%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n", $sum_begin[0][$j], $sum_begin[1][$j], $sum_begin[2][$j], $sum_begin[3][$j], $sum_begin[4][$j], $sum_begin[5][$j]);
}

for (my $i=0;$i<10;$i++){
  print $output "NA\tNA\tNA\tNA\tNA\tNA\n";
}


for (my $j=0;$j<$toGet;$j++){
  printf $output ("%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n", $sum_center[0][$j], $sum_center[1][$j], $sum_center[2][$j], $sum_center[3][$j], $sum_center[4][$j], $sum_center[5][$j]);
}

for (my $i=0;$i<10;$i++){
  print $output "NA\tNA\tNA\tNA\tNA\tNA\n";
}

for (my $j=0;$j<$toGet;$j++){
  printf $output ("%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n", $sum_end[0][$j], $sum_end[1][$j], $sum_end[2][$j], $sum_end[3][$j], $sum_end[4][$j], $sum_end[5][$j]);
}

for (my $i=0;$i<10;$i++){
  print $output "NA\tNA\tNA\tNA\tNA\tNA\n";
}

for (my $j=0;$j<$toGet;$j++){
  printf $output ("%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n", $sum_after[0][$j], $sum_after[1][$j], $sum_after[2][$j], $sum_after[3][$j], $sum_after[4][$j], $sum_after[5][$j]);
}

##########################################
# Output a plot (only if specified)

if ($outputPlot){
  my $RscriptPath = $scriptPath."autoGraph_begMidEnd.R";
  my $RfunctionGraphePath = $scriptPath."fonctionWIP.R";

  print "\nCreating plot...\n\n";

  system "Rscript $RscriptPath $outputCSV base_quality_positions.pdf $RfunctionGraphePath";

  print "The plot has been saved as \'base_quality_positions.pdf\'\n\n";
}

##########################################
# Print some informations if verbose is specified

if ($verbose){
  for (my $i=0;$i<=5;$i++){
    print "first interval : $count[$i] sequences\n";
  }
}

print "\nDone.\n\n";


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

##########################################
# END

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

29.O5.2017

=cut
