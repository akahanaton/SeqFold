#!/usr/bin/perl -w

########################################################################################
#mondify the value of $para in a parallel computing enviornment to speed up, e.g.,
#my $para = "qsub -cwd -V -l h_vmem=3G -l h_rt=6:00:00 -m ea -w e -b y";
#my $para = "bsub -M 3072000 -W 6:00";  
########################################################################################


use strict;
use warnings;
use Cwd;


my $para = ""; 
  

if (@ARGV < 3) {
	print "usage: perl sfold_wrapper.pl sfold_executable_file input_fasta_file sfold_output_directory\n";
	print "this command runs Sfold for a fasta file containing one or more sequences\n";
	print "\tsfold_executable_file     the path of Sfold executable file\n";
	print "\tinput_fasta_file          the sequence file in FASTA format\n";
	print "\tsfold_output_directory    Sfold output directory\n";
	exit(0);
}

my $Sfold = $ARGV[0];
my $inputFasta = $ARGV[1];
my $outDir = $ARGV[2];


#test if Sfold binary exist
if (! (-e $Sfold)) {
	print "ERROR: Sfold executable file does not exist ($Sfold)\n";
	exit(-1);
}

#test if input fasta file exist
if (! (-e $inputFasta)) {
	print "ERROR: input fasta file does not exist ($inputFasta)\n";
	exit(-1);
}

#create output directory if not exist
if (! (-e $outDir)) {
    system ("mkdir " . $outDir);
}

#print parameter
print "reading:          $inputFasta\n";


#####################################################################
#open input and output file
open(inputFastaFile, $inputFasta);


#read all sequences
my $i = 0;
while (my $line1 = <inputFastaFile>) {
	chomp($line1); 					
	my $line2 = <inputFastaFile>;   
	chomp($line2);
	$i++;

	my $seq = uc $line2;		    
	$seq =~ tr/U/T/;				
	my $len = length($seq);
	$line1 =~ />(.*)/;				
	my $seqName = $1;
	print "$seqName\n";
    	
	my $name = $outDir . "/" . $seqName . ".fa";
    open(file_local, ">$name");
    print file_local ">$seqName\n";
    print file_local "$seq\n";
    close(file_local);
    
    #Sfold sampling and clustering for each sequence
	my $cmd = $para . " " . $Sfold . " -a 1 -o " . $outDir . "/$seqName $name";
    system ($cmd);
}
