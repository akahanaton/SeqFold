SeqFold: Genome-scale reconstruction of RNA secondary structure integrating experimental data

A tool for RNA secondary structure prediction from experimental data. Given high-throughput data, it can reconstruct the secondary structures of the whole transcriptome, termed the RNA structurome. It outputs reconstructed RNA secondary structure as well as base-level accessibilities for each transcript. Input data types include parallel analysis of RNA structure (PARS), selective 2’-hydroxyl acylation analyzed by primer extension sequencing (SHAPE-Seq), fragmentation sequencing (FragSeq), and conventional SHAPE data.

For updates, please refer to http://www.stanford.edu/~zouyang/seqfold


-------------
Citation
-------------

Ouyang Z, Snyder MP, and Chang HY (2012) SeqFold: Genome-scale reconstruction of RNA secondary structure integrating high-throughput sequencing data. Genome Research.


-------------
License
-------------

Use of SeqFold is free for non-commercial research. Commercial users please contact the authors.


-------------
Prerequisites
-------------

1. Prepare data files from one of the following experiments.

PARS:
A tab-delimited file of read counts for RNase S1, and another for V1. One transcript per row. Format: transcript name (column 1), counts (column 2, semicolon separated, the raw number of reads obtained for each base).

SHAPE(-Seq):
A space-delimited file of SHAPE reactivities for each base of a transcript. One base per row. Format: base number (column 1), SHAPE reactivity (column 2).

FragSeq:
A folder with cutting score files (*.cutscores.ss.list) output by FragSeq_v0.0.1.

2. Prepare a sequence file in FASTA format.

3. Download and install Sfold following its description. Sfold can be downloaded here:
http://sfold.wadsworth.org/cgi-bin/index.pl


-------------
Install
-------------

1. Decompress the downloaded SeqFold package:
tar -zxvf seqfold.tar.gz

2. Change to the extracted directory:
cd seqfold


-------------
Run SeqFold
-------------

1. Generate structure preference profiles

PARS:
python pars2spp.py s1_file v1_file outfile_prefix

SHAPE(-Seq):
python shape2spp.py shape_file outfile_prefix

FragSeq:
python fragseq2spp.py path_to_cuttingscore_folder outfile_prefix

Each command will output a file with structure preference profiles. One transcript per row. Format: transcript name (column 1), structure preferences (column 2, semicolon separated, the preference of single-strandness for each base).

2. Generate sample structures and clusters for each transcript

perl sfold_wrapper.pl sfold_executable_file input_fasta_file sfold_output_directory

Note: In a parallel computing enviornment, one can speed up the runtime by mondifying the value of $para in sfold_wrapper.pl. Example: $para = "bsub -M 3072000 -W 6:00";

3. Generate RNA secondary structure predictions and base-level accessibilities

python seqfold.py sfold_output_directory structure_preference_profile

Optional parameters
-d    Output directory. Default: ./ 
-o    Prefix of output summary files. Default: out 
-f    Cutoff for sequences to be filtered with <= cutoff_frac fraction of sites having experimental data. Default: 0


-------------
Output Files
-------------

*.seqfold.ct
The predicted secondary structure of each transcript in CT format. One transcript a file. 

out.acc
The accessibility of each base of each transcript. One transcript per row. Format: transcript name (column 1), accessibilities (column 2, semicolon separated, the accessibility of each base) 


============================================================
last modified on October 18, 2012
