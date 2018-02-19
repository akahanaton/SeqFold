import re
import os
import operator
import glob

import sys
if len(sys.argv) < 2: 
    print 'python fragseq2spp.py path_to_cuttingscore_folder outfile_prefix'
    sys.exit(1)

cutoff = 0.6931472

path = sys.argv[1]
outfile_prefix = sys.argv[2]
cutoff = float(sys.argv[3])

spp_filename = outfile_prefix+'.spp'
spp_file = open(spp_filename, 'w')

i = 0
for infile_name in glob.glob( os.path.join(path, '*.cutscores.ss.list') ):
    print infile_name
    infile = open(infile_name,'r')
    infile_data = infile.readlines()
    infile.close()
    
    idx1 = infile_name.rindex('/')+1
    name = infile_name[idx1:]
    idx2 = name.index('.')
    name = name[:idx2]            
            
    pref = []            
    j = 0
    for line in infile_data:
        j = j+1
        if j == 1:
            continue
        field = re.split(' ', line[:-1])
        if field[0] != 'none' and float(field[0]) >= cutoff:
            pref.append('1')
        else:
            pref.append('NA')            
    pref = ';'.join(pref)
    spp_file.write(name+'\t'+pref+'\n')
spp_file.close()
        