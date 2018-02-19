#################################################################################
#SeqFold - A tool for genome-scale reconstruction of RNA secondary structure integrating experimental measurements

#Zhengqing Ouyang, Stanford Univerisity School of Medicine, zouyang@stanford.edu

#Please cite: Ouyang Z, Snyder MP, and Chang HY (2012) SeqFold: Genome-scale reconstruction of RNA secondary structure integrating high-throughput sequencing data. Genome Research. 

#Optional parameters
#Output directory. Will create a new directory if one does not already exist
out_dir = './' 
#Prefix of output summary files
out_prefix = 'out' 
#Cutoff for sequences to be filtered with <= cutoff_frac fraction of sites having experimental data
cutoff_frac = 0 
#################################################################################


import getopt
import os
import re
import glob
import numpy
from scipy import array

import sys
if len(sys.argv) < 2: 
    print 'python seqfold.py sfold_output_directory structure_preference_profile -d [Output directory: default ./] -o [Prefix of output summary files: default out] -f [Cutoff for sequences to be filtered with <= cutoff_frac fraction of sites having experimental data]'
    sys.exit(1)


sfold_dir = sys.argv[1]
spp_file_name = sys.argv[2]

#Optional parameters
try:
    opts, args = getopt.getopt(sys.argv[3:], 'd:o:f:')
except getopt.GetoptError, err:
    print str(err)
    usage()
    sys.exit(2)
for o,a in opts:
    if o == "-d":
		out_dir = a
    if o == "-o":
		out_prefix = a
    if o=="-f":
		cutoff_frac = float(a)

def structureMatch(spp_file_name, sfold_dir, out_dir, out_prefix):    
    acc_file_name = out_dir+'/'+out_prefix+'.acc'
    err_file_name = out_dir+'/'+out_prefix+'.err'
    acc_file = open(acc_file_name,'w')
    err_file = open(err_file_name,'w')

    spp_file = open(spp_file_name,'r')
    spp_data = spp_file.readlines()
    spp_file.close()
                    
    for line in spp_data:
        line_item = re.split('\t', line[:-1])       
        symbol = line_item[0]
        print symbol 
        
        path = sfold_dir + symbol + '/clusters/'
        
        cluster_file = path + 'ch.index.out'
        if os.path.exists(cluster_file) == False:
            err_file.write(symbol+'\t'+'Sfold clusters not exist'+'\n') 
            continue
        elif os.stat(cluster_file).st_size == 0:
            err_file.write(symbol+'\t'+'Sfold cluster file empty'+'\n') 
            continue
        
        bpdist_file = path + 'c01.bp.dist.from.ccentroid.out'
        if os.path.exists(bpdist_file) == False:
            err_file.write(symbol+'\t'+'Sfold clustering error'+'\n') 
            continue
        elif os.stat(bpdist_file).st_size == 0:
            err_file.write(symbol+'\t'+'Sfold clustering error'+'\n') 
            continue
        
        spp = []
        spp = line_item[1].split(';')
        
        site = []
        spp_site = []
        for j in range(0, len(spp)):
            if spp[j] != 'NA':
                site.append(j)
                spp_site.append(float(spp[j]))
        if float(len(spp_site))/float(len(spp)) <= cutoff_frac:
            err_file.write(symbol+'\t'+'fraction of SPP sites <= '+str(cutoff_frac)+'\n') 
            continue
        elif max(spp_site) > 1 or min(spp_site) < 0:
            err_file.write(symbol+'\t'+'SPP not within [0,1]')
            continue
        
        sfold_file_name = sfold_dir + symbol + '/bp.out'
        
        sfold_file = open(sfold_file_name,'r')
        sfold_data = sfold_file.readlines()
        sfold_file.close()
        sfold_up = {}   
        
        i = 0
        breaker = False
        for line1 in sfold_data:
            if line1.find('Structure') != -1:
                i = i + 1
                sfold_up[str(i)] = []
                for j in range(0, len(spp)):
                    sfold_up[str(i)].append(1)    
            else:
                line1 = line1.lstrip()
                space = re.compile(r'\s+')
                line1_item = space.split(line1[:-1])
                if int(line1_item[0]) > len(spp) or int(line1_item[1]) > len(spp):
                    err_file.write(symbol+'\t'+'sample structure length larger than SPP length'+'\n')
                    breaker = True
                    break 
                sfold_up[str(i)][int(line1_item[0])-1] = 0
                sfold_up[str(i)][int(line1_item[1])-1] = 0
        
        if breaker:
            continue
        
        if i < 1000:
            err_file.write(symbol+'\t'+'sample structures < 1000'+'\n')   
            continue
                          
        dist = {}
        for name in sfold_up.keys():
            dist[name] = 0
            for j in range(0, len(spp_site)):
                dist[name] = dist[name]+abs(spp_site[j]-sfold_up[name][site[j]])
        
        mind = min(dist.values())
        mind_set = []
        for name in sfold_up.keys():
            if dist[name] == mind:
                mind_set.append(name)
                        
        clust_map = {}
        clust_count = {}
        for sfold_file_name in glob.glob( os.path.join(path, '*.list') ):
            idx1 = sfold_file_name.rindex('/')+1
            name = sfold_file_name[idx1:]
            idx2 = name.index('.')
            name = name[:idx2]
            clust_count[name] = 0
            
            sfold_file = open(sfold_file_name,'r')
            sfold_data = sfold_file.readlines()
            sfold_file.close()            
            for line1 in sfold_data:
                line1_item = re.split(' ', line1[:-1])
                for j in range(0, len(line1_item)):
                    clust_map[line1_item[j]] = {}
                    clust_map[line1_item[j]] = name
                    
        for i in range(0, len(mind_set)):
            name = clust_map[mind_set[i]]
            clust_count[name] += 1
        
        keys = sorted(clust_count)    
        keys.sort(key=clust_count.get, reverse=True)  
        
        cluster_name = keys[0]
        
        num = 0.0    
        up = numpy.zeros(len(spp), float)      
        for strid in sfold_up.keys():
            if clust_map[strid] == cluster_name:
                num = num+1.0
                up = up + numpy.array(sfold_up[strid])
        up = up/num
        up = up.round(decimals=4)
        acc = map(str, list(up))
        acc = ';'.join(acc)
        
        acc_file.write(symbol+'\t'+acc+'\n')
        
        centroid_file_name = path + cluster_name + '.ccentroid.ct'
        centroid_file = open(centroid_file_name,'r')
        centroid_data = centroid_file.readlines()
        centroid_file.close()
        single_ct_file_name = out_dir+'/'+symbol+'.seqfold.ct'
        single_ct_file = open(single_ct_file_name,'w')
        single_ct_file.write('>'+symbol+'\n') 
        for line in centroid_data:
            line = line.lstrip()
            if line.startswith('SFOLD'):
                continue  
            single_ct_file.write(line)
        single_ct_file.close()
        
    acc_file.close()
    err_file.close()


if os.path.exists(out_dir) == False:
    os.system('mkdir '+out_dir) 
    

print 'Inferring structure and accessibility...'
structureMatch(spp_file_name, sfold_dir, out_dir, out_prefix)
print 'Done.'
