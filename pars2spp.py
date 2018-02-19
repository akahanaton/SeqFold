import re
import os
import operator
import math
import fisher

import sys
if len(sys.argv) < 3: 
    print 'python pars2spp.py s1_file v1_file outfile_prefix'
    sys.exit(1)

S1_file_name = sys.argv[1]
V1_file_name = sys.argv[2]
outfile_prefix = sys.argv[3]

cutoff_sp = 0.05

test_file_name = outfile_prefix+'.test'
spp_file_name = outfile_prefix+'.spp'
        
def parsTest(S1_file_name, V1_file_name, test_file_name):
    S1_file = open(S1_file_name,'r')
    S1_data = S1_file.readlines()
    S1_file.close()
    V1_file = open(V1_file_name, 'r')
    V1_data = V1_file.readlines()
    V1_file.close()
    S1_score = {}
    S1T = 0
    for line in S1_data:
        line_item = re.split('\t', line[:-1])
        name = line_item[0]
        score = line_item[1].split(';')
        score = map(int, score)
        S1_score[name] = score
        S1T = S1T+sum(score)
    V1_score = {}
    V1T = 0
    for line in V1_data:
        line_item = re.split('\t', line[:-1])
        name = line_item[0]
        score = line_item[1].split(';')
        score = map(int, score)
        V1_score[name] = score
        V1T = V1T+sum(score)
    logR = {}
    pval_S1_v = []
    pval_V1_v = []    
    sk = float((V1T+S1T))/float(2*S1T)
    vk = float((V1T+S1T))/float(2*V1T)
    genes = []
    for name in S1_score:
        if not (name in V1_score):
            continue
        genes.append(name)
        S1 = S1_score[name]
        V1 = V1_score[name]
        logR[name] = []
        for i in range(0, len(S1)):
            pval = fisher.pvalue(S1[i], (S1T-S1[i]), V1[i], (V1T-V1[i]))
            logR[name].append(round(math.log((sk*S1[i]+1)/(vk*V1[i]+1), 2), 4))
            pval_S1_v.append(pval.right_tail)
            pval_V1_v.append(pval.left_tail)
    fdr_S1_v = fdr_bh(pval_S1_v)
    fdr_V1_v = fdr_bh(pval_V1_v)
    
    test_file = open(test_file_name, 'w')
    test_file.write('#'+str(S1T)+'\n')
    test_file.write('#'+str(V1T)+'\n')
    j = 0
    jj = 0
    for j in range(0, len(genes)):
        name = genes[j]
        test_file.write('>'+name+'\n')
        for i in range(0, len(S1_score[name])):
            test_file.write(str(i+1)+'\t'+str(S1_score[name][i])+'\t'+str(V1_score[name][i])+'\t'+str(logR[name][i])+'\t'+str(fdr_S1_v[jj])+'\t'+str(fdr_V1_v[jj])+'\n')
            jj = jj+1
        j = j +1
    test_file.close()


def fdr_bh(pv):
    if not pv:
        return []
    m = len(pv)
    args, pv = zip(*sorted(enumerate(pv), None, operator.itemgetter(1)))
    if pv[0] < 0 or pv[-1] > 1:
        raise ValueError("p-values must be between 0 and 1")
    qvalues = m * [0]
    mincoeff = pv[-1]
    qvalues[args[-1]] = mincoeff
    for j in xrange(m-2, -1, -1):
        coeff = m*pv[j]/float(j+1)
        if coeff < mincoeff:
            mincoeff = coeff
        qvalues[args[j]] = mincoeff
    return qvalues


def getSPP(test_file_name, spp_file_name, cutoff_sp):
    test_file = open(test_file_name,'r')
    test_data = test_file.readlines()
    test_file.close()
    spp_file = open(spp_file_name, 'w')
    i = 0
    count_total = 0
    count_unpair = 0
    count_pair = 0
    for line in test_data:
        line_item = re.split('\t', line[:-1])
        if line.startswith('#'):
            continue
        elif line.startswith('>'):
            if i > 0:
                spp = ';'.join(spp)
                spp_file.write(name+'\t'+spp+'\n')
            name = line_item[0][1:]
            spp = []
            i = i+1
        else:
            count_total = count_total+1
            if (float(line_item[4]) < cutoff_sp):
                spp.append('1')
                count_unpair = count_unpair+1
            elif (float(line_item[5]) < cutoff_sp):
                spp.append('0')
                count_pair = count_pair+1
            else:
                spp.append('NA')
    spp = ';'.join(spp)
    spp_file.write(name+'\t'+spp+'\n')
    spp_file.close()
    print 'Total bases: '+str(count_total)
    print 'Unpaired bases: '+str(count_unpair)
    print 'Paired bases: '+str(count_pair)


parsTest(S1_file_name, V1_file_name, test_file_name)
getSPP(test_file_name, spp_file_name, cutoff_sp)
os.system('rm '+test_file_name)
    