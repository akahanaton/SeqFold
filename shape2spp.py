import re
import sys

if len(sys.argv) < 2: 
    print 'python shape2spp.py shape_file outfile_prefix'
    sys.exit(1)

shape_filename = sys.argv[1]
shape_file = open(shape_filename, 'rU')
shape_data = shape_file.readlines()
shape_file.close()
    
outfile_prefix = sys.argv[2]

spp_filename = outfile_prefix+'.spp'
spp_file = open(spp_filename, 'w')

pref = []
for line in shape_data:
    field = re.split(' ', line[:-1])
    if field[1] == 'NA':
        pref.append('NA')
    else:
        pref.append(float(field[1]))
print len(pref)

pref_site = []
for i in range(0, len(pref)):
    if pref[i] != 'NA' and pref[i] > -999:
        pref_site.append(pref[i])
print len(pref_site)

pref_site.sort(reverse=True)
idx1 = int(len(pref_site)*0.02)
idx2 = int(len(pref_site)*0.1)
scale = sum(pref_site[idx1:idx2])/(idx2-idx1)
print str(idx1)+' '+str(idx2)+' '+str(scale)

for i in range(0, len(pref)):
    if pref[i] != 'NA':
        if pref[i] <= -999:
            pref[i] = 'NA'
        else:
            pref[i] = round(pref[i]/scale, 4)
            if pref[i] > 1:
                pref[i] = 1
            elif pref[i] < 0:
                pref[i] = 0
   
pref = map(str, pref)
pref_line = ';'.join(pref)
spp_file.write(outfile_prefix+'\t'+pref_line+'\n')
spp_file.close()        
