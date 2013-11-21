'''
Reads the reactome human protein interactions tsv, then prints the ensemble
reaction pairs that have both ends within the ensemble ids passed in as a 
list in the first arg.
-stlee
'''
import sys


# read conversion tables:
inf = open("probe2ensembl.txt", 'r')
ensembl2p = {}

for line in inf:
    d = line.split()
    if len(d) < 3:
        continue
    ensembl2p[int(d[2][5:-1])] = d[1]
inf.close()


inf = open("probe2entrez.txt", 'r')
p2entrez = {}

for line in inf:
    d = line.split()
    if len(d) < 3:
        continue
    p2entrez[d[1]] = int(d[2][1:-1])
inf.close()

convert = {}
for k, v in ensembl2p.iteritems():
    if v in p2entrez:
        convert[k] = p2entrez[v]


#print convert

# read actual interactions:
inf = open("homo_sapiens.interactions.txt", 'r')

m = {}

for line in inf:
    d = line.split('\t')

    if len(d) < 5:
        continue

    if len(d[1]) < 1 or len(d[4]) < 1: # no ensembl tags...
        continue

    f = [eid[12:] for eid in d[1].split('|')]
    t = [eid[12:] for eid in d[4].split('|')]
    
    for f_eid in f:
        for t_eid in t:
            try:
                #m[(int(f_eid), int(t_eid))] = 1
                m[(f_eid, t_eid)] = 1
            except:
                continue

#print m
#print len(m)

inf.close()

# read gene list
inf = open(sys.argv[1], 'r')
s = set()
for line in inf:
    s.add(int(line.strip()))

inf.close()

# generate interaction pairs
outf = open(sys.argv[2], 'w')
i = 0
for k, v in m.keys():
    if k == v:
        continue

    try:
        ik = int(k)
        iv = int(v)
    
        ck = convert[ik]
        cv = convert[iv]
    except:
        #print i, "whoops"
        i += 1
        continue

    #if k in s and v in s and int(k) in convert and int(v) in convert:
    #    print outf.write(str(convert[k]) +'\t' + str(convert[v]) + '\n')
    #    i += 1

    if ck in s and cv in s:
        outf.write(str(ck) + '\t' + str(cv) + '\n')
        i += 1

#print i
