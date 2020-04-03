

import sys
import gzip
import random
import numpy as np
import random

from multiprocessing import Pool

def samplevalues(f):
    fin = open(f, 'rt')
    vals = []
    for line in fin:
        try:
            val = float(line.strip())
            vals.append(val)
        except:
            print('error: ' + line)  # if randomly selected and not an NA
    fin.close()
    return(np.array(vals))


def makeScores( paramList ):
    n1, n2, reps, vals = paramList
    scores = []
    med_diff = []
    n = len(vals)
    s = set(range(0,n))
    for i in range(0,reps):
        idx = random.sample(s, int(n1))
        jdx = random.sample(s.difference(set(idx)), int(n2))
        subset1 = vals[idx]
        subset2 = vals[jdx]
        m1 = np.median(subset1)
        m2 = np.median(subset2)
        a1 = np.median([abs(m1 - x) for x in subset1])
        b1 = np.median([abs(m2 - x) for x in subset2])
        if a1 > 0 or b1 > 0:
            val = (m2-m1) / np.sqrt( 1.4826*a1 + 1.4826*b1 )
            scores.append(val)
        else:
            print('zero in divisor error')

    return(scores)



def main():
    print("*****")
    filein = sys.argv[1]  # has the sample size for n1 and n2 per study / subtype
    fileout = sys.argv[2]

    fin = open(filein, 'r')
    fout = open(fileout,'w')

    print('working on: ' + filein)
    print('writing to: ' + str(fout))

    fin.readline()  # skip header
    lines = fin.read().strip().split('\n')

    p = Pool(7)

    for li in lines:
        bits = li.strip().split(',')
        filename = 'permSamples2/samples_'+bits[0]+'.csv'
        vals = samplevalues(filename)
        n1 = float(bits[1])
        n2 = float(bits[2])
        print(filename + '  ' + str(n1) + '   ' + str(n2))
        params = [(n1,n2, 1000, vals) for i in range(0,1000)]
        scores = p.map(makeScores, params) #(n1, n2, 1000000, vals)
        for i in range(0,len(scores)):
            for j in scores[i]:
                fout.write(bits[0] + ',' + bits[1] + ',' + bits[2] + ',' + str(j) + '\n')

    fout.close()


# goal is to get a distribution of scores if the phenotype was random
# with sample sizes that match each test using samples from those studies.


if __name__ == '__main__':
    main()
