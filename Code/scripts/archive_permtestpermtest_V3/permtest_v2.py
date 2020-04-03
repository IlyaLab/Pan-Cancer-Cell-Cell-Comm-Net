

import sys
import gzip
import random
import numpy as np
import random

from multiprocessing import Pool

def samplevalues(f, p, mode):
    fin = gzip.open(f, 'rt')
    vals = []
    cnts = 0
    for line in fin:
        if cnts == 1000000:
            print(str(len(vals)))
            cnts = 0
        else:
            cnts += 1
        if random.random() < p:  ## p% of positive values
            bits = line.strip().split('\t')
            if bits[8] != 'NA':
                if mode == 'pos':
                    x = float(bits[8])
                    if x > 0.0:
                        vals.append(bits[8])  # randomly selected, and positive, and not NA
                else:
                    vals.append(float(bits[8]))  # if randomly selected and not an NA
    fin.close()
    return(vals)


def makeScores(pi, total, reps, vals):

    scores = []
    n = len(vals)
    s = set(range(0,n))
    for i in range(0,reps):
        idx = random.sample(s, pi)
        jdx = random.sample(s.difference(set(idx)), (total-pi))
        subset1 = vals[idx]
        subset2 = vals[jdx]
        m1 = np.median(subset1)
        m2 = np.median(subset2)
        a1 = np.median([abs(m1 - x) for x in subset1])
        b1 = np.median([abs(m2 - x) for x in subset2])
        scores.append( (m2-m1) / np.sqrt( pow(a1,2)+pow(b1,2) ))
    return(scores)


def makeScores2(n1, n2, reps, vals):

    scores = []
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
            if random.random() < 0.001:
                print(val)
        else:
            print('zero in divisor error')

    return(scores)



def main():
    print("*****")
    filein = sys.argv[1]
    sampin = sys.argv[2]
    fileout = sys.argv[3]

    fin = open(filein, 'r')
    sin = open(sampin, 'r')
    fout = open(fileout,'w')

    print('working on: ' + filein)
    print('with      : ' + sampin)
    print('writing to: ' + str(fout))

    fin.readline()  # skip header
    lines = fin.read().strip().split('\n')
    

    vals = np.array([ float(x) for x in vals if x != 'NA' and x != ' ' and x != ''])
    print('sampled ' + str(len(vals)) + ' values')

    #for pi in pts:
    #    print("    " + str(pi))
    #    scores = makeScores2(2204, 3720, 100000, vals)

    scores = makeScores2(4267.0, 4324.0, 200, vals)
    for si in scores:
        fout.write(str(si) + '\n')
    fout.close()

    print('done')
    print('max score: ' + str(max(scores)))
    print('min score: ' + str(min(scores)))
    print('median score: ' + str(np.median(scores)))
    print('*****')

if __name__ == '__main__':
    main()
