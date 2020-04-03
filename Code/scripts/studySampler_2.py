

#import gzip
import sys
import random 


# get the study for each barcode
tcga = open('PanCancer/PanCancer_Clinical_Extra_Endpoints.csv', 'r').read().strip().split('\n')

# get the number of samples in each study.
tcgaCount = dict()
studyList = []
sc = open('sampleCounts/studyCounts.csv','r').read().strip().split('\n')
for si in sc:
	print('si: ' + si)
	bits = si.split(',')
	tcgaCount[bits[0]] = float(bits[1])
	studyList.append(bits[0])

tcgaDict = dict()   # look up from barcodes to studies
for i in range(0,len(tcga)):
	bits = tcga[i].split(',')
	tcgaDict[bits[1]] = bits[2]


# open each output file
fileDict = dict()
for si in studyList:
	fileDict[si] = open('permSamples2/samples_'+si+'.csv', 'w')


# creates our dict of prob thresholds.
probT = dict()
for si in studyList:
	n = float(tcgaCount[si])
	probT[si] = 1.00 / (n*10.00)   # selects 10% of weights

dat = 'ready'
while dat != '':
	try:
		dat = sys.stdin.readline()
		bits = dat.split(',')
		study = tcgaDict[ bits[0][0:12] ]
		value = bits[3].strip()
		if random.random() < probT[study] and value != 'NA':
			fileDict[study].write(value + '\n')
	except:
		print('exception: ' + dat)

for si in studyList:
	fileDict[si].close()

print('done!')
