

#import gzip
import sys
import random 


# get the study for each barcode
tcga = open('PanCancer/five_signature_mclust_ensemble_results.tsv', 'r').read().strip().split('\n')

# get the number of samples in each study.s
subtypeCount = dict()
subtypeList = []
sc = open('sampleCounts/subtypeCounts.csv','r').read().strip().split('\n')
for si in sc:
	print('si: ' + si)
	bits = si.split(',')
	subtypeCount[bits[0]] = float(bits[1])
	subtypeList.append(bits[0])


tcgaDict = dict()   # look up from barcodes to studies
for i in range(0,len(tcga)):
	bits = tcga[i].split('\t')
	barcode = bits[0].replace('.', '-')
	tcgaDict[barcode] = bits[2]


# open each output file
fileDict = dict()
for si in subtypeList:
	fileDict[si] = open('permSamples3/samples_'+si+'.csv', 'w')


# creates our dict of prob thresholds.
probT = dict()
for si in subtypeList:
	n = float(subtypeCount[si])
	probT[si] = 1.00 / (n*10.00)   # selects 10% of weights


dat = 'ready'
while dat != '':
	try:
		dat = sys.stdin.readline()
		bits = dat.split(',')
		study = tcgaDict[ bits[0] ]
		value = bits[3].strip()
		if random.random() < probT[study] and value != 'NA':
			fileDict[study].write(value + '\n')
	except:
		print('exception: ' + dat)

for si in subtypeList:
	fileDict[si].close()

print('done!')
