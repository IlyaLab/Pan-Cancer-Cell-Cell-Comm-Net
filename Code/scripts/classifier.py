import sys
from sklearn.model_selection import cross_val_score

def openFiles(fileListPath):
	# given the list of results files with 
	# what samples are found in each file 
	# (which is on each row)
	# return a dictionary that takes a sample index
	# and returns what file to read.
	sampleLookup = dict()
	fans = open(fileListPath,'r').read().strip().split()
	files = fans[0]
	rangs = fans[1]	
	for i in range(0,len(files)):
		rsPts = eval(rangs[i])
		for j in range(rsPts[0], rsPts[1]):
			sampleLookup[j] = files[i]
	return(sampleLookup)


# create the data set
def createX(fileMap,idx):
	# fileMap takes an index number and returns a file name
	# idx is the index to EB++ data column
	return(0)


# subset features.

# do CV classification

def main():
	if (len(sys.argv) < 2):
		print('resultFiles phenoFile') 
		return(0)
	print('running')
	resultFileListPath = sys.argv[1]
	phenoFile          = sys.argv[2]

	pin = open(phenoFile,'r').read().strip().split('\n')
	idx = pin[0].strip().split('\t')
	pheno = pin[1].strip().split('\t')

	print('idx')
	print(idx[0:5])

	print('pheno')
	print(pheno[0:5])

	fileMap = openFiles(resultFileListPath)        # maps sample idx to file name
	fils2rd = list(set([fileMap[i] for i in idx])) # get list of files to read
	
	X = createX(fils2rd,idx)

	#clf = svm.SVC(kernel='linear', C=1)
	#scores = cross_val_score(clf, X, Y, cv=5)
	return(0)



if __name__ == "__main__":
	main()




