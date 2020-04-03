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
def createX(fileMap,idx,fileListPath):
	# fileMap takes an index number and returns a file name
	# idx is the index to EB++ data column
	# For files containing samples:
		# Need to read a line from each of the files
		# Subset that set of lines to match idx
		# Each result file has column header with sample IDXs
		# get first row of each file we need
		# subset features.
	#minSample = np.min(idx)
	#maxSample = np.max(idx)
	#fs = list(set( [fileMap[i] for i in range(minSample,maxSample+1))  # get the file names to open
	#fins = [open(fi,'r') for fi in fs]
	fans = open(fileListPath,'r').read().strip().split()
 	fins = [  open(fi,'r') for fi in fans[0]  ] 
	firstRows = [fi.readline().strip().split('\t') for fi in fins]   # this is the sample indices
	X = []
	# for each list of samples in input files
	#    if that sample is in the idx
	#		then save it's index (location) *within* the file. 
	fidx = [  [j for j,k in enumerate(jdxs) if k in idx] for jdxs in firstRows ]  # for each list of samples, if in idx, save it.
	# should have a list of lists, each is index *within* resPar files.
	for (i in range(0,100):
		rows = [fi.readline().strip().split('\t') for fi in fins] # get next row
		vals = [float(x) for x in rows]
		nums = [np.array(x) for x in rows]
		keep = [ns[fi] for (ns,fi) in zip(nums,fidx)]  # should have list of np arrays...
		valx = np.concatenate(keep)  
		# Need to put valx into correct order.
		# if valx passes some test
		X.append(valx)
	return(X)



# do CV classification

def main():
	if (len(sys.argv) < 2):
		print('resultFiles phenoFile') 
		return(0)
	print('running')
	resultFileListPath = sys.argv[1]
	phenoFile          = sys.argv[2]

	pin   = open(phenoFile,'r').read().strip().split('\n')
	idx   = [int(x) for x in pin[0].strip().split('\t')]
	pheno = [int(x) for x in pin[1].strip().split('\t')]

	print('idx')
	print(idx[0:5])

	print('pheno')
	print(pheno[0:5])

	fileMap = openFiles(resultFileListPath)        # maps sample idx to file name
	fils2rd = list(set([fileMap[i] for i in idx])) # get list of files to read
	
	X = createX(fils2rd,idx)

	clf = svm.SVC(kernel='linear', C=1)
	scores = cross_val_score(clf, X, pheno, cv=5)
	print(scores)
	return(scores)



if __name__ == "__main__":
	main()




