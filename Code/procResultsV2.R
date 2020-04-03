

# processing the results

load('accessories_v2.rda')

files = list.files('Results_v2/')

save(files, file='result_files_v2.rda')

for (i in 1:length(files)) {
	load( paste0('Results_v2/',files[i]) )
	print(paste0('working on: ', files[i]) )
	if (length(res1$Vals) == nrow(ccs)) {
		# then we have matching vals and edges
		df <- data.frame(Barcode = res1$Barcode, SampleID=ti, EdgeID = 1:1062719, EdgeWt = round(res1$Vals, digits=10))
		fileout = paste0('Dataframes/table', i, '.csv')
		write.table(df, file=fileout, sep=',', quote=F, row.names=F, col.names=F)
		system(paste0("gzip ", fileout))
	} else {
		print ( paste0('error on ', i) )
	}
}


system('zcat Dataframes/*.csv.gz | gzip > all_results_v2.csv.gz')
