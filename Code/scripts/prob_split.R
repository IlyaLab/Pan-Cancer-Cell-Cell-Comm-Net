

#doppio(~/cytokine_network)% ls /titan/cancerregulome9/workspaces/users/dgibbs/CytokineNetworkData/Prob/ | head
#
#prob_aDC_1.rda
#prob_aDC_2.rda
#prob_aDC_3.rda
#prob_aDC_4.rda
#prob_Adipocytes_1.rda
#prob_Adipocytes_2.rda
#prob_Adipocytes_3.rda
#prob_Adipocytes_4.rda
#prob_astrocytes_1.rda

load('Data/p_distr_v2.rda')

cellNames <- names(condECDF)

names(condECDF[[cellNames[1]]])
#[1] "0%_25%"   "25%_50%"  "50%_75%"  "75%_100%"

for (ci in cellNames) {
	
	for (di in 1:4) {
		
	   x <- condECDF[[ci]][[di]]

	   save(x, file= paste0('/titan/cancerregulome9/workspaces/users/dgibbs/CytokineNetworkData/Prob_v2/prob_', ci ,'_', di ,'.rda')) 
	}

}


