workdir = '/titan/cancerregulome9/workspaces/users/dgibbs/CytokineNetworkData'

scriptList = ['scripts/permtest_v3_sample1_glados.py',
              'scripts/permtest_v3_sample2_glados.py',
              'scripts/permtest_v3_sample3_glados.py']


sampleList = ['sampleCounts/cn_n_study_pfi.csv',
'sampleCounts/cn_n_study_stage.csv',
'sampleCounts/cn_n_subtype_pfi.csv',
'sampleCounts/cn_n_subtype_stage.csv'
]

outfileList = ['nullStatsPool/cn_n_study_pfi_',
'nullStatsPool/cn_n_study_stage_',
'nullStatsPool/cn_n_subtype_pfi_',
'nullStatsPool/cn_n_subtype_stage_'
]


for fileCnt in range(1,101):
    for i in range(0,3):
        for j in range(0,4):
            cmd = ' '.join(['python3', scriptList[i], sampleList[j], outfileList[j]+'sample'+str(i+1)+'_'+str(fileCnt)+'.csv'])
            print(cmd)
