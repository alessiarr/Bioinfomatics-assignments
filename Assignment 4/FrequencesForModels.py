#Assignment 4, Bioinformatics
##get CpG Islands from chromosome 22 and relative frequences for inside and outside models

import pandas as pd
import random

def getFrequences(file, fastaFile):
	fasta = open(fastaFile, 'r')
	lines = fasta.readlines()
	#print(lines)
	chr22 = ''
	for line in lines:
		if line[0] != '<':
			chr22 += line[:-1]
	fasta.close 
    
	CpGIslands = pd.read_table(file)
	CpGChr22Start = CpGIslands.loc[lambda CpGIslands: CpGIslands['chr'] == 'chr22']['start'].reset_index(drop=True) #select the columns corresponding to 'start' and 'end' of chr22
	#print(chr22)
	CpGChr22End = CpGIslands.loc[lambda CpGIslands: CpGIslands['chr'] == 'chr22']['end'].reset_index(drop=True)
	
	relativeFreqInside = pd.DataFrame(columns=['A','C','G','T'],index=['A','C','G','T','tot'])
	relativeFreqInside.iloc[:] = 0
	relativeFreqOutside = pd.DataFrame(columns=['A','C','G','T'],index=['A','C','G','T','tot']) 
	relativeFreqOutside.iloc[:] = 0

	for i in range(len(CpGChr22Start)):
		length = CpGChr22End.iat[i] - CpGChr22Start.iat[i]
		island = chr22[CpGChr22Start.iat[i]:CpGChr22End.iat[i]+1].upper()
		#print(length)
		start = random.randint(0,len(chr22))
		nonIsland = chr22[start:start+length+1].upper()
		for j in range(1,length):
			relativeFreqInside.at[island[j-1],island[j]]+=1 #increase frequency of the dimer j-1/j (j-1 on the rows and j on the columns)
			if nonIsland[j] != 'N' and nonIsland[j-1] != 'N':
				relativeFreqOutside.at[nonIsland[j-1],nonIsland[j]]+=1
	for i in range(0,4):
		relativeFreqInside.iat[4,i] = sum([relativeFreqInside.iat[j,i] for j in range(0,4)])
		relativeFreqOutside.iat[4,i] = sum([relativeFreqOutside.iat[j,i] for j in range(0,4)])
    
	#print(relativeFreqInside,relativeFreqOutside)
	x = random.randint(0,len(CpGChr22Start))
	y = random.randint(0,len(chr22)) 
	CpGExample = chr22[CpGChr22Start[x]:CpGChr22End[x]+1] #will be used to test the model
	nonCpGExample = chr22[y:y+len(CpGExample)+1] #will be used to test the model
	while 'N' in nonCpGExample: #so that the model will never get tested with a sequence containing N
		y = random.randint(0,len(chr22)) 
		nonCpGExample = chr22[y:y+len(CpGExample)+1]
	
	
	relativeFreqInside.to_csv('freqInside.txt', sep='\t') 
	relativeFreqOutside.to_csv('freqOutside.txt', sep='\t')
	with open('CpGExample.txt','a') as f:	
		f.write(CpGExample) #store in a txt file the relative frequences of the nucleotides in CpG islands
	with open('nonCpGExample.txt','a') as f: 
		f.write(nonCpGExample) #store in a txt file the relative frequences of the nucleotides in CpG islands
	return     

    
getFrequences('model-based-cpg-islands-hg19.txt','chr22.fa')