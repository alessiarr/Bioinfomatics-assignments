##assignment 2 Bioinformatics, part1
import pandas as pd
import random
from BruteForceSCS import * #if we want to compare the results of the two algorithms

def findOverlap(string1, string2):
    #suffixes of string1 and prefixes of string2
    score = 0    
    for i in range(3,min(len(string1),len(string2))): #index suffix
        if string1[len(string1)-i:] == string2[:i]:
            score = i
    return score 

def getGenome(length=1000):
    genome = ''.join(random.choice('ACGT') for i in range(length))
    return genome

def getSubstrings(seq,length=100):
    L = []
    for i in range(len(seq)-length+1):
        L.append(seq[i:i+length])
    return set(L)

def greedySCS(strings : set):
	matrixSCS = pd.DataFrame(columns=[i for i in strings], index=[j for j in strings]) #matrix representing the overlap scores
    #step 1: pairwise suffix-prefix alignments
	for s1 in strings:
		for s2 in strings:
			if s1 == s2:
				matrixSCS.loc[s1,s2] = -1 #so that it will never be the maximum value
			else:
				matrixSCS.loc[s1,s2] = findOverlap(s1,s2)
	print(matrixSCS)
	maxScore = int((matrixSCS.max().values).max()) 
	#step 2: merge overlapping sequences
	while maxScore >0:
		for i in matrixSCS.index:
			for j in matrixSCS.columns:
				if matrixSCS.loc[i,j] == maxScore:
					x = i
					y = j
		mergedString = x + y[maxScore:]
		matrixSCS.drop([x,y],axis=0,inplace=True)
		matrixSCS.drop([x,y],axis=1,inplace=True)
		matrixSCS[mergedString] = -1 #add the new column of merged sequence
		matrixSCS.loc[mergedString] = -1 #add the new row of merged sequences
		#print(matrixSCS)
        #RECOMPUTE SCORES FOR THE NEW ROW AND COLUMN
		for i in matrixSCS.index.drop(mergedString):
			matrixSCS.loc[mergedString,i] = findOverlap(mergedString,i)
			matrixSCS.loc[i,mergedString] = findOverlap(i, mergedString)
		#print(matrixSCS)
		maxScore = int((matrixSCS.max().values).max()) 
	return 'The contigs are:', [i for i in matrixSCS.index]
        
#print(greedySCS(getSubstrings(getGenome(length=14), length=6)))
'''string = getGenome(length=14)
print(string)
stringset = getSubstrings(string,length = 6)
print(stringset)
print(BruteForceSCS(stringset)) 
print(greedySCS(stringset))'''
print(greedySCS({'ATCGGA','TACCCA', 'AGCTAC', 'CGGATT', 'TTGCTA'}))


