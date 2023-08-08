#Bioinformatics, Assignment 4
##produce the inside and outside models

import pandas as pd
import math 
import random

def getModels(file1,file2):
	freqInside = pd.read_csv(file1, sep='\t', index_col=0)
	#print(freqInside)
	freqOutside = pd.read_csv(file2, sep='\t',index_col=0)
	totInside = sum(freqInside.iat[4,i] for i in range(0,4))
	totOutside = sum(freqOutside.iat[4,i] for i in range(0,4))
	for i in range(0,4): #through the columns of the table
		for j in range(0,5): #through the rows of the table 
			freqInside.iat[j,i]/totInside #this gives the absolute frequency, that is also the P(X and Y)
			freqOutside.iat[j,i]/totOutside
		
	transitionInside = pd.DataFrame(columns=['A','C','G','T'],index=['A','C','G','T'] )
	transitionOutside = pd.DataFrame(columns=['A','C','G','T'],index=['A','C','G','T'] )
	initialStateInside = pd.Series(data=[freqInside.iat[4,i] for i in range(0,4)], index=['A','C','G','T']) # P(A), P(C), P(G) and P(T)
	initialStateOutside = pd.Series(data=[freqOutside.iat[4,i] for i in range(0,4)], index=['A','C','G','T'])
	for i in range(0,4): #rows
		for j in range(0,4): #columns
			transitionInside.iat[i,j] = freqInside.iat[i,j]/initialStateInside[i] #produce the conditional probabilities P(X|Y) = P(X and Y)/P(Y) where Y is the label of the row (that is the x(i-1) nucleotide) and X the label of the columns (the x(i) nucleotide)
			transitionOutside.iat[i,j] = freqOutside.iat[i,j]/initialStateOutside[i]
	
	print('INSIDE MODEL:\n',transitionInside,'\n\nOUTSIDE MODEL:\n',transitionOutside)
	return (transitionInside, initialStateInside, transitionOutside, initialStateOutside)

def testModel(query, transitionInside,initialStateInside,transitionOutside,initialStateOutside):
    probInside = initialStateInside[query[0]]
    probOutside = initialStateOutside[query[0]]
    for i in range(1,len(query)):
        probInside = probInside*transitionInside.at[query[i-1],query[i]]
        probOutside = probOutside*transitionOutside.at[query[i-1],query[i]]
        
    score = math.log(probInside/probOutside)
    if score > 0:
        print('Score: ',score,', the sequence is predicted to be a CpG island')
    else:
        print('Score: ', score, ', the sequence is predected not to be a CpG Island')
       

p = getModels('freqInside.txt', 'freqOutside.txt')
f = open('CpGExample.txt','r')
CpGExample = f.read()
f.close
f = open('nonCpGExample.txt','r')
nonCpGExample = f.read()
f.close
testModel(CpGExample.upper(),p[0],p[1],p[2],p[3]) #the model is tested with a real CpG Island sequence
testModel(nonCpGExample.upper(),p[0],p[1],p[2],p[3]) #the model is tested with a sequence randomly taken from chromosome 22 (likely is not a CpG Island)
testModel(''.join(random.choice('ACGT') for i in range(300)),p[0],p[1],p[2],p[3]) #random query sequence