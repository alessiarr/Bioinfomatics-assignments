#assignment 3, Bioinformatics
### BWT and reversal


import random


def getGenome(length=1000):
    genome = ''.join(random.choice('ACGT') for i in range(length))
    return genome+'$'
    
    
def BWT(s, specialCh):
	rotations = []
	for i in range (len(s)):
		s = s[len(s)-1] + s[:len(s)-1]
		rotations.append(s)	
	rotations.sort()
	offsetArray = []
	for e in rotations:
		i = 1
		while e[-i] != specialCh:
			i+=1
		offsetArray.append(i-1)
	#print(offsetArray)
	#print(rotations)
	
	return (''.join(map(lambda x: x[-1], rotations)), ''.join(map(lambda x: x[0], rotations)), offsetArray) #get both the last and first column of the bwm 

def indexing(column):
	#print(column, len(column))
	indexedColumn = list()
	for i in range(len(column)):
		j = column[i] + str(column[:i+1].count(column[i]))
		indexedColumn.append(j)
	#print(indexedColumn)
	return indexedColumn


def reversal(l,f):
	f = indexing(f)
	l = indexing(l)
	#print(l)
	t = [f[0][0]]    #t is the initial string, the special character is the last character of the initial string
	indexColumn = 0
	for i in range(1,len(l)):
		t.append(l[indexColumn][0]) #remove the index from each character
		#print(t)
		indexColumn = f.index(l[indexColumn]) #find the correspondent element in the first column   
	return ''.join(t[::-1])

    
'''genome = getGenome(20)
columns = (BWT(genome,'$'))[:2] #the third output is the offset array, we do not need it here
print(genome)
print(columns[0]) #this represents the transform
print(reversal(columns[0],columns[1]))'''