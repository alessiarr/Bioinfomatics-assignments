#Assignment 3, Bioinformatics
##string matching 

from BurrowsWheelerTransform import *
import random

def getGenome(length=1000):
    genome = ''.join(random.choice('ACGT') for i in range(length))
    return genome+'$'
    
    
def findMatches(genome,p): #where p is the string to match
	columns = (BWT(genome, '$'))[:2]
	f = indexing(columns[1])
	l = indexing(columns[0])
	
	matches = 0
	for i in range(len(f)):
		j = 1
		if f[i][0] == p[-j]:
			ind = i
			j+=1
			count = True
			while count and j!=len(p)+1:
				if l[ind][0] == p[-j]:
					ind = f.index(l[ind])
					j+=1
				else:
					count = False
		if j == len(p)+1:
			matches += 1
	
	return f'The string P matches {matches} times'

genome = getGenome(20)	#the third output is the offset array, do not need it here
p1 = ''.join(random.choice('ACGT') for i in range(5))
p2 = ''.join(random.choice('ACGT') for i in range(5))
#matches of p1
index = random.sample(range(20),1)
genome = genome[:index[0]] + p1 + genome[index[0]:]
print(genome,p1)
print(findMatches(genome, p1))

#matches of p2
index = random.sample(range(len(genome)),2)
index.sort()
genome = genome[:index[0]] + p2 + genome[index[0]:index[1]] + p2 + genome[index[1]:]
print(genome,p2)
print(findMatches(genome, p2))
