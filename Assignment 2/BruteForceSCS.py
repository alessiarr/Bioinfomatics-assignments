##assignment 2 Bioinformatics, part2
import itertools
import random

def getGenome(length=1000):
    genome = ''.join(random.choice('ACGT') for i in range(length))
    return genome

def getSubstrings(seq,length=100):
    L = []
    for i in range(len(seq)-length+1):
        L.append(seq[i:i+length])
    return set(L)
    
def BruteForceSCS(strings):
	characters = 0
	for string in strings:
		for a in range(len(string)):
			characters+=1
	bestSuperstring = 'A'*characters #create a string to contain the worst-case scenario, i.e. the final superstring has lenght equal at the sum of the characters of the strings
	
	for permutation in itertools.permutations(strings):	#all the possible combinations are produced
		#print(permutation)
		shortestCS = permutation[0]
		for i in range(1,len(permutation)):
			score = 0
			for j in range(1,len(permutation[i])): #range up to the n-1 character of the string (= all the prefexes)
				if shortestCS[len(shortestCS)-j:] == permutation[i][:j]: #suffix-prefix overlap
					score = j
			
			shortestCS = shortestCS + permutation[i][score:] #takes into account both the merging and the concatenation (= i is 0)
            
		if len(bestSuperstring) > len(shortestCS):
			bestSuperstring = shortestCS
	#print(bestSuperstring)
	return bestSuperstring


'''string = getGenome(length=14)
print(string)
stringset = getSubstrings(string,length = 6)
print(stringset)
print(BruteForceSCS(stringset))'''
#print(BruteForceSCS({'ATCGGA','TACCCA', 'AGCTAC', 'CGGATT', 'TTGCTA'}))
