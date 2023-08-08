###Smith and Waterman algorithm for sequences alignment
#the user shoud type on the terminal the name of the file, the two strings, the scores for match, mismatch and gap
#the algorithm returns just one of the possible optimal alignments (if more than one is present)

import pandas as pd
import sys

'''mismatch = -1
match = 1
gap = -2'''

def align(matrix : pd.core.frame.DataFrame, x : int, y : int, match, mismatch):
	if (matrix.columns[y] == matrix.index[x]):
		score = match 
	else:
		score = mismatch
	return score
	
def SmithWaterman(string1 : str, string2 : str, match : int, mismatch : int, gap : int):
	matrix = pd.DataFrame(columns = [''] + [y for y in string2], index = [' ']+[x for x in string1])
	
	#step 1. Matrix initialization
	for i in range(0, len(string1)+1):
		matrix.iat[i,0] = 0 
	for j in range(0, len(string2)+1):
		matrix.iat[0,j] = 0
	
	#step 2. Matrix filling
	for i in range(1,len(string1)+1):
		for j in range(1, len(string2)+1):
			v1 = matrix.iat[i-1,j-1] + align(matrix, i, j, match, mismatch)
			v2 = matrix.iat[i-1, j] + gap
			v3 = matrix.iat[i, j-1] + gap
			v4 = 0
			
			matrix.iat[i,j] = max(v1, v2, v3, v4)
	#step 3. Traceback
	maxScore = int((matrix.max().values).max())
	for y in range(len(string2)-1,0,-1):
		for x in range(len(string1)-1,0,-1):
			if matrix.iat[x,y] == maxScore:
				j = y
				i = x
				break
                
	print(matrix)
	#j = matrix.columns[matrix[i][:]==maxScore]
	#print(j,i)
	globalScore = matrix.iat[i,j]
	align1 = []
	align2 = []
	
	while i>0 and j>0:
		if ((matrix.iat[i-1,j-1] == 0 and matrix.iat[i-1,j] == 0) and matrix.iat[i,j-1] == 0):
			align1.append(matrix.index[i])
			align2.append(matrix.columns[j])
			break        
		elif matrix.iat[i,j] == (matrix.iat[i-1, j-1]+align(matrix,i,j, match, mismatch)):
			align1.append(matrix.index[i])
			align2.append(matrix.columns[j])
			i-=1
			j-=1
		elif matrix.iat[i,j] == (matrix.iat[i-1, j] + gap):
			align1.append(matrix.index[i])
			align2.append('-')
			i-=1
		elif matrix.iat[i,j] == (matrix.iat[i, j-1] + gap):
			align1.append('-')
			align2.append(matrix.columns[j])
			j-=1
		else:       #this else stands for all the cases in which the actual value of the cell was negative, therefore none of the neighbors could have produced it
			break  
            	
	return f"\n{' '.join(align1[::-1])} \n{' '.join(align2[::-1])} \nglobal score = {globalScore}" #here the lists are transformed into strings

if __name__ =='__main__':
	if len(sys.argv) !=6:
		print('Usage: %s string1 string2 match mismatch gap' % sys.argv[0])
		sys.exit()
	str1 = sys.argv[1]
	str2 = sys.argv[2]
	m = int(sys.argv[3])
	mism = int(sys.argv[4])
	g = int(sys.argv[5])
    
	print(SmithWaterman(str1,str2, m, mism, g))
    
