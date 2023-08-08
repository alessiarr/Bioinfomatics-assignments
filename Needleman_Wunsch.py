###Needleman and Wunsch altorithm for sequences alignment 
#the user shoud type on the terminal the name of the file, the two strings, the scores for match, mismatch and gap
#the algorithm returns just one of the possible optimal alignments (if more than one is present)

import pandas as pd
import sys
#from IPython.display import display


'''mismatch = -1
match = 1
gap = -2'''

def align(matrix : pd.core.frame.DataFrame, x : int, y : int, match, mismatch):
	if (matrix.columns[y] == matrix.index[x]):
		score = match
	else:
		score = mismatch
	return score
	
def NeedlemanWunsch(string1 : str, string2 : str, match : int, mismatch : int, gap : int):
	matrix = pd.DataFrame(columns = [' ']+[d for d in string2], index = [' ']+[c for c in string1])
	#return matrix

	#step 1. Matrix initialization
	scoreRow = 0
	scoreColumn = 0
	for i in range(0,len(string1)+1):
		matrix.iat[i,0] = scoreRow
		scoreRow+= gap
	for j in range(0,len(string2)+1):
		matrix.iat[0,j] = scoreColumn
		scoreColumn += gap
	#return matrix
	##at this point the matrix should have the first row and column filled

	#step 2. Matrix filling
	for i in range(1, len(string1)+1):
		for j in range(1, len(string2)+1):
			v1 = matrix.iat[i-1,j-1] + (align(matrix, i, j, match, mismatch))
			v2 = matrix.iat[i-1,j] + gap
			v3 = matrix.iat[i,j-1] + gap
			matrix.iat[i,j] = max(v1, v2, v3)
	#print(matrix)
	##at this point all the cells should have a score assigned to
	
	#step 3. Traceback
	i = len(string1)
	j = len(string2)
	globalScore = matrix.iat[i,j]
	align1 = []
	align2 = []

	while i>0 and j>0:
		matrix.style.applymap(lambda x: 'color: red' if matrix.column == j and matrix.index == i else 'background-color:white')
		if matrix.iat[i,j] == (matrix.iat[i-1, j-1]+align(matrix,i,j, match, mismatch)):
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
	
	if i>0:
		for i in range(i,0,-1):
			align1.append(matrix.index[i])
			align2.append('-')		
	if j>0:
		for j in range(j,0,-1):
			align1.append('-')
			align2.append(matrix.columns[j])
	print(matrix)
	return f"\n{' '.join(align1[::-1])} \n{' '.join(align2[::-1])} \nglobal score = {globalScore}"

if __name__ =='__main__':
	if len(sys.argv) !=6:
		print('Usage: %s string1 string2 match mismatch gap' % sys.argv[0])
		sys.exit()
	str1 = sys.argv[1]
	str2 = sys.argv[2]
	m = int(sys.argv[3])
	mism = int(sys.argv[4])
	g = int(sys.argv[5])
   
	print(NeedlemanWunsch(str1,str2, m, mism, g))
    
