import pandas as pd
import numpy as np
import sys

def get_edit_dist(a,b,b_lens): 
	#weighted by length of segment
	d = 0
	t = 0
	for j in range(len(a)):
		max_d = max([b_lens[a[j]], b_lens[b[j]]])
		if a[j]!=b[j]: 
			d+=max_d
		t = t + max_d

	return float(d)/t

def needleman_wunsch(list1, list2, match=1, mismatch=-1, gap=-2):
	"""
	Perform Needleman-Wunsch alignment on two lists of strings.

	Args:
		list1 (list of str): The first list of strings.
		list2 (list of str): The second list of strings.
		match (int): Score for a match.
		mismatch (int): Penalty for a mismatch.
		gap (int): Penalty for a gap.

	Returns:
		tuple: Aligned lists and alignment score.
	"""
	n = len(list1)
	m = len(list2)

	# Initialize scoring matrix
	score_matrix = np.zeros((n + 1, m + 1), dtype=int)
	
	# Initialize traceback matrix
	traceback = np.zeros((n + 1, m + 1), dtype=int)

	# Fill the scoring and traceback matrices
	for i in range(1, n + 1):
		score_matrix[i][0] = gap * i
	for j in range(1, m + 1):
		score_matrix[0][j] = gap * j

	for i in range(1, n + 1):
		for j in range(1, m + 1):
			match_score = score_matrix[i - 1][j - 1] + (match if list1[i - 1] == list2[j - 1] else mismatch)
			delete = score_matrix[i - 1][j] + gap
			insert = score_matrix[i][j - 1] + gap
			
			score_matrix[i][j] = max(match_score, delete, insert)

			# Traceback: 1=Diagonal, 2=Up, 3=Left
			if score_matrix[i][j] == match_score:
				traceback[i][j] = 1
			elif score_matrix[i][j] == delete:
				traceback[i][j] = 2
			else:
				traceback[i][j] = 3

	# Traceback to get alignment
	aligned_list1 = []
	aligned_list2 = []
	i, j = n, m

	while i > 0 or j > 0:
		if i > 0 and j > 0 and traceback[i][j] == 1:
			aligned_list1.append(list1[i - 1])
			aligned_list2.append(list2[j - 1])
			i -= 1
			j -= 1
		elif i > 0 and (j == 0 or traceback[i][j] == 2):
			aligned_list1.append(list1[i - 1])
			aligned_list2.append('...')
			i -= 1
		elif j > 0:  # Ensure j > 0 before appending from list2
			aligned_list1.append('...')
			aligned_list2.append(list2[j - 1])
			j -= 1

	return list(reversed(aligned_list1)), list(reversed(aligned_list2)), score_matrix[n][m]

def main():

	t = pd.read_csv(sys.argv[1], header=0, sep="\t")
	t_blens = pd.read_csv(sys.argv[2], header=0, sep="\t")
	b_lens = {}
	for i, r in t_blens.iterrows():
		b_lens[str(r['bID_d'])] = r["l"]
	b_lens['...'] = 0
	outrows = []
	for i,s1 in enumerate(t['hap_struc_d']):
		for j,s2 in enumerate(t['hap_struc_d']):
			a,b,score=needleman_wunsch(s1.split('-'), s2.split('-'))
			s=get_edit_dist(a,b,b_lens)
			outrows.append({"contig_a":t['contig'][i],
							"contig_b":t['contig'][j],
							"hap1":s1,
							"hap2":s2,
							"aln1":"-".join(a),
							"aln2":"-".join(b),
							"edit_d":s,
							"score":score})
	t = pd.DataFrame(outrows)
	t.to_csv(sys.argv[3], sep="\t", index=False)

if __name__ == '__main__':

	main()