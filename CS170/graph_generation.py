import numpy as np 
from itertools import combinations, permutations
import random

"""		   2
		/	  \ 
	  1			3
	 /      /     \ 
	0      /	   7
	\     /       /
	 4           6
	  \         /
	  	    5 
The diamond's vertices are numbered like so"""


def diamond_fill(adj_matrix):
	"""Fills adjacency matrix with edge weights proposed in paper
	Vertex numbering follows clockwise movement. Goes in cycles of 0 to 7 (mod 8)"""
	counter = 0
	while (counter < k): #can make this into a for-loop as well
		for i in range(0,8):
			for j in range(0,8):
				x = counter*diamond + i #counter gives us which ``kth'' diamond we're at
				y = counter*diamond + j
				if (x == y):
					adj_matrix[x,y] = 0
					adj_matrix[y,x] = 0
				elif (i == west and j == sw_mid) or (i == sw_mid and j == west):
					adj_matrix[x,y] = 0
				elif (i == east and j == ne_mid) or (i == ne_mid and j == east):
					adj_matrix[x,y] = 0
				else:
					adj_matrix[x,y] = 1
		counter+=1
	#print(adj_matrix.item((0,0)))
	#print(adj_matrix)

def east_west(adj_matrix):
	"""Fills in the adj_matrix for the east-west connection between diamonds
	Temporarily hard-coded version"""
	for i in range(0, vertices):
		for j in range(0, vertices):
			if (abs(i-j)) == 1:
				if ((i%diamond == 7 or i%diamond == 0) and (j%diamond == 7 or j%diamond == 0)):
					#print((i,j))
					adj_matrix[i,j] = 1
	# adj_matrix[vertices-1, 0] = 1
	# adj_matrix[0, vertices-1] = 1
	#temporarily commenting the block out because we might not need to return to our start point 

def isolate(adj_matrix, edge_val=50):
	"""Isolates N_1 vertex as proposed in paper
	edge_val is currently set to 50 because we're reserving a bigger integer for other purposes
	mentioned in the paper but I haven't gotten to it."""
	N = set([north+8*x for x in range(k)])
	S = set([south+8*x for x in range(k)])
	N_sub = N - {north}
	NS = N_sub.union(S)
	edge_pairs = permutations(NS,2)
	for item in edge_pairs:
		#print(item)
		x = item[0]
		y = item[1]
		adj_matrix[x,y] = 0
	for vertex in NS:
		adj_matrix[north, vertex] = edge_val
		adj_matrix[vertex, north] = edge_val

def random_color_assignment(vertices):
	"""Randomly assigns the string of red and blues
	Might not be necessary as the conditions right now might generate no valid paths"""
	not_random, random_str = "", ""
	for i in range(vertices):
		if i%2==1:
			not_random+="R"
		else:
			not_random+="B"
	counter = 0
	while counter < vertices:
		r = random.randrange(1)
		if r%2 == 0:
			random_str+="R"
		else: 
			random_str+="B"
		counter+=1
	return not_random

if __name__ == '__main__':	
	vertices = 16	
	diamond = 8 #8 vertices in a diamond
	k = vertices//diamond #how many diamond circuits we'll have

	north = 2
	south = 5
	west = 0
	east = 7
	sw_mid = 4
	ne_mid = 3

	max_edge_val = 100

	adj_matrix = np.empty(shape=[vertices, vertices])
	adj_matrix.fill(0)

	diamond_fill(adj_matrix)
	isolate(adj_matrix)
	east_west(adj_matrix)
	colors = random_color_assignment(vertices)
	adj_matrix = adj_matrix.astype(int)
	with open("1.in", "w") as f:
		f.write(str(vertices))
		f.write('\n')
		for row in adj_matrix:
			for elem in row:
				f.write(str(elem)+ " ")
			f.write('\n')
		f.write(colors)
	f.close()

"""TODO: Check what the paper means to 'set all the other edges to 2M' (#3) 
Check about the isolating function. I am not getting (n-2 choose k)-k+1 edges so I'm not sure what the subgraph means
"""