import numpy as np
import random
from collections import defaultdict

def process_input(inputFile):
	#process input file and return adjacency matrix
	f = open(inputFile, "r")
	adj_matrix = {}
	for line in f:
		st = line.split(":")
		node1 = int(st[0])		
		st2 = st[1].strip().split(",")
		neighbors =list(map(int,st2))
		adj_matrix[node1]= [n for n in neighbors]	
	return adj_matrix

def generate_random_graph(nvertices,nedges):
	adj_matrix = defaultdict(list)
	vertices = set()
	edges = set()
	visited = []
	for i in range (1,nvertices-1):
		vertices.add(i)
	#print(vertices)
	for i in range (1,nvertices-1):
		for j in range(i+1,nvertices-1):
			edges.add((i,j))
	#print(edges)
	currentVertex = 1
	vertices.remove(currentVertex)
	count = 0;
	while vertices:
		vertex = random.sample(vertices,1).pop()
		vertices.remove(vertex)		
		adj_matrix[vertex].append(currentVertex)
		adj_matrix[currentVertex].append(vertex)
		if currentVertex < vertex: edges.remove((currentVertex,vertex))
		else: edges.remove((vertex,currentVertex)) 
		currentVertex = vertex	
		count+=1
	#add 0-1 and last to next to last as source and sink
	adj_matrix[0].append(1)
	adj_matrix[1].append(0)
	adj_matrix[nvertices-1].append(nvertices-2)
	adj_matrix[nvertices-2].append(nvertices-1)
	count+=2
	edge_list = random.sample(edges,nedges-count)
	for (u,v) in edge_list:
		adj_matrix[u].append(v)
		adj_matrix[v].append(u)
	return adj_matrix


def random_resistant(edge_list,min,max):
	resistants = {}
	for (u,v) in edge_list:
		res = random.randint(min, max)
		resistants[(u,v)] = res
		resistants[(v,u)] = res
	return resistants

def dfs_tree(adj_matrix,start_vertex):
	tree = [start_vertex]
	tree_edges = []
	dfs_visit(adj_matrix,start_vertex,tree,tree_edges)
	return (tree_edges,tree)

def dfs_visit(adj_matrix,v,tree,tree_edges):
	for w in adj_matrix[v]:
		if w not in tree:
			tree.append(w)
			tree_edges.append((v,w))
			dfs_visit(adj_matrix,w,tree,tree_edges)
def get_edge_list_undirected(adj_matrix):
	edge_list = []
	for v,neighbors in adj_matrix.items():
		for neighbor in neighbors:
			if v<neighbor:
				edge_list.append((v,neighbor))
	return edge_list

def find_path_tree(tree,start,end,path=[]):
	path = path + [start]
	if start == end:		
		return path
	for (u,v) in tree:
		if u == start and v not in path:
			newpath = find_path_tree(tree,v,end,path)
			if newpath:
				return newpath
	return None

def get_cycle(tree,visited,edge_list):	
	external_edge = []
	cycles = []
	for (x,y) in edge_list:
		if (x,y) not in tree and (y,x) not in tree:
			external_edge.append((x,y))
	
	for (x,y) in external_edge:
		idx_x = visited.index(x)
		idx_y = visited.index(y)
		if idx_x > idx_y:
			x,y = y,x

		path = find_path_tree(tree,x,y)
		cycle = []		
		for i in range(0,len(path)-1):
			(u,v) = (path[i],path[i+1])
			cycle.append((u,v))
		cycle.append((y,x))
		cycles.append(cycle)
	return cycles

def calculate_error(flows,pressure,resistants):
	max_error = 0.0
	for (u,v),flow in flows.items():
		flows_prediction = (np.absolute(pressure[u]-pressure[v]))**0.5/resistants.get((u,v))
		diff = np.absolute(flows_prediction-flow)
		if(diff > max_error):
			max_error = diff
	return max_error