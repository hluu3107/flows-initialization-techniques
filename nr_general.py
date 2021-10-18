import numpy as np
import networkx as nx
import random
from itertools import *
import matplotlib.pyplot as plt
from gurobipy import *

def randomGraph(n,p):
	edges = combinations(range(n), 2)
	G = nx.Graph()
	G.add_nodes_from(range(n))
	if p <= 0: return G
	if p >= 1: return nx.complete_graph(n, create_using=G)
	for _, node_edges in groupby(edges, key=lambda x: x[0]):
		node_edges = list(node_edges)
		random_edge = random.choice(node_edges)
		G.add_edge(*random_edge)
		for e in node_edges:
			if random.random() < p:
				G.add_edge(*e)
	return G

def iter_newton(A,x,psi,imax=100,tol=1e-6):
	iter = 1
	for i in range(imax):
		J = jacobian(A,x,psi)
		Y = function(A,x,psi)

		if np.linalg.matrix_rank(J) < len(x):
			#print("not full rank")
			inv = np.linalg.inv(J)
			x = np.add(np.dot(inv,-Y),x)
		else:
			dx = np.linalg.solve(J,-Y)
			x = np.add(x,dx)
		iter+=1
		if np.linalg.norm(Y) < tol:
			break
	print("numberof iterations ", iter)
	return x

def ssquare(a):
	return np.multiply(np.absolute(a),a)

def flows(A,x,psi):
	q = np.add(psi,np.dot(A,x))
	return q

def function(A,x,psi):
	return np.dot(A.transpose(),ssquare(flows(A,x,psi)))

def jacobian(A,x,psi):
	abs_flows = np.absolute(flows(A,x,psi))
	diag_flows = np.diag(abs_flows)
	J = 2*np.dot(np.dot(A.transpose(),diag_flows),A)
	return J

def lngp(nodes,edges,outflow):
	bound = sum(outflow[key] for key in outflow if outflow[key]>0)
	m = Model('netflow')
	m.Params.LogToConsole = 0
	#vars
	flows = m.addVars(edges,name="flows",lb = -bound, ub = bound)
	maxflow = m.addVar(name="maxflow")
	#aux vars for abs
	absflows = m.addVars(edges,name="absflows")

	#flow conservations constraints
	flows_conservation = {v: m.addConstr(flows.sum(v, '*') == outflow.get(v, 0) + flows.sum('*', v)) for v in nodes}
	#max abs constrains
	for (u,v) in edges:
		#m.addConstr(maxflow >= flows[u,v]) 
		m.addConstr(absflows[u,v]== abs_(flows[u,v]))
		m.addConstr(maxflow >= absflows[u,v]) 
	
	m.setObjective(maxflow, GRB.MINIMIZE)
	m.update()
	m.optimize()
	return m.getAttr('x',flows)

def lnpg2(nodes,edges,outflow):
	
	bound = sum(outflow[key] for key in outflow if outflow[key]>0)
	epsilon = 1e-4
	print(epsilon)
	m = Model('netflow2')
	#m.Params.LogToConsole = 0
	#vars
	flows = m.addVars(edges,name="flows",lb = -bound, ub = bound)
	#aux vars for abs
	absflows = m.addVars(edges,name="absflows")
	#flow conservations constraints
	flows_conservation2 = {v: m.addConstr(flows.sum(v, '*') == outflow.get(v, 0) + flows.sum('*', v)) for v in nodes}
	#max abs constrains
	for (u,v) in edges:
		#m.addConstr(maxflow >= flows[u,v]) 
		m.addConstr(absflows[u,v]== abs_(flows[u,v]))
		m.addConstr(absflows[u,v] >= epsilon) 

	m.setObjective(quicksum(absflows), GRB.MINIMIZE)
	m.update()
	m.optimize()
	return m.getAttr('x',flows)

def genGraph(minn,maxn, q):
	nnodes = random.randint(minn,maxn)
	G = randomGraph(nnodes,q)
	nedges = G.number_of_edges()
	cycles = nx.cycle_basis(G)
	adj_matrix=nx.convert.to_dict_of_lists(G)
	basis = []
	lookup = {}
	counter = 0
	ncycles = len(cycles)
	A = np.zeros((nedges,ncycles))
	for u,neighbors in adj_matrix.items():
		for v in neighbors:
			if u<v:
				lookup[(u,v)] = counter
				counter+=1
	counter = 0
	for cycle in cycles:
		cur = []
		for i in range(0,len(cycle)-1):
			u,v = cycle[i],cycle[i+1]
			cur.append((u,v))
			if u<v:
				A[lookup[(u,v)]][counter]=1
			else:
				A[lookup[(v,u)]][counter]=-1
		first,last = cycle[-1],cycle[0]
		cur.append((first,last))
		if first < last:
			A[lookup[(first,last)]][counter]=1
		else:
			A[lookup[(last,first)]][counter]=-1
		counter+=1
		basis.append(cur)
	#print(basis)
	
	return G, adj_matrix, A, lookup

def init2(G,lookup, source,target, initflow=1):
	epsilon = initflow*1e-6
	all_paths = list(nx.all_simple_paths(G,source,target,cutoff=G.number_of_nodes()//2))
	shortest = list(nx.bidirectional_shortest_path(G,source,target))
	#other_paths = all_paths.difference(sps)
	#print(shortest)
	short_paths = []
	other_paths = []
	for p in all_paths:
		path = []
		for i in range(0,len(p)-1):
			u,v = p[i],p[i+1]
			path.append((u,v))
		if len(path)==len(shortest)-1:
			short_paths.append(path)
		else:
			other_paths.append(path)

	#print(other_paths)
	#print(short_paths)

	other_value = 0 if not other_paths else epsilon/len(other_paths)
	short_value = (1-epsilon)/len(short_paths)
	psi = [0]*len(lookup)
	
	for sp in short_paths:
		for (u,v) in sp:
			if u<v:
				psi[lookup[(u,v)]] += short_value
			else:
				psi[lookup[(v,u)]] -= short_value

	for path in other_paths:
		for(u,v) in path:
			if u<v:
				psi[lookup[(u,v)]] += other_value
			else:
				psi[lookup[(v,u)]] -= other_value
	return psi


def main():
	G, adj_matrix,A, lookup = genGraph(15,15, 0.3)
	k = len(lookup) - len(adj_matrix)+1
	source, sink = 0, len(adj_matrix)-1
	flow_value = 1

	nodes = [n for n in adj_matrix.keys()]
	edges = []
	outflow = {source: flow_value, sink: -flow_value}
	for u,neighbors in adj_matrix.items():
		for v in neighbors:
			if u<v:
				edges.append((u,v))

	#send some epsilon flow on every edge
	small_flow = flow_value/len(nodes)
	outflow_small = {source: small_flow, sink: -small_flow}
	flows1 = lnpg2(nodes,edges,outflow_small)
	#print(flows1)
	psi = [0]*len(lookup)

	for (u,v),flow in flows1.items():
		if u>v: u,v = v,u
		psi[lookup[(u,v)]] = flow
	
	#send the rest of the flow on the shortest path
	shortest = list(nx.bidirectional_shortest_path(G,source,sink))
	leftover_flow = flow_value-small_flow
	for i in range(0,len(shortest)-1):
		u,v = shortest[i],shortest[i+1]
		if u<v:
			psi[lookup[(u,v)]] += leftover_flow
		else:
			psi[lookup[(v,u)]] -= leftover_flow

	x = [0]*k
	r1 = iter_newton(A,x,psi)
	print(flows(A,r1,psi))

	#minimize the max flow on every edge
	flows2 = lngp(nodes,edges,outflow)
	#print(flows2)
	psi2 = [0]*len(lookup)
	for (u,v),flow in flows2.items():
		if u>v: u,v = v,u
		psi2[lookup[(u,v)]] = flow
	x = [0]*k
	r2 = iter_newton(A,x,psi2)
	print(flows(A,r2,psi2))
	plt.figure(figsize=(8,5))
	nx.draw(G, node_color='lightblue', with_labels=True, node_size=500)
	plt.show()
def main1():
	#2 cycles
	epsilon = 1e-4
	#adj_matrix = {0:[1,2,3], 1:[0,2], 2:[0,1,3], 3:[0,2]}
	# A = np.transpose(np.array([[1,1,-1,0,0],[0,0,-1,1,1]]))
	# psi = np.array([epsilon,epsilon,1-epsilon,0,0])
	# #psi = np.array([1/3,1/3,1/3,1/3,1/3])
	# x = np.array([0,0])

	#3cycles
	#A = np.transpose(np.array([[-1,-1,1,0,0,0],[0,0,-1,-1,1,0],[0,1,0,1,0,-1]]))
	#adj_matrix = {1:[2,3,4], 2:[1,3,4], 3:[1,2,4], 4:[1,3,2]}
	#psi = np.array([-1/3,1/3,0,1/3,1/3,1/3])
	#psi = np.array([0,0,epsilon,epsilon,1-epsilon,0])
	#x = np.array([0,0,0])
	
	#plt.figure(figsize=(8,5))
	#nx.draw(G, node_color='lightblue', with_labels=True, node_size=500)
	#plt.show()
	

	G, adj_matrix,A, lookup = genGraph(5,10, 0.5)
	k = len(lookup) - len(adj_matrix)+1
	
	source, sink = 0, len(adj_matrix)-1
	flow_value = 1
	outflow = {source: flow_value, sink: -flow_value}
	nodes = [n for n in adj_matrix.keys()]
	edges = []
	
	for u,neighbors in adj_matrix.items():
		for v in neighbors:
			if u<v:
				edges.append((u,v))
	flows = lngp(nodes,edges,outflow)
	
	psi = [0]*len(lookup)
	for (u,v),flow in flows.items():
		if u>v: u,v = v,u
		psi[lookup[(u,v)]] = flow
	x = [0]*k
	iter_newton(A,x,psi)


if __name__ == '__main__':
    main()
