import numpy as np
from helper import *
from hardy_cross import inital_guess1
import time

def ssquare(a):
	return np.multiply(np.absolute(a),a)

def get_flows(A,x,psi):
	flows = np.add(psi,np.dot(A,x))
	return flows

def f(A,x,psi):
	return np.dot(A.transpose(),ssquare(get_flows(A,x,psi)))

def fp(A,x,psi):
	abs_flows = np.absolute(get_flows(A,x,psi))
	diag_flows = np.diag(abs_flows)
	J = 2*np.dot(np.dot(A.transpose(),diag_flows),A)
	#print(J)
	return J

def error(a,b):
	return np.sum(np.absolute(np.subtract(a,b)))

def nr(A,x,psi,x_init=1,tol=1e-3,iter=20):
	i = 0
	while True:
		F = f(A,x,psi)
		Fp = fp(A,x,psi)
		xold = x
		#check if Jacobian is singular
		if np.isfinite(np.linalg.cond(Fp)):
			delta = np.linalg.solve(Fp,-F)
			x = np.add(x,delta)
		else:
			print("singular Jacobian")
			x = 1/2*x
		e = error(F,f(A,x,psi))
		e1 = np.sum(np.absolute(x))
		if i==1:
			print('i: ', i, ' error: ', round(e,4)/(x_init))
			#print('norm x1: ', e1)
		else:
			print('i: ', i, ' error: ', round(e,4))
			#print('norm x0: ', e1, '\n')
			
		#print('flows: ', np.round(get_flows(A,x,psi),2),'\n')
		i+=1
		if e < tol or i > iter:
			break
	return x

def graph_nr(adj_matrix,x_init=1):
	source = 0
	sink = len(adj_matrix)-1
	(tree,visited) = dfs_tree(adj_matrix,0)
	edge_list = get_edge_list_undirected(adj_matrix)
	cycles = get_cycle(tree,visited,edge_list)

	lookup = {}
	counter = 0
	for i in range(len(cycles)):
		cycle = cycles[i]
		for (u,v) in cycle:
			if u> v:
				u,v = v,u
			if (u,v) not in lookup:
				lookup[(u,v)] = counter
				counter+=1

	nedges = len(lookup)
	k = len(cycles)
	A = np.zeros((nedges,k))
	#random flows
	#psi = np.random.randint(low=1,high=4,size=(nedges))
	psi = np.full(nedges,0)
	x = np.full(k,x_init)
	
	#x = np.random.randint(low=-100,high=100,size=(k))
	#random flows on one inflows
	# psi = np.zeros(nedges)
	# flows=inital_guess1(tree,visited,edge_list,source,sink,10)
	# for (u,v),f in flows.items():
	# 	if u>v: u,v = v,u
	# 	if (u,v) in lookup:
	# 		psi[lookup[(u,v)]] = f

	for i in range(len(cycles)):
		cycle = cycles[i]
		for (u,v) in cycle:
			for (u,v) in cycle:
				if u>v:
					edge = lookup[(v,u)]
					A[edge][i] = -1
				else:
					edge = lookup[(u,v)]
					A[edge][i] = 1
	sol = nr(A,x,psi,x_init)
	return sol

def generate_random_nr(n,m):
	adj_matrix = generate_random_graph(n,m)
	sol = graph_nr(adj_matrix,10)
	return sol
	
def main():
	A = np.transpose(np.array([[-1,-1,1,0,0,0],[0,0,-1,-1,1,0],[0,1,0,1,0,-1]]))
	psi = np.array([1,1,0,1,0,2])
	x = np.array([1.38,1,0.927])
	initf = f(A,x,psi)
	initj = fp(A,x,psi)
	inv = np.linalg.inv(initj)
	print(inv)
	alpha = np.max(np.abs(initf))
	m = float('inf')
	for i in range(len(initj)):
		s = initj[i][i] - sum(abs(initj[i][j]) for j in range(len(initj)) if j!=i)
		#print(s)
		m = min(m,s)
	beta = np.max((np.abs(inv)).sum(axis=1, dtype='float'))
	h = alpha * beta * beta * 32
	print(alpha)
	print(beta)
	print(h)
	ans = nr(A,x,psi)
	print(ans)
	print(get_flows(A,ans,psi))
if __name__ == '__main__':
    main()





