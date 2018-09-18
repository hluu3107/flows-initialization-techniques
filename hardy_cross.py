import numpy as np
from helper import *
import time

def hardy_cross_adjustment(cycle,flows,resistants):
	numerator = 0.0
	denominator = 0.0	
	for (u,v) in cycle:
		r = resistants.get((u,v))
		flow = 0.0
		direction = 0

		if (u,v) in flows.keys():
			flow = flows.get((u,v))
			direction = 1
		else:
			flow = flows.get((v,u))

		if direction==1:
			numerator += r*flow**2
		else:
			numerator -= r*flow**2
		denominator += 2*r*flow
	flow_adjustment = numerator/denominator
	return flow_adjustment

def hardy_cross_update(cycle,flows,flow_adjustment):
	for (u,v) in cycle:
		f = 0.0
		
		direction = 0
		if (u,v) in flows.keys():
			f = flows.get((u,v))
			direction = 1
		else:
			f = flows.get((v,u))
		
		if direction == 1:
			f = f-flow_adjustment
		else:
			f = f+flow_adjustment

		if f < 0:
			f = np.absolute(f)
			direction = 1 if direction == 0 else 0
			if (u,v) in flows.keys():
				flows.pop((u,v),None)
				flows[(v,u)] = f
			else:
				flows.pop((v,u),None)
				flows[(u,v)] = f
		
		if direction == 1:
			flows[(u,v)] = f
		else:
			flows[(v,u)] = f

#generate initial guess by finding a path using dfs
def inital_guess1(tree,visited,edge_list,source,sink,initial_flow):
	path = find_path_tree(tree,source,sink)
	flows = {}
	for i in range(0,len(path)-1):
		flows[(path[i],path[i+1])] = initial_flow
	
	for (u,v) in edge_list:
		if (u,v) not in flows.keys() and (v,u) not in flows.keys():
			flows[(u,v)] = 0.0
	#print(f'flows is {flows}')
	return flows

def solve_potentials(adj_matrix,tree,flows,resistants):
	A = np.zeros([len(adj_matrix),len(adj_matrix)],dtype=int)
	B = np.zeros([len(adj_matrix)])
	B[-1] = 0
	A[-1,-1] = 1
	count = 0	
	for (u,v) in tree:
		resistant = resistants.get((u,v))
		if((u,v) in flows):
			flow = flows.get((u,v))
			A[count,u] = 1
			A[count,v] = -1
			B[count] = resistant*flow*flow
		else:
			flow = flows.get((v,u))
			A[count,u] = -1
			A[count,v] = 1
			B[count] = resistant*flow*flow
		count+=1
	potentials = np.linalg.solve(A, B)
	return potentials

def printHardy(flow_values,node_potentials):
	print('\nNode potentials:')
	for i in range(0,len(node_potentials)):
		print(f"{i}: {node_potentials[i]:.2f}")
	print('\nFlow:')
	for (i, j) in flow_values.keys():
		flow = flow_values[i, j]
		prediction = max(0, node_potentials[i] - node_potentials[j])**0.5
		print(f"{i} -> {j}: {flow:.2f}; (= {prediction:.2f})" if prediction > 0 else ")")

def hardy_cross(adj_matrix,initial_flow_value,resistants,timeout):
	start_time = time.time()
	(tree,visited) = dfs_tree(adj_matrix,0)
	edge_list = get_edge_list_undirected(adj_matrix)
	cycles = get_cycle(tree,visited,edge_list)
	
	iteration = 0;
	max_iteration = 1e6;
	tol = 1E-4

	initial_flow = initial_flow_value
	tol_error = initial_flow*0.01
	source = 0
	sink = len(adj_matrix)-1

	flows = inital_guess1(tree,visited,edge_list,source,sink,initial_flow)
	#print(flows)
	#resistants = {(0,1):1.0,(1,0):1.0,(0,2):1.0,(2,0):1.0,(1,2):100000.0,(2,1):100000.0,(1,3):1.0,(3,1):1.0,(2,3):1.0,(3,2):1.0}
	#flows = {(0,2):5.0,(2,3):5.0,(0,1):5.0,(1,3):5.0,(1,2):0.0}
	#flows = {(0,1):1000,(1,2):500,(1,4):500,(2,3):1000,(4,3):1000,(3,5):2000,(5,2):500,(5,4):500,(5,6):1000}
			
	flags = np.zeros(len(cycles))
	iteration = 1;
	
	stop_before = False
	#Using time as stopping criteria
	while(time.time()-start_time<=timeout):
		i=0
		for cycle in cycles:
			flow_adjustment = hardy_cross_adjustment(cycle,flows,resistants)
			if flow_adjustment == 0:
				flags[i] = 1
			else:
				flags[i] = 0
			i = i+1
			hardy_cross_update(cycle,flows,flow_adjustment)
		if any(flag == 0 for flag in flags):
			continue
		else:				
			stop_before = True
		iteration = iteration + 1
	#potentials = solve_potentials(adj_matrix,tree,flows,resistants)

	# for iteration in range(1,max_iteration):
	# 	#print(f'iteration {iteration}: flow is {flows}')
	# 	i = 0
	# 	for cycle in cycles:
	# 		flow_adjustment = hardy_cross_adjustment(cycle,flows,resistants)
	# 		#print(f'adjustment for cycle {i} is {flow_adjustment}')
			
	# 		if flow_adjustment == 0 or np.absolute(flow_adjustment) <= tol:
	# 			flags[i] = 1
	# 		else:
	# 			flags[i] = 0

	# 		i+=1
	# 		hardy_cross_update(cycle,flows,flow_adjustment)			
		
	# 	if any(flag == 0 for flag in flags):
	# 		continue
	# 	else:
	# 		potentials = solve_potentials(adj_matrix,tree,flows,resistants)
	# 		max_error = calculate_error(flows,potentials,resistants)
	# 		if max_error < tol_error:
	# 			break
	# 		else:
	# 			continue

	# potentials = solve_potentials(adj_matrix,tree,flows,resistants)
	# max_error = calculate_error(flows,potentials,resistants)

	#print(f'adj matrix is {adj_matrix}')
	#print(f'cycles are {cycles}')
	#print(f'stop after {iteration} iterations and max error is {max_error}')
	#print(f'final flow is {flows}')
	
	return flows
