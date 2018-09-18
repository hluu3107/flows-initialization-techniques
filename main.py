from hardy_cross import *
from optimization import *
from helper import *
import time

def main():
	inputFile = "test"
	#adj_matrix = process_input(inputFile)
	adj_matrix = generate_random_graph(50,500)
	edge_list = get_edge_list_undirected(adj_matrix)
	resistants = random_resistant(edge_list,1,1)
	source = 0
	sink = len(adj_matrix)-1
	(tree,visited) = dfs_tree(adj_matrix,source)
	
	#print('done generate graph')
	initial_flow_value = 10
	
	
	start_time2 = time.time()
	(flow_values,potentials2) = solve_flow(adj_matrix, initial_flow_value)
	end_time2 = time.time()
	

	print(f'Hardy cross')
	timeout = end_time2-start_time2 #second
	start_time1 = time.time()
	flows1 = hardy_cross(adj_matrix,initial_flow_value,resistants,timeout)
	end_time1 = time.time()
	potentials1 = solve_potentials(adj_matrix,tree,flows1,resistants)
	error1 = calculate_error(flows1,potentials1,resistants)
	print(f'Error 1: {error1} + time spend: {end_time1-start_time1:.2f}')
	#printHardy(flows1,potentials1)
	
	print(f'Optimization')
	flows2 = {}
	for (u,v) in flows1.keys():
		flows2[(u,v)] = flow_values.get((u,v))
	#printResult(flows2,potentials2)
	error2 = calculate_error(flows2,potentials2,resistants)
	print(f'Error 2: {error2} + time spend: {end_time2-start_time2:.2f}')
if __name__ == '__main__':
    main()