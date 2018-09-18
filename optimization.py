from gurobipy import *
import time

def solve_flow(adj_matrix, initial_flow_value):
	source, sink = 0, len(adj_matrix)-1
	flow_value = initial_flow_value
	outflow = {source: flow_value, sink: -flow_value}
	nodes = [n for n in adj_matrix.keys()]
	arcs = []
	for u,neighbors in adj_matrix.items():
		for neighbor in neighbors:
			arcs.append((u,neighbor))
	# Create optimization model
	m = Model('netflow')
	m.setParam('OutputFlag', False)
	
	# Create variables
	flow = m.addVars(arcs)
	# Constrains
	constraints = {
    	v: m.addConstr(flow.sum(v, '*') >= outflow.get(v, 0) + flow.sum('*', v))
    	for v in nodes
	}	
    
    # set objective to sum of flows cubed (piecewise-linear approximation)
	# http://www.gurobi.com/documentation/7.5/refman/py_model_setpwlobj.html
	m.setObjective(0, GRB.MINIMIZE)
	low, high, n_pieces = 0, flow_value, int(initial_flow_value*1e3)
	xs = list(low + (high-low)*i/n_pieces for i in range(0, n_pieces+1))
	ys = list(x**3/3 for x in xs)

	for arc in arcs:
		m.setPWLObj(flow[arc], xs, ys)
	#m.setParam("TimeLimit")
	m.update()
	m.optimize()
	if m.status == GRB.Status.OPTIMAL:
		flow_values = m.getAttr('x', flow)
		node_potentials = m.getAttr('Pi', constraints)	
		return flow_values,node_potentials			
	else:
		print("Terminate by violation: ", m.status)
	

def printResult(flow_values,node_potentials):
	print('\nNode potentials:')
	for v in node_potentials.keys():
		print(f"{v}: {node_potentials[v]:.2f}")
	print('\nFlow:')
	for (i, j) in flow_values.keys():
		flow = flow_values[i, j]
		prediction = max(0, node_potentials[i] - node_potentials[j])**0.5
		print(f"{i} -> {j}: {flow:.2f}; (= {prediction:.2f})" if prediction > 0 else ")")