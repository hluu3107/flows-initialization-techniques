from gurobipy import *

nodes = "abcdef"

flow_value = 1

source, sink = 'a', 'f'

arcs = '''
ab ac
ba bd
ca cd ce
db dc df
ec ef
fd fe'''

nodes = list(nodes)
arcs = list(tuple(arc) for arc in arcs.split() if arc)

outflow = {source: flow_value, sink: -flow_value}

# Create optimization model
m = Model('netflow')
m.setParam('OutputFlag', False)

# Create variables
flow = m.addVars(arcs)

# Flow conservation constraints
constraints = {
    v: m.addConstr(flow.sum(v, '*') == outflow.get(v, 0) + flow.sum('*', v))
    for v in nodes
}

for arc in arcs:
    m.addConstr(flow.get(arc) >= 0)

# set objective to sum of flows cubed (piecewise-linear approximation)
# http://www.gurobi.com/documentation/7.5/refman/py_model_setpwlobj.html

m.setObjective(0, GRB.MINIMIZE)

low, high, n_pieces = 0, flow_value, 1000
xs = list(low + (high-low)*i/n_pieces for i in range(0, n_pieces+1))
ys = list(x**3 for x in xs)
for arc in arcs:
    m.setPWLObj(flow[arc], xs, ys)

# Compute optimal solution

m.update()
m.optimize()

m.write("model.Lp")
# Print solution
if m.status == GRB.Status.OPTIMAL:
    flow_values = m.getAttr('x', flow)
    node_potentials = m.getAttr('Pi', constraints)

    print('\nNode potentials:')
    for v in nodes:
        print(f"{v}: {node_potentials[v]:.2f}")


    print('\nFlow:')
    for i, j in arcs:
        flow = flow_values[i, j]

        prediction = max(0, node_potentials[i] - node_potentials[j])**0.5
        flows_final[(i,j)] = flow_values[i, j]
        print(f"{i} -> {j}: {flow:.2f}; (= {prediction:.2f}",
              f"* {flow/prediction:.2f})" if prediction > 0 else ")")
else:
    print("Unexpected exit status from Gurobi", m.status)


# expected output:

# Node potentials:
# a: 2.07
# b: 1.54
# c: 1.07
# d: 1.00
# e: 0.54
# f: 0.00

# Flow:
# a -> b: 0.42; (= 0.73 * 0.58)
# a -> c: 0.58; (= 1.00 * 0.58)
# b -> a: 0.00; (= 0.00 )
# b -> d: 0.42; (= 0.73 * 0.58)
# c -> a: 0.00; (= 0.00 )
# c -> d: 0.15; (= 0.27 * 0.58)
# c -> e: 0.42; (= 0.73 * 0.58)
# d -> b: 0.00; (= 0.00 )
# d -> c: 0.00; (= 0.00 )
# d -> f: 0.58; (= 1.00 * 0.58)
# e -> c: 0.00; (= 0.00 )
# e -> f: 0.42; (= 0.73 * 0.58)
# f -> d: 0.00; (= 0.00 )
# f -> e: 0.00; (= 0.00 )