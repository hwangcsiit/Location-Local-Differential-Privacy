from scipy.optimize import brent, fmin, minimize
import numpy as np
from collections import defaultdict

def Optgroup(origin):
    data = []
    for item in origin:
        data.append(item[0])
    b0 = (0, 46)  
    b1 = (0, 46)  
    b2 = (0, 46) 
    bnds = (b0, b1, b2)

    con1 = {'type': 'ineq', 'fun': constraint1}
    con2 = {'type': 'ineq', 'fun': constraint2}
    cons = ([con1, con2])

    x0 = np.array([46,30,2])  # set the initial value; very important to set
    res = minimize(objFun, x0, method='SLSQP', bounds=bnds, constraints=cons)
    final_group = defaultdict(list)

    for loc1 in data:
        G_1=[]
        G_2=[]
        G_3=[]
        for loc2 in data:
            group = comparestring(loc1, loc2, res.x)
            if group == 1:
                G_1.append(loc2)
            if group == 2:
                G_2.append(loc2)
            if group == 3:
                G_3.append(loc2)
        final_group[loc1].append(G_1)
        final_group[loc1].append(G_2)
        final_group[loc1].append(G_3)

    return final_group


def objFun(x):
	dataor = ("domain.txt")
	data = []
	for item in dataor:
		data.append(item[0])
	data = list(set(data))

	func = 0

	for loc1 in data:
		GLCP = {}
		for loc2 in data:
			group = comparestring(loc1, loc2, x)
			if group in GLCP:
				GLCP[group]+=1
			else:
				GLCP[group]=1

		for i in range(3):
			if i+1 in GLCP:
				func = func + x[i]*GLCP[i+1]
			else:
				func = func 
	return -func 


def constraint1(x):
    return x[0]-x[1]+1

def constraint2(x):
    return x[1]-x[2]+1

def comparestring(loc1, loc2, beta):
	com = 0
	group = 0
	for i in range(len(loc1)):
		if loc1[i]==loc2[i]:
			com = com + 1
		else:
			break
	for i in range (3):
		if com >= beta[i]:
			group = i+1
			break
	return group
