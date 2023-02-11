import matplotlib.pyplot as plt
import numpy as np
import math 
import random
from numpy import linalg as LA
import Partition


def readall(fname):
    with open(fname) as f:
        all_lines = f.readlines()
    all_lines = [x.strip().split() for x in all_lines]
    result=[]
    for line in all_lines:
        result.append([x.strip(',') for x in line])
    return result

def Hadarmard_init(k):
    H = [None] * k
    for row in range(k):
        H[row] = [True] * k
        
# Initialize Hadamard matrix of order n.
    i1 = 1
    while i1 < k:
        for i2 in range(i1):
            for i3 in range(i1):
                H[i2+i1][i3]    = H[i2][i3]
                H[i2][i3+i1]    = H[i2][i3]
                H[i2+i1][i3+i1] = not H[i2][i3]
        i1 += i1
    return H

def data_to_matrix(filename):
    file=open(filename)
    lines=file.readlines()
    rows=len(lines)
    datamat1=np.ones((rows,46))
    datamat = datamat1.astype(np.str)
    row=0
    for line in lines:
        line=line.strip().split(', ')
        datamat[row,:]=line[:]
        row+=1
    return datamat

def DE(domain, user,p,q):
    if(round(np.random.uniform(0, 1), 10)< p-q):
        return user
    else:
        sample = random.sample(domain,1)
        return sample[0][0]

def DEwhole(domain,data,epsilon):
	#Using the corresponding epsilon of DE 
    values = list()
    newdomain = []
    for item in domain:
        newdomain.append(item[0])
    #compute p, q with epsilon, d domain size
    d = 372
    p = (np.e)**epsilon/((np.e)**epsilon+d-1)
    #each user generates a value. preserve it with possibility p, change it to other value with possibility q
    q = 1/((np.e)**epsilon+d-1)
    for i in range(len(data)):
        values.append(DE(domain, data[i][0],p,q))
    fres = list([0 for i in range(372)])
    #the server computes the number of results of each item
    for i in range(len(values)):
        index = newdomain.index(values[i])
        fres[index] += 1
       
    #compute the unbaised report
    for i in range(len(fres)):
        fres[i] = (fres[i]- len(values)*q)/(p-q)  
    return fres

def SRR(user, GLCP, Pro, domain, m):
    r=random.uniform(0, 1)
    s = 0
    j = 0
    while j<m:
        s = s + Pro[j]*len(GLCP[j])
        if r <= s:
            sample = random.sample(GLCP[j],1)
            return sample
        else:
            j=j+1            

def SRRwhole(domain,data, GLCP, epsilon, size):
    newdomain = []
    for item in domain:
        newdomain.append(item[0])
    Pro_list = readall('length.txt')
    index_list = []
    Fpro_list = []
    for item in Pro_list:
        index_list.append(item[0])
        Fpro_list.append([item[1],item[2],item[3]])
    values = list()
    #compute probability
    m = 3 #revise the parameter m value 
    All_Pro = {}
    for i in range(len(data)):
        Pro =[]
        R = 0
        index = index_list.index(data[i][0])
        Fpro = Fpro_list[index]

        for j in range(2,m+1):
            R = R + (j-1)*float(Fpro[j-1])
        q = (m-1)/((m-1)*size*math.exp(epsilon)-(math.exp(epsilon)-1)*R)
        dif = q*(math.exp(epsilon)-1)/(m-1)
        for j in range(1,m+1):
            Pro.append(q+(m-j)*dif)
        if data[i][0] not in All_Pro:
            All_Pro[data[i][0]] = Pro
        sample=SRR(data[i], GLCP[data[i][0]], Pro, newdomain, m)
        values.append(sample)  
    #construct the transform
    transform = np.ones((size,size))
    for i in range(size):
        templist = GLCP[newdomain[i]]
        for j in range(m):
            grouplist = templist[j]
            for item in grouplist:
                tempindex = newdomain.index(item)
                transform[i,tempindex] = Pro[j]
    fres = list([0 for i in range(size)])
    #the server computes the number of results of each item
    for i in range(len(values)):
        index = newdomain.index(values[i][0])
        fres[index] += 1
    newfre = EM(newdomain, data, transform, fres)
    return newfre


def com_truefre(domain, data):
    newdomain = []
    fres = list([0 for i in range(372)])
    for item in domain:
        newdomain.append(item[0])
    for i in range(len(data)):
        index = newdomain.index(data[i][0])
        fres[index] += 1
    print("truefre:", fres)
    return fres

def l1(fre1,fre2):  
    fre = list()
    for i in range(len(fre1)):
        fre.append(np.linalg.norm(fre1[i]-fre2[i]))
    return sum(fre)

def EM(domain, data, transform, fre):
    #generate the C_x for each input
    Hadamard_matrix = functions.Hadarmard_init(int(2**9))
    index = 1
    
    Y = []
    for i in range(len(domain)):
        print(i)
        total = 0
        hl = Hadamard_matrix[i+1]
        key = domain[i]
        value = []
        for j in range(len(domain)):
            if hl[j] is True:
                value.append(domain[j])
                total = total + fre[j]
        Y.append(total)

        #Construct the linear equations
        X = np.ones((372, 372))

        for i in range(len(domain)):
            temp = domain[i]
            hl = Hadamard_matrix[i+1]
            for j in range(len(domain)):
                p = 0
                if i == j:               
                    for index in range(len(domain)):
                        if hl[index] is True:
                            p = p + transform[i,index]
                else:
                    for index in range(len(domain)):
                        if hl[index] is True:
                            p = p + transform[j,index]
                X[i,j] = p
    fre3 = np.linalg.inv(X).dot(Y)
    for i in range(len(fre)):
    	if fre3[i]<0:
    		fre3[i] = 0
    summ = sum(fre3)
    for i in range(len(fre3)):
        fre3[i] = fre3[i]/summ*len(data)
    return fre3

def EML(n, ns_hist, transform, max_iteration=1000, loglikelihood_threshold=3):
    theta = np.ones(n) / float(n)
    theta_old = np.zeros(n)
    r = 0
    sample_size = n
    old_logliklihood = 0

    while LA.norm(theta_old - theta, ord=1) > 1 / sample_size and r < max_iteration:
        theta_old = np.copy(theta)
        X_condition = np.matmul(transform, theta_old)
       
        TMP = transform.T / X_condition

        P = np.copy(np.matmul(TMP, ns_hist))
        P = P * theta_old

        theta = np.copy(P / sum(P))

        logliklihood = np.inner(ns_hist, np.log(np.matmul(transform, theta)))
        imporve = logliklihood - old_logliklihood

        if r > 1 and abs(imporve) < loglikelihood_threshold:
            break

        old_logliklihood = logliklihood

        r += 1
    return theta

'''
Example of GLCP format
GLCP = {'0111100100000110010001010100101001110000110011': [['0111100100000110010001010100101001110000110011'], ['0000000100000100111011111111000110101001010100','0111100100000100111011111100100011111010001111'], ['0111100100000100111011111100100011111010000010']], 
'0111100100000100111011111100100011111010000010': [['0111100100000100111011111100100011111010000010'], ['0000000100000100111011111111000110101001010100','0111100100000100111011111100100011111010001111'], ['0111100100000110010001010100101001110000110011']], 
'0000000100000100111011111111000110101001010100': [['0000000100000100111011111111000110101001010100'], ['0000000100000100111011111111000110101001010100'], ['0000000100000100111011111111000110101001010100','0111100100000100111011111100100011111010001111']],
'0111100100000100111011111100100011111010001111': [['0111100100000100111011111100100011111010001111'], ['0000000100000100111011111111000110101001010100'], ['0000000100000100111011111111000110101001010100','0111100100000100111011111100100011111010001111']]}
'''
