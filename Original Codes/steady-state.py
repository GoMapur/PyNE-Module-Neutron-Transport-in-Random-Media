# %THIS CODE SOLVES THE STEADY-STATE, MONOENERGETIC TRANSPORT EQUATION 
# %IN A HOMOGENEOUS MEDIUM WITH ISOTROPIC SCATTERING AND ISOTROPIC SOURCE
# %IN ROD AND SLAB (SN) GEOMETRIES. IT USES 2-POINT CENTRAL DIFFERENCING
# %(ORDER H^2) WITH 3-POINT FORWARD/BACKWARD AT THE BOUNDARIES.

# %MAIN OUTPUTS: transmission "TR"; reflection "RL"
# %			  scalar flux "SCAL"

# %RICHARD VASQUES & NITIN KUMAR YADAV

import numpy as np
import matplotlib.pyplot as plt
import os, sys
import random

def config_file(): #only supports homogenous right now. When refactoring code to clean out repeat functions, make generic
	folder = os.path.join(os.getcwd(), os.path.dirname(__file__))
	config = os.path.join(folder, "steady-state.config") 
	result = []
	#add error checking later
	with open(config, 'rb') as f:
		for line in f:
			line = line.strip()
			if len(line) == 0:
				continue
			if line[0] == "%":
				continue
			line = line[1:len(line)-1].split(",")
			if len(line) == 9:
				if line[0]:
					if line[0] == "SLAB":
						rod_slab = 1
					rod_slab = 0
				if line[1]:
					if rod_slab:
						N = int(line[1])
					else:
						N = 2
				if line[2]:
					T = float(line[2])
				if line[3]:
					n = int(line[3])
				if line[4]:
					et = float(line[4])
				if line[5]:
					cc = float(line[5])
					es = cc * et
				if line[6]:
					yo = float(line[6])
				if line[7]:
					y_ = float(line[7])
				if line[8]:
					Q = float(line[8])
			elif len(line) == 12:
				if line[0]:
					T = float(line[0])
				if line[1]:
					n = float(line[1])
				if line[2]:
					yo = float(line[2])
				if line[3]:
					y_ = float(line[3])
				if line[4]:
					m1 = float(line[4])
				if line[5]:
					m2 = float(line[5])
				if line[6]:
					Et1 = float(line[6])
				if line[7]:
					cc = float(line[7])
				if line[8]:
					Q1 = float(line[8])
				if line[9]:
					Et2 = float(line[9])
				if line[10]:
					cc2 = float(line[10])
				if line[11]:
					Q2 = float(line[11])
				print [T, n, yo, y_, m1, m2, Et1, cc, Q1, Et2, cc2, Q2]
			else:
				print "Malformed input data, please check."
				break
			result.append([rod_slab, N, T, n, et, cc, es, yo, y_, Q])
	return result

def homogenous_config_manual():
	print 'Default is ROD geometry'
	rod_slab = int(raw_input('Enter 1 if you want to change to SLAB Geometry:	'))

	if rod_slab == 1:
		N = 1
		while N % 2 == 1: # N has to be even
			N = int(raw_input('Enter the number of Discrete ordinate directions:	'))
		print 'Solving problem in SLAB Geometry'
		print 
	else:
		rod_slab = 0
		print 'Solving problem in ROD Geometry'
		print 
		N = 2
	T = float(raw_input('Enter the total length of the system:	'))
	n = 1
	while n < 4:
		n = int(raw_input('Enter at how many points you want to calculate:	'))
	et = float(raw_input('Enter the total cross section Sigma_t:	'))
	cc = 10
	while (cc > 1) or (cc < 0):
		cc = float(raw_input('Enter the Scattering Ratio c (between 0 and 1):	'))
	es = cc * et
	yo = float(raw_input('Enter the boundary value in the positive direction:	'))
	y_ = float(raw_input('Enter the boundary value in the negative direction:	'))
	Q = float(raw_input('Enter the homogeneous isotropic source:	'))
	return [[rod_slab, N, T, n, et, cc, es, yo, y_, Q]]

def frange(x, y, jump):
	result = []
	while x <= y:
		result.append(x)
		x += jump
	result.append(y)
	return result

def plot(x, SCAL, save, name):
	try:
		if len(x) > SCAL.shape[0]:
			x.pop(-1)
		
		plt.plot(x, SCAL, 'b')
		if save:
			plt.savefig(name +".png")
		else:
			plt.show()
		plt.clf()
	except Exception as e:
		print e
		print len(x), SCAL.shape
	return
	# not compatible with virtualenv

def save_csv(matrix, name):
	np.savetxt("csv/" + name + ".csv", matrix, delimiter=",")

def homogeneous(var, save=False):
	rod_slab = var[0]
	N = var[1]
	T = var[2]
	n = var[3]
	et = var[4]
	cc = var[5]
	es = var[6]
	yo = var[7]
	y_ = var[8]
	Q = var[9]

	M = n*N
	h = T/n

	A = np.zeros((M,M))
	B = np.full((1, M), Q/2)

	# % GAUSS-LEGENDRE QUADRATURE

	beta = range(1,N)
	beta = np.divide(beta, np.sqrt(np.subtract(np.multiply(np.power((range(1,N)), 2), 4), 1) ))
	#power 2, times 4, minus 1, sqrt each term, element wise divide
	#should it return a vector or a value? When testing with random numbers in Matlab, it looks like value. But steps after imply that its a vector
	x, w = np.linalg.eig(np.add(np.diag(beta, -1), np.diag(beta, 1))) #check indexing later
	u = x #already in the right form
	wt = np.multiply(np.power(w[0,:].T, 2), 2)
	#first row (index 0), transpose, square, multiply
	#central differences
	#finite volume
	#finite differences
	#things to search
	if rod_slab != 1:
		u[0] = -1
		u[1] = 1

	# % Diagonal Block of matrix up to N/2
	for t in range(1,(N/2) + 1):
		s = ((t - 1) * n) - 1 # minus one to account for indexing
		ut = u[t-1]
		A[s+1, s+1] = (-11 * ut) / (6*h) + et - es * wt[t-1] / 2
		A[s+1, s+2] = 3 * ut / h 
		A[s+1, s+3] = -3 * ut / (2*h)
		A[s+1, s+4] = ut/(3*h)
		
		for i in range(2, n):
			A[s+i, s+i-1] = -ut/(2*h)
			A[s+i, s+i] = (et-es*wt[t-1]/2)
			A[s+i, s+i+1] = ut/(2*h)
		A[s+n, s+n-1] = -ut	/ (2*h)
		A[s+n, s+n] = (et-es*wt[t-1]/2)
		B[0, s+n] = -ut * y_ / (2*h) + Q/2
		# % Remaining Blocks in same direction up to N/2
		l = t
		if (l == 1) and (N > 2):
			for p in range(l+1, (N/2) + 1):
				S = ((p-1) * n) - 1
				for i in range(1, n+1):
					A[s+i, S+i] = -es * wt[p-1] / 2
		elif (l > 1) and (N > 2):
			for p in range(1,l): 
				S = ((p-1) * n) - 1
				for i in range(1, n+1):
					A[s+i, S+i] = -es * wt[p-1] / 2
			for p in range(l+1, (N/2) + 1):
				S = ((p-1) * n) - 1
				for i in range(1, n+1):
					A[s+i, S+i] = -es * wt[p-1] / 2
		# % Blocks from N/2 to N
		a = 0
		for p in range(1, (N/2) + 1):
			S = ((N/2 + p - 1) * n) - 1
			for i in range(2, n+1):
				A[s+i, S+i-1] = -es * wt[(N/2+p) - 1] / 2
			a = a + (es * wt[(N/2+p) - 1] * yo/2)
		B[0,s+1] = a + Q/2
	# % Diagonal Block of matrix from N/2+1 to N
	for t in range(N/2 + 1, N+1):
		s = ((t-1) * n) - 1
		A[s+1, s+1] = (et-es*wt[t-1]/2)
		A[s+1, s+1+1] = u[t-1] / (2 * h)
		B[0,s+1] = u[t-1] * yo / (2 * h) + Q/2
		for i in range(2, (n-1) + 1):
			A[s+i, s+i-1] = -u[t-1] / (2*h)
			A[s+i, s+i] = (et-es*wt[t-1]/2)
			A[s+i, s+i+1] = u[t-1]/(2*h)

		A[s+n, s+n-3] = -u[t-1] / (3*h)
		A[s+n, s+n-2] = 3 * u[t-1] / (2*h)
		A[s+n, s+n-1] = -3 * u[t-1] / h;
		A[s+n, s+n] = (11 * u[t-1] / (6*h) + et - es * wt[t-1] / 2)
		# % Remaining Blocks in same direction up to N
		l = t
		if (l==N/2 + 1) and (N>2):
			for p in range(l+1, N+1):
				S = ((p-1) * n) - 1
				for i in range(1, n+1):
					A[s+i, S+i] = -es * wt[p-1] / 2
		elif (1 > N/2 + 1) and (N>2):
			for p in range(N/2 + 1, (l-1) + 1):
				S = ((p-1) * n) - 1
				for i in range(1, n+1):
					A[s+i, S+i] = -es * wt[p-1] / 2
			for p in range(l+1, N+1):
				S = ((p-1) * n) - 1
				for i in range(1, n+1):
					A[s+i, S+i] = -es * wt[p-1] / 2
		# % Blocks from 1 to N/2
		a = 0
		for p in range(1, (N/2) + 1):
			S = ((p-1) * n) - 1
			for i in range(1, (n-1) + 1):
				A[s+i, S+i+1] = -es * wt[p-1] / 2
			a = a + (es * wt[p-1] * y_ / 2)
		B[0,s+n] = a + Q / 2

	save_csv(A, "A")
	save_csv(B, "B")
	#may be some accuracy solving issues algorithm wise with numpy 
	# refer to http://stackoverflow.com/questions/25001753/numpys-linalg-solve-and-linalg-lstsq-not-giving-same-answer-as-matlabs
	# below are different attempts to see if any difference RESULT: no difference
	# AX = B.T
	# A1 = np.linalg.inv(A)
	# np.savetxt("Ainv.csv", A1, delimiter=",")
	# X = np.dot(A1, B.T)
	# np.savetxt("B.csv", B.T, delimiter=",")
	# np.savetxt("X_manual.csv", X, delimiter=",")
	# np.savetxt("X_solve.csv", X1, delimiter=",")
	# X2 = np.linalg.lstsq(A, B.T)[0]
	# np.savetxt("X_lstsq.csv", X2, delimiter=",")
	# print A1 == A
	# print np.linalg.cond(A1)

	X = np.linalg.solve(A, B.T) #checked values compared to the other methods of solving and it seems to be accurate enough

	Y = np.zeros((M+N,1))

	# % Adding boundary conditions to the array
	i = 1
	for t in range(1, (N/2) + 1):
		s = (t-1)
		for j in range(i, (t*n+s) + 1): #might have conflict here with the s-1
			Y[j-1] = X[(j-s)-1]
		Y[(j+1) - 1] = y_
		i = j + 2
	i -= 1
	for t in range((N/2) + 1, N+1):
		Y[(i+1) -1] = yo
		i += 1
		for j in range(i+1, (t*n+t) + 1):
			Y[j - 1] = X[(j-t) - 1]
		i = j

	RL = 0

	for t in range(1, (N/2) + 1):
		s = (t-1) * n 
		RL = RL + abs(wt[t-1] * u[t-1] * Y[(s+t) - 1]) # might have conflict here
	TR = 0
	for t in range((N/2) + 1, N + 1):
		s = (t-1) * n
		S = (t - N/2 - 1) * n
		k = N/2 + 1
		TR = TR + wt[t-1] * u[t-1] * Y[(s+n+t) - 1] # might have a conflict here
	m = n + 1
	SCAL = np.zeros((m, 1))
	for t in range(1, N + 1):
		s = (t-1) * m
		for i in range(1, m + 1):
			SCAL[i-1] = SCAL[i-1] + wt[t-1] * Y[(s+i) - 1]

	save_csv(Y, "Y")
	save_csv(SCAL, "SCAL")

	x = frange(0, T, h)

	#forming name for the file save
	name = str(rod_slab) + "_" + str(T) + "_" + str(N) + "_" + str(n) + "_" + str(et) + "_" + str(cc) + "_" + str(es) + "_" + str(yo) + "_" + str(y_) + "_" + str(Q)
	plot(x, SCAL, save, name)

#extra order does matter
#check h ascending order

def SN_per_bench_solver(T,m1,m2,n,N,Es1,Es2,Et1,Et2,yo,y_,Q1,Q2,u,wt,a):
	interval = T / n
	m12 = [0, m1, m2] #which material, so index of this matters
	mm = 0
	s = 0
	i = 0
	x = None
	if a > 1:
		if a < m2/interval + 2:
			i += 1
			x1 = (a-1) * interval
			s += x1
			tmp = [] #can shorten to declare immediately but verbose for readability
			tmp.append(s)
			tmp.append(2)
			x = np.asarray([tmp]) #x is offset to account for material 
		else:
			i += 1
			x1 = (a - m2/interval - 1) * interval
			s += x1;
			tmp = []
			tmp.append(s)
			tmp.append((mm % 2) + 1)
			x = np.asarray([tmp])
			mm += 1
	while s < T:
		i += 1
		x1 = m12[(mm % 2) + 1] #whats the point of this? isnt this either 1 or 2?
		s = s + x1
		tmp = []
		if s <= T:
			tmp.append(s)
		else:
			tmp.append(T) #check here
		tmp.append((mm % 2) + 1)
		if x == None:
			x = np.asarray([tmp])
		else: #check logic here
			x = np.vstack((x, np.asarray([tmp])))
		mm += 1 
	H = T / n
	n1 = 1
	i = 1
	j = 1
	extra = [0] * n
	t1 = i * H

	save_csv(x, "xsn")


	tmp = [0, int(x[j-1 , 2-1])]
	L = np.asarray([tmp])
	L1 = []
	L1.append(tmp)
	h = [] #is h just a list?
	if t1 == x[j-1, 1-1]:
		extra[i-1] = 1
		h.append(x[j-1, 1-1] / 2)
		#n1+1 case
		tmp = [h[n1-1], int(x[j-1, 2-1])]
		L1.append(tmp)
		tmp = np.asarray([tmp])
		L = np.vstack((L, tmp))
		h.append(h[n1-1])
		#n1+2 case
		tmp = np.asarray([x[j-1,1-1], 3]) #are the values of x guaranteed to be ints???
		L1.append([x[j-1,1-1], 3])
		L = np.vstack((L, tmp))
		n1 += 2
		i += 1
		t1 = i * H
	else:
		while t1 <= x[j-1,1-1]:
			h.append(H) #h(n1) = H
			tmp = [i * H]
			if tmp[0] == x[j-1,1-1]:
				tmp.append(3)
			else:
				tmp.append(int(x[j-1,2-1]))
			L1.append(tmp)
			L = np.vstack((L, np.asarray([tmp])))
			n1 += 1
			i += 1
			t1 = i * H
	j = 2
	if x[1-1,1-1] != T:
		while i <= n:
			while t1 <= x[j-1,1-1]:
				tmp = [i * H]
				if tmp[0] == x[j-1,1-1]:
					tmp.append(3)
				else:
					tmp.append(int(x[j-1,2-1]))
				L1.append(tmp)
				L = np.vstack((L, np.asarray([tmp])))
				h.append(L1[n1+1-1][1-1] - L1[n1-1][1-1])
				n1 += 1 #n1 isnt really necessary to keep the index
				i += 1
				t1 = i * H
			j += 1
			if L1[n1-1][1-1] == T and L1[n1-1][2-1] == 3 and L1[n1-1-1][2-1] == 3:
				i -= 1
				extra[i-1] += 1
				tmp = []
				tmp.append((x[j-1-1, 1-1] + x[j-2-1, 1-1])/2)
				tmp.append(x[j-1-1, 2-1])
				tmp[1] = int(tmp[1])
				L = np.vstack((L, np.asarray([tmp])))
				L1[n1-1] = tmp
				h[n1-1-1] = L1[n1-1][1-1] - L1[n1-1-1][1-1] #still within range
				n1 += 1
				tmp = [x[j-1-1, 1-1], 3]
				tmp[1] = int(tmp[1])
				L1.append(tmp)
				L = np.vstack((L, np.asarray([tmp])))
				h.append(L1[n1-1][1-1] - L1[n1-1-1][1-1]) #adds to h[] (Matlab) not indexing
				i += 1
	save_csv(L1, "L")
	n1 -= 1
	L = L1
	save_csv(np.asarray(h), "h")
	#L tells you which material, L's index doesnt matter, x's index doesnt matter, Es, Et, Q does
	A = np.zeros((N*n1,N*n1))
	B = np.zeros((N*n1,1))
	Et12 = [0, Et1, Et2]
	Es12 = [0, Es1, Es2] 
	Q12 = [0, Q1, Q2]


	#% Diagonal Block of matrix up to N/2.................................
	for t in range(1, N/2 + 1):
		s = (t-1) * n1
		A[s+1-1,s+1-1] = -u[t-1] * (1 / h[1-1]) + Et12[L[1-1][2-1]] - Es12[L[1-1][2-1]] * wt[t-1]/2
		A[s+1-1,s+2-1] = u[t-1] * (1 / h[1-1])
		B[s+1-1, 0] = Q12[L[1-1][2-1]]/2
		for i in range(2, n1-1+1):
			if L[i-1][2-1] == 3:
				if i == n1 - 1:
					A[s+i-1,s+i-1] = -u[t-1] * (1 / h[i-1] + 1 / (h[i-1] + h[i+1-1])) + Et12[L[i+1-1][2-1]] - Es12[L[i+1-1][2-1]] * wt[t-1]/2
					A[s+i-1,s+i+1-1] = u[t-1] * (1/h[i-1] + 1/h[i+1-1])
					B[s+i-1,0] = u[t-1] * y_ * (h[i-1] / (h[i+1-1] * (h[i-1] + h[i+1-1]))) + Q12[L[i+1-1][2-1]]/2
				else:
					A[s+i-1,s+i-1] = -u[t-1] * (1/h[i-1] + 1/(h[i-1] + h[i+1-1] )) + Et12[L[i+1-1][2-1]] - Es12[L[i+1-1][2-1]] * wt[t-1]/2
					A[s+i-1,s+i+1-1] = u[t-1] * (1/h[i-1] + 1/h[i+1-1])
					A[s+i-1,s+i+2-1] = -u[t-1] * (h[i-1]/(h[i+1-1] * (h[i-1] + h[i+1-1])))
					B[s+i-1, 0] = Q12[L[i+1-1][2-1]]/2
			else:
				A[s+i-1,s+i-1-1] = -u[t-1] * h[i-1]/(h[i-1-1]*(h[i-1-1]+h[i-1]))
				A[s+i-1,s+i-1] = u[t-1] * (h[i-1] - h[i-1-1])/(h[i-1] * h[i-1-1]) + Et12[L[i-1][2-1]] - Es12[L[i-1][2-1]] * wt[t-1]/2
				A[s+i-1,s+i+1-1] = u[t-1] * h[i-1-1] / (h[i-1] * (h[i-1-1] + h[i-1]))
				B[s+i-1,0] = Q12[L[i-1][2-1]] /2
		A[s+n1-1,s+n1-1-1] = -u[t-1] * h[n1-1] / (h[n1-1-1] * (h[n1-1-1] + h[n1-1]))
		A[s+n1-1,s+n1-1] = u[t-1] * (h[n1-1] - h[n1-1-1]) / (h[n1-1] * h[n1-1-1]) + Et12[L[n1-1][2-1]] - Es12[L[n1-1][2-1]] * wt[t-1]/2
		B[s+n1-1,0] = -u[t-1] * y_ * h[n1-1-1] / (h[n1-1] * (h[n1-1-1] + h[n1-1] )) + Q12[L[n1-1][2-1]]/2


		#% Blocks from N/2 to N........................................
		a = 0
		for p in range(1, N/2 + 1):
			S = (N/2 + p - 1) * n1
			for i in range(2, n1+1):
				if L[i-1][2-1] == 3:
					A[s+i-1,S+i-1-1] = -Es12[L[i+1-1][2-1]] * wt[N/2+p-1]/2
				else:
					A[s+i-1,S+i-1-1] = -Es12[L[i-1][2-1]] * wt[N/2+p -1]/2
			a += (Es12[L[1-1][2-1]] * wt[N/2+p-1] * yo/2)
		B[s+1-1,0] = B[s+1-1,0] + a
		
	#% Diagonal Block of matrix from N/2+1 to N.........................
	for t in range(N/2+1, N+1):
		s = (t-1) * n1
		A[s+1-1,s+1-1] = u[t-1] * (h[2-1] - h[1-1]) / (h[2-1] * h[1-1]) + Et12[L[1-1][2-1]] - Es12[L[1-1][2-1]] * wt[t-1]/2
		A[s+1-1,s+1+1-1] = u[t-1] * h[1-1] / (h[2-1] * (h[1-1] + h[2-1]))
		B[s+1-1,0] = u[t-1] * yo * h[2-1] / (h[1-1] * (h[1-1] + h[2-1])) + Q12[L[1-1][2-1]] / 2
		for i in range(2, n1-1+1):
			if L[i+1-1][2-1] == 3:
				if i == 2:
					A[s+i-1,s+i-1] = u[t-1] * (1/h[i-1] + 1/(h[i-1] + h[i-1-1])) + Et12[L[i-1][2-1]] - Es12[L[i-1][2-1]] * wt[t-1]/2
					A[s+i-1,s+i-1-1] = -u[t-1] * (1/h[i-1] + 1/h[i-1-1])
					B[s+i-1,0] = -u[t-1] * yo * (h[i-1]/(h[i-1-1] * (h[i-1] + h[i-1-1]))) + Q12[L[i-1][2-1]]/2
				else:
					A[s+i-1,s+i-1] = u[t-1] * (1/h[i-1] + 1/(h[i-1] + h[i-1-1])) + Et12[L[i-1][2-1]] - Es12[L[i-1][2-1]] * wt[t-1]/2
					A[s+i-1,s+i-1-1] = -u[t-1] * (1/h[i-1] + 1/h[i-1-1])
					A[s+i-1,s+i-2-1] = u[t-1] * (h[i-1] / (h[i-1-1] * (h[i-1] + h[i-1-1])))
					B[s+i-1,0] = Q12[L[i-1][2-1]]/2
			else:
				A[s+i-1,s+i-1-1] = -u[t-1] * h[i+1-1] / (h[i-1] * (h[i-1] + h[i+1-1]))
				A[s+i-1,s+i-1] = u[t-1] * (h[i+1-1] - h[i-1]) / (h[i+1-1] * h[i-1]) + Et12[L[i+1-1][2-1]] - Es12[L[i+1-1][2-1]] * wt[t-1]/2
				A[s+i-1,s+i+1-1] = u[t-1] * h[i-1] / (h[i+1-1] * (h[i-1] + h[i+1-1]))
				B[s+i-1,0] = Q12[L[i+1-1][2-1]]/2
		A[s+n1-1,s+n1-1] = u[t-1] * (1/h[n1-1]) + Et12[L[n1-1][2-1]] - Es12[L[n1-1][2-1]] * wt[t-1]/2
		A[s+n1-1,s+n1-1-1] = -u[t-1] * (1/h[n1-1])
		B[s+n1-1,0] = Q12[L[n1-1][2-1]]/2

		
		#% Blocks from 1 to N/2........................................
		a = 0
		for p in range(1, N/2+1):
			S = (p-1) * n1
			for i in range(1,n1-1+1):
				if L[i+1-1][2-1] == 3: 
					A[s+i-1,S+i+1-1] = -Es12[L[i-1][2-1]] * wt[p-1]/2
				else:
					A[s+i-1,S+i+1-1] = -Es12[L[i+1-1][2-1]] * wt[p-1]/2
			a += (Es12[L[n1-1][2-1]] * wt[p-1] * y_/2)

		B[s+n1-1, 0] += a
	Z = np.linalg.solve(A, B)
	return Z, n1, B, L, A, extra

def periodic_config_manual(): #redundant, clean up to make it generic later
	print 'Default is ROD geometry'
	rod_slab = int(raw_input('Enter 1 if you want to change to SLAB Geometry:	'))

	if rod_slab == 1:
		N = 1
		while N % 2 == 1: # N has to be even
			N = int(raw_input('Enter the number of Discrete ordinate directions:	'))
		print 'Solving problem in SLAB Geometry'
		print 
	else:
		rod_slab = 0
		print 'Solving problem in ROD Geometry'
		print 
		N = 2
	T = float(raw_input('Enter the total length of the system:	'))
	n = 1
	while n < 4:
		n = int(raw_input('Enter at how many points you want to calculate:	'))
	et = float(raw_input('Enter the total cross section Sigma_t:	'))
	cc = 10
	while (cc > 1) or (cc < 0):
		cc = float(raw_input('Enter the Scattering Ratio c (between 0 and 1):	'))
	es = cc * et
	yo = float(raw_input('Enter the boundary value in the positive direction:	'))
	y_ = float(raw_input('Enter the boundary value in the negative direction:	'))
	Q = float(raw_input('Enter the homogeneous isotropic source:	'))
	return [[rod_slab, N, T, n, et, cc, es, yo, y_, Q]]

def periodic(var):
	#T and n have to be floats for math
	rod_slab = var[0]
	N = var[1] #has to be even add checking functions to config file function later
	T = var[2]
	n = var[3]
	yo = var[4]
	y_ = var[5]
	m1 = var[6]
	m2 = var[7]
	Et1 = var[8]
	cc = var[9]
	Q1 = var[10]
	Et2 = var[11]
	cc2 = var[12]
	Q2 = var[13]

	beta = range(1,N)
	beta = np.divide(beta, np.sqrt(np.subtract(np.multiply(np.power((range(1,N)), 2), 4), 1) ))
	x, w = np.linalg.eig(np.add(np.diag(beta, -1), np.diag(beta, 1))) #check indexing later
	u = x #already in the right form
	wt = np.multiply(np.power(w[0,:].T, 2), 2)
	if rod_slab != 1:
		u[0] = -1
		u[1] = 1

	#ROD case
	# N = 2
	# wt = [1, 1]
	# u = [-1, 1]
	#Add ROD SLAB code to deal with weights, weights is 1xn mat

	Es1 = cc * Et1
	Ea1 = Et1 - Es1
	Es2 = cc2 * Et2
	Ea2 = Et2 - Es2

	a = 1
	reflec = 0
	reflec2 = 0
	transm = 0
	transm2 = 0
	SF = np.zeros((n+1,1))
	SF2 = np.zeros((n+1,1))
	cond = 0

	total = (m1+m2) / (T/n)

	while (a < (total+1)):
		print "Problem " + str(a) + " of " + str(total)
		Z,n1,B,L,A, extra = SN_per_bench_solver(T,m1,m2,n,N,Es1,Es2,Et1,Et2,yo,y_,Q1,Q2,u,wt,a)
		save_csv(np.asarray(extra), "extra")
		save_csv(L, "L")
		save_csv(A, "A")
		save_csv(B, "B")
		save_csv(Z, "Z")
		#% Adjusting points..........................
		X = np.zeros((n*N,1))
		for i in range(1, (N/2) + 1):
			X[(i-1)*n+1 - 1] = Z[(i-1)*n1+1 - 1]
			k = 2  #not sure why this is needed. 
			for j in range(2, n + 1):
				k += extra[j-1-1]
				# k is correct
				#there is some trickery going on with k.
				X[(i-1)*n+j - 1] = Z[(i-1)*n1+k - 1]
				k += 1
		for i in range(N/2 + 1, N + 1):
			k = 1 
			for j in range(1, n + 1):
				k += extra[j-1]
				X[(i-1)*n+j - 1] = Z[(i-1)*n1+k - 1]
				k += 1
		
		Y = np.zeros((N*(n+1),1));
		#% Adding boundary conditions...............
		i = 1
		for t in range(1, N/2 + 1):
			s = t-1
			for j in range(i, t*n+s + 1):
				Y[j - 1] = X[j-s - 1]
			Y[j+1 - 1] = y_
			i = j + 2
		i -= 1
		for t in range(N/2+1, N+1):
			Y[i+1 - 1] = yo
			i += 1
			for j in range(i+1, t*n+t + 1):
				Y[j - 1] = X[j-t - 1]
			i = j

		save_csv(X, "X")
		save_csv(Y, "Y")
		#% Calculating Reflection, Transmission and Scalar Flux......	 
		RL = 0
		for t in range(1, N/2 + 1):
			s = (t-1) * n
			RL += abs(wt[t - 1] * u[t - 1] * Y[s+t - 1])
		TR = 0
		for t in range(N/2 + 1, N + 1):
			s = (t-1) * n
			TR += wt[t - 1] * u[t - 1] * Y[s+n+t - 1]
		m = n + 1
		SCAL = np.zeros((m,1))
		for t in range(1, N + 1):
			s = (t-1) * m
			for i in range(1, m + 1):
				SCAL[i - 1] = SCAL[i - 1] + wt[t -1] * Y[s+i - 1]
		
		reflec += RL
		reflec2 += (RL)**2
		transm += TR
		transm2 += (TR)**2
		for i in range(1, n+1 + 1):
			SF[i - 1] = SF[i - 1] + SCAL[i - 1]
			SF2[i - 1] = SF2[i - 1] + (SCAL[i - 1])**2
		a += 1
	a -= 1
	SF = np.divide(SF, a)
	SF2 = np.divide(SF2, a)

	kk = T/n
	#% p is line along x-axis.
	p = frange(0, T, kk)
	p.pop()

	save_csv(SF, "SF")
	plt.plot(p,SF,'r')
	plt.plot(p,SF2,'b')
	plt.show()

def non_homogenous(var):
	rod_slab = var[0]
	N = var[1]
	T = var[2]
	m1 = var[3]
	m2 = var[4]
	n = var[5]
	N = var[6]
	Es1 = var[7]
	Es2 = var[8]
	Et1 = var[9]
	Et2 = var[10]
	yo = var[11]
	y_ = var[12]
	Q1 = var[13]
	Q2 = var[14]
	u = var[15]
	wt = var[16]
	randseed = var[17]
	med = var[18]
	a = var[19]

	# %Building realization----------
	# %EXP. Random Medium....

	#random python docs: https://docs.python.org/2/library/random.html
	if not med:
		random.seed(randseed)
		m12 = [0, m1, m2] #refers to which material, pad with 0 so index 1 == material 1 (slightly less confusing this way)
		xx == 0.0
		xx = random.random() #python random gives a float in [0.0, 1.0) 
		#turns out both Matlab and python use deterministic randomgenerators
		if xx <= m1/(m1+m2):
			mm = 0
		else:
			mm = 1
		s = 0
		i = 0
		x = None #is this a guaranteed loop?
		while s < T:
			i += 1
			x1 = random.expovariate(1/m12[mm%2 + 1]) #exprnd returns a 1x1
			#why does x1 matter? its not being used in the loop unless it changes the previous value for random need to check matlab docs
			s += x1
			tmp = []
			if s <= T:
				tmp.append(s)
			else:
				tmp.append(T)
			tmp.append(mm%2 + 1)
			if not x:
				x = np.asarray([tmp])
			else:
				x = np.vstack((x, np.asarray([tmp])))
			mm += 1
		randseed = random.getstate() #use setstate() later
	else:
		interval = T/n
		m12 = [0, m1, m2]
		mm = 9
		i = 0
		s = 0
		if a > 1:
			if a < m1/interval+2: #check order of op here
				i += 1
				x1 = (a-1) * interval
				s += x1
				x = np.asarray([[s, 2]])
			else: # not sure if this is inner or outer if CHECK 
				i += 1
				x1 = (a - m1/interval - 1) * interval
				s += x1
				x = np.asarray([s, mm%2 + 1])
				mm += 1
		while s < T:
			i += 1
			x1 = m12[mm%2 + 1]
			s += 1
			tmp = []
			if s <= T:
				tmp.append(s)
			else:
				tmp.append(T)
			tmp.append(mm%2 + 1)
			x = np.vstack((x, np.asarray([tmp])))
			mm += 1
	# %Adding interfaces and extra points inside layers--------------
	H = T/n
	n1 = 1
	i = 1
	j = 1
	extra = [] #using a list will be easier for extra follow periodic code
	h1 = [] # ditto comment above
	t1 = i * H
	L = np.asarray([[0, x[j-1, 2-1]]])
	if t1 > x[j-1, 1-1]:
		extra[i-1] = 2
		h[n1-1] = x[j-1, 1-1]/2
		L[n1+1-1,1-1] = h[n1-1]
		L[n1+1-1,2-1] = x[j-1, 2-1]
		h[n1+1-1]=h[n1-1]
		L[n1+2-1,1-1] = x[j-1, 1-1]
		L[n1+2-1,2-1] = 3
		n1 += 2
	else:
		while t1 <= x[j-1, 1-1]:
			h[n1-1] = H




def main(sysargs):
	# default state is file configuration
	
	# if not len(sysargs):
	# 	save = False
	# elif sysargs[0].lower() == "manual":
	# 	var = homogenous_config_manual()
	# elif sysargs[0].lower() == "save":
	# 	save = True
	# var = config_file()
	# for v in var:
	# 	homogenous(v, save)
	
	#lots of redudant code, need to make them generic later

	periodic([0, 2, 20.0, 100, 0.0, 0.0, 1.0, 1.0, 1.0, 0.5, 1.0, 0.0, 0.0, 0.0])

if __name__ == '__main__':
	main(sys.argv[1:])













