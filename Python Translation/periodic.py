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
		print 'Problem ' + str(a) + ' of ' + str(total)
		Z,n1,B,L,A, extra = SN_per_bench_solver(T,m1,m2,n,N,Es1,Es2,Et1,Et2,yo,y_,Q1,Q2,u,wt,a)
		save_csv(np.asarray(extra), 'extra')
		save_csv(L, 'L')
		save_csv(A, 'A')
		save_csv(B, 'B')
		save_csv(Z, 'Z')
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

		save_csv(X, 'X')
		save_csv(Y, 'Y')
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

	save_csv(SF, 'SF')
	plt.plot(p,SF,'r')
	plt.plot(p,SF2,'b')
	plt.show()

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

	save_csv(x, 'xsn')


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
	save_csv(L1, 'L')
	n1 -= 1
	L = L1
	save_csv(np.asarray(h), 'h')
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
