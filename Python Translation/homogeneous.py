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

	save_csv(A, 'A')
	save_csv(B, 'B')
	#may be some accuracy solving issues algorithm wise with numpy
	# refer to http://stackoverflow.com/questions/25001753/numpys-linalg-solve-and-linalg-lstsq-not-giving-same-answer-as-matlabs
	# below are different attempts to see if any difference RESULT: no difference
	# AX = B.T
	# A1 = np.linalg.inv(A)
	# np.savetxt('Ainv.csv', A1, delimiter=',')
	# X = np.dot(A1, B.T)
	# np.savetxt('B.csv', B.T, delimiter=',')
	# np.savetxt('X_manual.csv', X, delimiter=',')
	# np.savetxt('X_solve.csv', X1, delimiter=',')
	# X2 = np.linalg.lstsq(A, B.T)[0]
	# np.savetxt('X_lstsq.csv', X2, delimiter=',')
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

	save_csv(Y, 'Y')
	save_csv(SCAL, 'SCAL')

	x = frange(0, T, h)

	#forming name for the file save
	name = str(rod_slab) + '_' + str(T) + '_' + str(N) + '_' + str(n) + '_' + str(et) + '_' + str(cc) + '_' + str(es) + '_' + str(yo) + '_' + str(y_) + '_' + str(Q)
	plot(x, SCAL, save, name)
