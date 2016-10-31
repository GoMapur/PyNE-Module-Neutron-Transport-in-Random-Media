# Codes for Nuclear Engineering Department @ Berkeley
# This is a summer project, translation and packaging past codes for 1-D
# Nutron transport simulations.
# Note this file is highly extensible so I suggest it could be used in future
# similar simulation designs, plz improve anything if you'd like, it would be
# my honor.
# Mingjian Lu, July 2016

import numpy as np
import matplotlib.pyplot as plt
from OneDGridModel import *
from OneDGridModelSolver import *
import functools
import itertools
from decimal import *
getcontext().prec = 14

class Model_1D_Benchmark():
    def __init__(self, total_len, point_num, boundary_cond, materials, gauss_discrete_direction_num):
        self.total_len = total_len
        self.gauss_discrete_direction_num = gauss_discrete_direction_num
        self.point_num = point_num
        self.boundary_cond = boundary_cond
        self.materials = materials

    def benchmark(self):
        raise NotImplementedError

class Model_1D_Stochastic_Finite_Step_Benchmark(Model_1D_Benchmark):
    def __init__(self, total_len, point_num, boundary_cond, materials, gauss_discrete_direction_num = 2):
        Model_1D_Benchmark.__init__(self, float(total_len), point_num, boundary_cond, materials, gauss_discrete_direction_num)

    def benchmark_once(self):
        # NOTE: Gauss-Legendre is taken care of inside the solver instead of the
        #       benchmark, It has an internal cache to dealwith redundant
        #       calculation
        # TODO: Need to add boundary points to the solution to make it symmetric
        #       and complete.
        grid_model = Stochastic_Gird(self.total_len, self.boundary_cond, self.materials)
        solver = Model_1D_Stochastic_Finite_Step_Solver(grid_model, self.point_num, gauss_discrete_direction_num = self.gauss_discrete_direction_num)
        solution = solver.solve_scalar_flux()
        # solver.plot_scalar_flux()
        return solution

    def benchmark(self):
        # TODO: This fixed number should changed to be a indicator
        #       showing when should we stop
        iter_times = 1
        solution_list = list(itertools.starmap(self.benchmark_once, [()] * iter_times))
        sum_of_sol = functools.reduce( (lambda x, y: np.array(x) + np.array(y)), solution_list )
        avg_sol = np.array(sum_of_sol) / float(iter_times)
        plt.plot(frange(0, self.total_len, self.total_len / self.point_num), avg_sol)
        plt.show()


class Model_1D_Periodic_Finite_Step_Benchmark(Model_1D_Benchmark):
    def __init__(self, total_len, point_num, boundary_cond, materials, gauss_discrete_direction_num = 2):
        Model_1D_Benchmark.__init__(self, float(total_len), point_num, boundary_cond, materials, gauss_discrete_direction_num)

    def benchmark(self):
        grid_model = Periodic_Grid(self.total_len, self.boundary_cond, self.materials)
        N = self.gauss_discrete_direction_num
        T = self.total_len
        n = self.point_num
        yo = self.boundary_cond[0]
        y_ = self.boundary_cond[1]
        m1 = self.materials[0].thickness()
        m2 = self.materials[1].thickness()

        mat1, mat2 = grid_model.materials[0], grid_model.materials[1]
        T = grid_model.len()
        m1,m2 = mat1.thickness(), mat2.thickness()
        Es1,Es2 = mat1.scattering_section(), mat2.scattering_section()
        Et1,Et2 = mat1.cross_section(), mat2.cross_section()
        yo,y_ = grid_model.left_boundary_condition(), grid_model.right_boundary_condition()
        Q1,Q2 = mat1.source(), mat2.source()

        a = 1
        reflec = 0
        reflec2 = 0
        transm = 0
        transm2 = 0
        SF = np.zeros((n+1,1))
        SF2 = np.zeros((n+1,1))
        cond = 0

        total = (m1+m2) / (T/n)

        while (a < (total+2)):
            print('Problem ' + str(a) + ' of ' + str(total+1))
            solver = Model_1D_Periodic_Solver(grid_model, self.point_num, self.gauss_discrete_direction_num, a)
            u,wt = solver.gauss_u(), solver.gauss_weight()
            Z,n1,B,L,A, extra = solver.solve()
            # save_csv(np.asarray(extra), 'extra')
            # save_csv(L, 'L')
            np.savetxt("A_re"+str(a)+".csv", np.asarray(A))
            # save_csv(Z, 'Z')
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

            # save_csv(X, 'X')
            # save_csv(Y, 'Y')
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

        # save_csv(SF, 'SF')
        plt.plot(p,np.array(SCAL),'g')
        # plt.plot(p,SF,'r')
        # plt.plot(p,SF2,'b')
        plt.show()

class Model_1D_Homogeneous_Finite_Step_Benchmark(Model_1D_Benchmark):
    def __init__(self, total_len, point_num, boundary_cond, materials, gauss_discrete_direction_num = 2):
        Model_1D_Benchmark.__init__(self, float(total_len), point_num, boundary_cond, materials, gauss_discrete_direction_num)

    def benchmark(self):
        rod_slab = 1
        N = self.gauss_discrete_direction_num
        T = self.total_len
        n = self.point_num
        et = self.materials[0].cross_section()
        cc = self.materials[0].scattering_ratio()
        es = et * cc
        yo = self.boundary_cond[0]
        y_ = self.boundary_cond[1]
        Q = self.materials[0].source()

        M = n*N
        h = T/n

        A = np.zeros((M,M))
        B = np.full((1, M), Q/2)

        # % GAUSS-LEGENDRE QUADRATURE
        u = np.polynomial.legendre.leggauss(N)[0]
        wt = np.polynomial.legendre.leggauss(N)[1]

        #first row (index 0), transpose, square, multiply
        #central differences
        #finite volume
        #finite differences
        #things to search

        # TODO: seperate rod/slab init cases
        if rod_slab == 0 and N == 2:
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
            A[s+n, s+n-1] = -ut    / (2*h)
            A[s+n, s+n] = (et-es*wt[t-1]/2)
            B[0, s+n] = -ut * y_ / (2*h) + Q/2
            # % Remaining Blocks in same direction up to N/2
            l = t
            if (l == 1) and (N > 2):
                for p in range(l+1, (N/2) + 1):
                    S = ((p-1) * n) - 1
                    for i in range(1, n+1):
                        A[s+i, S+i] = -es * wt[p-1] / 2
                        # print str(0) + " == " + str(s+i) + " " + str(S+i) + " " + str(-es * wt[p-1] / 2)
            elif (l > 1) and (N > 2):
                for p in range(1,l):
                    S = ((p-1) * n) - 1
                    for i in range(1, n+1):
                        A[s+i, S+i] = -es * wt[p-1] / 2
                        # print str(1) + " == " + str(s+i) + " " + str(S+i) + " " + str(-es * wt[p-1] / 2)
                for p in range(l+1, (N/2) + 1):
                    S = ((p-1) * n) - 1
                    for i in range(1, n+1):
                        A[s+i, S+i] = -es * wt[p-1] / 2
                        # print str(2) + " == " + str(s+i) + " " + str(S+i) + " " + str(-es * wt[p-1] / 2)
            # % Blocks from N/2 to N
            a = 0
            for p in range(1, (N/2) + 1):
                S = ((N/2 + p - 1) * n) - 1
                for i in range(2, n+1):
                    A[s+i, S+i-1] = -es * wt[(N/2+p) - 1] / 2
                    # print str(3) + " == " + str(s+i) + " " + str(S+i-1) + " " + str(-es * wt[(N/2+p) - 1] / 2)
                    # print str(s+i) + " " + str(S+i-1) + " " + str(-es) + " " + str(N/2+p) + " " + str(wt[(N/2+p) - 1])
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
            elif (l > N/2 + 1) and (N>2):
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

        # save_csv(A, 'A')
        # save_csv(B, 'B')
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

        np.savetxt("A_h.csv", np.asarray(A))
        np.savetxt("B_h.csv", np.asarray(B.T))

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

        # save_csv(Y, 'Y')
        # save_csv(SCAL, 'SCAL')

        x = frange(0, T, h)

        #forming name for the file save
        name = str(rod_slab) + '_' + str(T) + '_' + str(N) + '_' + str(n) + '_' + str(et) + '_' + str(cc) + '_' + str(es) + '_' + str(yo) + '_' + str(y_) + '_' + str(Q)

        plt.plot(x, SCAL)
        plt.show()

def frange(start, end=None, inc=None):
    """A range function, that does accept float increments..."""
    import math

    if end == None:
        end = start + 0.0
        start = 0.0
    else: start += 0.0 # force it to be a float

    if inc == None:
        inc = 1.0
    count = int(math.ceil((end - start) / inc))

    L = [None,] * (count + 1)

    L[0] = start
    for i in xrange(1,count + 1):
        L[i] = L[i-1] + inc
    return L
