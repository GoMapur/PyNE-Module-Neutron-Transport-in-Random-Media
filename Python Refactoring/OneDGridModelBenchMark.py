# Codes for Nuclear Engineering Department @ Berkeley
# This is a summer project, translation and packaging past codes for 1-D
# Nutron transport simulations.
# Note this file is highly extensible so I suggest it could be used in future
# similar simulation designs, plz improve anything if you'd like, it would be
# my honor.
# Mingjian Lu, July 2016

import numpy as np
from OneDGridModel import *
from OneDGridModelSolver import *

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
        return solution

    def benchmark(self):
        # TODO: This fixed number should changed to be a indicator
        #       showing when should we stop
        # while i < 100000:
        #     solution = self.benchmark_once()
        print(self.benchmark_once())

class Model_1D_Periodic_Finite_Step_Benchmark(Model_1D_Benchmark):
    def __init__(self, total_len, point_num, boundary_cond, materials, gauss_discrete_direction_num = 2):
        Model_1D_Benchmark.__init__(self, total_len, point_num, boundary_cond, materials, gauss_discrete_direction_num)

    def benchmark(self):
        grid_model = Periodic_Grid(total_len, boundary_cond, materials)
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
            print('Problem ' + str(a) + ' of ' + str(total))
            solver = Model_1D_Periodic_Solver(grid_model, self.point_num, self.gauss_discrete_direction_num, a)
            Z,n1,B,L,A, extra = solver.solve(T,m1,m2,n,N,Es1,Es2,Et1,Et2,yo,y_,Q1,Q2,u,wt,a)
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

class Model_1D_Homogeneous_Finite_Step_Benchmark(Model_1D_Benchmark):
    def __init__(self, total_len, point_num, boundary_cond, materials, gauss_discrete_direction_num = 2):
        Model_1D_Benchmark.__init__(self, total_len, gauss_discrete_direction_num, point_num, boundary_cond, materials)

    def benchmark(self):
        return
