# TODO: 1. Documents codes
#       2. Add parameter docs
#       3. Unit testing
#       4. Class docs

# Codes for Nuclear Engineering Department @ Berkeley
# This is a summer project, translation and packaging past codes for 1-D
# Nutron transport simulations.
# Note this file is highly extensible so I suggest it could be used in future
# similar simulation designs, plz improve anything if you'd like, it would be
# my honor.
# Mingjian Lu, July 2016

# TODO: 1. Debug

import numpy as np

class Spatial_Point():
    def __init__(self, place, material = None, isRequired = True):
        self.mmaterial = None
        self.xx = place
        self.rek = isRequired

    def isInterface(self):
        return self.material is None

    def material(self):
        return self.mmaterial

    def x(self):
        return self.xx

    def required():
        return self.rek

class Solution_Point():
    def __init__(self, corresponding_point, index):
        self.cp = corresponding_point
        self.index = index
        self.dir_val = {}

    def add_dir_val(self, dir, val):
        self.dir_val[dir] = val

    def spatial_point(self):
        return self.cp

class Model_1D_Numerical_Solver():
    cache = {}

    def __init__(self, grid, gauss_discrete_direction_num = 2, total_point_num = -1, discretization_stepsize = -1):
        """ NOTE: please clean up the discrete_direction_num before call this
            super constructor. Eg, for Gauss Legendre Quadrature, plz ensure it
            is even number.
            TODO: Add slab geometry
        """
        self.n = total_point_num
        self.step_size = discretization_stepsize
        self.discrete_direction_num = gauss_discrete_direction_num / 2 * 2
        self.grid = grid
        self.mesh = []
        # Note: local cache is for solved solutions in case redundant calculation
        self.local_cache = {}
        # Calculate GAUSS-LEGENDRE QUADRATURE
        if self.discrete_direction_num in Model_1D_Numerical_Solver.cache:
            self.guass_legendre = Model_1D_Numerical_Solver.cache[discrete_direction_num]
        else:
            self.guass_legendre = np.polynomial.legendre.leggauss(self.discrete_direction_num)
            Model_1D_Numerical_Solver.cache[self.discrete_direction_num] = self.guass_legendre

    def point_num(self):
        self.__check_n()
        return self.n

    def step_size(self):
        self.__check_s()
        return self.step_size

    def avg_step_size(self):
        self.__check_n()
        return float(self.grid.len()) / float(self.n)

    def gauss_u(self):
        return self.guass_legendre[0]

    def gauss_weight(self):
        return self.guass_legendre[1]

    def __check_n(self):
        if self.n == -1:
            raise Exception("The grid needs to be generated since you are using a nontrivial grid.")

    def __check_s(self):
        if self.n == -1:
            raise Exception("This is a complex grid so there is no fixed step size, plz use avg_step_size method instead.")

    def check_init_state(self):
        return len(self.mesh) > 1

    def add_point(self, x):
        raise NotImplementedError

    def solve(self):
        raise NotImplementedError

class Model_1D_Stochastic_Finite_Step_Solver(Model_1D_Numerical_Solver):
    """ This is the finite step method solver for the stochastic model.
        Note theortically finite volumn method should perform better than this.
        1. The mesh generated by this class needs revision: TODO: 10/7/2016
        2. abstract general method to super class
    """
    def __init__(self, grid_model, base_point_num = -1, base_step_size = -1, gauss_discrete_direction_num = 2):
        assert (base_point_num != -1 or base_step_size != -1) and not (base_point_num != -1 and base_step_size != -1), "Please either specify base point number or base step size."
        Model_1D_Numerical_Solver.__init__(self, grid = grid_model, gauss_discrete_direction_num = (gauss_discrete_direction_num / 2) * 2)
        if base_point_num != -1:
            self.base_step_size = grid_model.len() / float(base_point_num)
            self.base_point_num = base_point_num
        else:
            self.base_step_size = base_step_size
            self.base_point_num = round(grid_model.len() / base_step_size)
        # Start constructing the mesh, note we will add additional points on
        # interface and inside intercal if the basic points number and step
        # size is not enought to cover all the intervals
        # TODO: BUG may exist, need to debug, Is n * step_size guaranteed to
        #       be the same as total length?
        self.mesh += [Spatial_Point(0.0, None)]
        last_required_point = 0.0
        for interval in grid_model:
            while self.mesh[-1].x() < interval.right():
                next_base_point = min(self.base_step_size + last_required_point, grid_model.len())
                if next_base_point < interval.right():
                    # If adding a point does not cause exceeding the interface
                    self.mesh += [Spatial_Point(next_base_point, interval.material())]
                    last_required_point = next_base_point
                elif next_base_point >= interval.right():
                    if interval.left() == self.mesh[-1].x():
                        self.mesh += [Spatial_Point(interval.mid_point(), interval.material(), False), Spatial_Point(interval.right(), isRequired = next_base_point == interval.right())]
                    else:
                        self.mesh += [Spatial_Point(interval.right(), isRequired = next_base_point == interval.right())]
                    if next_base_point == interval.right():
                        last_required_point = next_base_point
                # print(str(self.mesh[-1].x()) + "/" + str(interval.right()))

        self.mesh_interval_len = [(self.mesh[i+1].x() - self.mesh[i].x()) for i in range(len(self.mesh) - 1)]
        self.n = len(self.mesh)

    def solve(self):
        """ Build matrix and solve, more docs later
            Currently using three points to simulate. We shall try more points
            later
            # TODO: Add abstraction to simplify the calculation
            #       Docs and param explainations
        """
        mesh_point_num = len(self.mesh) - 1
        matrix_size = discrete_direction_num * mesh_point_numb
        A = [[0.0 for _ in range(matrix_size)] for __ in range(matrix_size)]
        B = [0.0 for _ in range(matrix_size)]
        # Begin constructing the matrix, start by iterating through directions
        u = self.gauss_u()
        wt = self.gauss_weight()
        h = self.mesh_interval_len

        # Upper half
        for dir_index in range(self.discrete_direction_num / 2):
            dir_submatrix_index = dir_index * mesh_point_num
            # Deal with edge case, in which left part does not exist
            # TODO: Test different points, 2,3,4,5, make this more flexible
            A[dir_submatrix_index][dir_submatrix_index] = -u[dir_index] / h[0] + self.mesh[0].material().cross_section() - self.mesh[0].material().scattering_section() * wt[dir_index] / 2.0
            A[dir_submatrix_index][dir_submatrix_index+1] = u[dir_index] / h[0]
            B[dir_submatrix_index] = self.mesh[0].material().source() / 2.0

            A[dir_submatrix_index + self.n - 1][dir_submatrix_index + self.n - 2] = -u[dir_index] * h[-1]/h[-2] * (h[-1] + h[-2])
            A[dir_submatrix_index + self.n - 1][dir_submatrix_index + self.n - 1]= u[dir_index] * (h[-1] - h[-2])/(h[-1]*h[-2]) + self.mesh[-1].material().cross_section() - self.mesh[-1].material().scattering_section() * wt[dir_index] / 2.0
            B[dir_submatrix_index + self.n - 1] = self.mesh[-1].material().source() / 2.0 - u[dir_index] * h[-2]/h[-1] * 1.0 / (h[-1]+h[-2]) * self.grid.right_boundary_condition()
            # The first and last points are already taken care of
            for spatial_point_index in range(1, len(self.mesh[:-1])):
                cur_index = dir_submatrix_index + spatial_point_index
                spatial_point = self.mesh[spatial_point_index]
                h_ = h[spatial_point_index]
                _h = h[spatial_point_index - 1]
                if spatial_point.isInterface():
                    next_mat = self.mesh[spatial_point_index + 1].material()
                    h__ = h[spatial_point_index + 1]
                    th = h_ + h__
                    hh = 1.0/h_ + 1.0/h__
                    dh = h_/h__
                    A[cur_index][cur_index] = -u[dir_index] * (1.0/h_ + 1.0/th) + next_mat.cross_section() - next_mat.scattering_section() * wt[dir_index] / 2.0
                    A[cur_index][cur_index + 1] = u[dir_index] * hh
                    if spatial_point_index == mesh_point_num - 1:
                        # If the point is the last third point, then it does not have enough points to do the calculation, we need to deal with this using boundary condition
                        B[cur_index] = next_mat.source() / 2.0 + u[dir_index] * dh * 1.0/th * self.grid.right_boundary_condition()
                    else:
                        A[cur_index][cur_index + 2] = -u[dir_index] * dh * 1.0/th
                        B[cur_index] = next_mat.source() / 2.0
                else:
                    cur_mat = spatial_point.material()
                    th = h_ + _h
                    dh = h_/_h
                    ddh = h_ - _h
                    mh = h_ * _h
                    A[cur_index][cur_index - 1] = -u[dir_index] * dh * 1.0/th
                    A[cur_index][cur_index] = u[dir_index] * ddh/mh + cur_mat.cross_section() - cur_mat.scattering_section() * wt[dir_index] / 2.0
                    A[cur_index][cur_index + 1] = u[dir_index] * 1.0/dh * 1.0/th
                    B[cur_index] = cur_mat.source() / 2.0
            # TODO: Var_naming confirm with Richard
            for tmp_pt_index in range(mesh_point_num):
                for tmp_dir_index in range(self.discrete_direction_num / 2):
                    if tmp_dir_index == dir_index:
                        continue
                    negative_dir_in_matrix_index = tmp_dir_index * mesh_point_num
                    matrix_entry_index_row = dir_submatrix_index + tmp_pt_index
                    matrix_entry_index_col = dir_submatrix_index + negative_dir_in_matrix_index
                    if self.mesh[tmp_pt_index].isInterface():
                        A[matrix_entry_index_row][matrix_entry_index_col] = -self.mesh[tmp_pt_index+1].scattering_section() * wt[tmp_dir_index] / 2.0
                    else:
                        A[matrix_entry_index_row][matrix_entry_index_col] = -self.mesh[tmp_pt_index].scattering_section() * wt[tmp_dir_index] / 2.0

            for tmp_dir_index in range(self.discrete_direction_num / 2):
                positive_dir_in_matrix_index = (self.discrete_direction_num / 2 + tmp_dir_index) * mesh_point_num
                for tmp_pt_index in range(1, mesh_point_num):
                    matrix_entry_index_row = dir_submatrix_index + tmp_pt_index
                    matrix_entry_index_col = positive_dir_in_matrix_index + tmp_pt_index - 1
                    if self.mesh[tmp_pt_index].isInterface():
                        A[matrix_entry_index_row][matrix_entry_index_col] = -self.mesh[tmp_pt_index+1].scattering_section() * wt[self.discrete_direction_num / 2 + tmp_dir_index] / 2.0
                    else:
                        A[matrix_entry_index_row][matrix_entry_index_col] = -self.mesh[tmp_pt_index].scattering_section() * wt[self.discrete_direction_num / 2 + tmp_dir_index] / 2.0
                B[dir_submatrix_index] += self.mesh[0].material().scattering_section() * wt[self.discrete_direction_num / 2 + tmp_dir_index] * self.grid.right_boundary_condition() / 2.0

        # Second half, which is basically the same, thus plz refactorizing this part
        for dir_index in range(self.discrete_direction_num / 2, self.discrete_direction_num):
            dir_submatrix_index = dir_index * mesh_point_num
            # Deal with edge case, in which left part does not exist
            # TODO: Test different points, 2,3,4,5, make this more flexible
            A[dir_submatrix_index][dir_submatrix_index] = u[dir_index] * (h[1] - h[0])/(h[0] * h[0]) + self.mesh[0].material().cross_section() - self.mesh[0].material().scattering_section() * wt[dir_index] / 2.0
            A[dir_submatrix_index][dir_submatrix_index+1] = u[dir_index] * h[0]/h[1] * (h[0]+h[1])
            B[dir_submatrix_index] = self.mesh[0].material().source() / 2.0 + u[dir_index] * h[1]/h[0] * 1.0 / (h[0]+h[1]) * self.grid.left_boundary_condition()

            A[dir_submatrix_index + self.n - 1][dir_submatrix_index + self.n - 2] = u[dir_index] / h[-1]
            A[dir_submatrix_index + self.n - 1][dir_submatrix_index + self.n - 1]= -u[dir_index] / h[-1] + self.mesh[-1].material().cross_section() - self.mesh[-1].material().scattering_section() * wt[dir_index] / 2.0
            B[dir_submatrix_index + self.n - 1] = self.mesh[-1].material().source() / 2.0
            # The first and last points are already taken care of
            for spatial_point_index in range(1, len(self.mesh[:-1])):
                cur_index = dir_submatrix_index + spatial_point_index
                spatial_point = self.mesh[spatial_point_index]
                h_ = h[spatial_point_index]
                _h = h[spatial_point_index - 1]
                if spatial_point.isInterface():
                    prev_mat = self.mesh[spatial_point_index - 1].material()
                    __h = h[spatial_point_index - 2]
                    th = _h + __h
                    hh = 1/_h + 1/__h
                    dh = _h/__h
                    A[cur_index][cur_index] = u[dir_index] * (1.0/_h + 1.0/th) + prev_mat.cross_section() - prev_mat.scattering_section() * wt[dir_index] / 2.0
                    A[cur_index][cur_index + 1] = -u[dir_index] * hh
                    if spatial_point_index == 2:
                        # If the point is the first third point, then it does not have enough points to do the calculation, we need to deal with this using boundary condition
                        B[cur_index] = prev_mat.source() / 2.0 - u[dir_index] * dh * 1.0/th * self.grid.left_boundary_condition()
                    else:
                        A[cur_index][cur_index + 2] = u[dir_index] * dh * 1.0/th
                        B[cur_index] = prev_mat.source() / 2.0
                else:
                    cur_mat = spatial_point.material()
                    th = h_ + _h
                    dh = h_/_h
                    ddh = h_ - _h
                    mh = h_ * _h
                    A[cur_index][cur_index - 1] = -u[dir_index] * dh * 1.0/th
                    A[cur_index][cur_index] = u[dir_index] * ddh/mh + cur_mat.cross_section() - cur_mat.scattering_section() * wt[dir_index] / 2.0
                    A[cur_index][cur_index + 1] = u[dir_index] * 1.0/dh * 1.0/th
                    B[cur_index] = cur_mat.source() / 2.0
            # TODO: Var_naming confirm with Richard
            for tmp_pt_index in range(mesh_point_num):
                for tmp_dir_index in range(self.discrete_direction_num / 2, self.discrete_direction_num):
                    if tmp_dir_index == dir_index:
                        continue
                    positive_dir_in_matrix_index = tmp_dir_index * mesh_point_num
                    matrix_entry_index_row = dir_submatrix_index + tmp_pt_index
                    matrix_entry_index_col = dir_submatrix_index + positive_dir_in_matrix_index
                    if self.mesh[tmp_pt_index + 1].isInterface():
                        A[matrix_entry_index_row][matrix_entry_index_col] = -self.mesh[tmp_pt_index].scattering_section() * wt[tmp_dir_index] / 2.0
                    else:
                        A[matrix_entry_index_row][matrix_entry_index_col] = -self.mesh[tmp_pt_index + 1].scattering_section() * wt[tmp_dir_index] / 2.0

            for tmp_dir_index in range(self.discrete_direction_num / 2):
                negative_dir_in_matrix_index = tmp_dir_index * mesh_point_num
                for tmp_pt_index in range(mesh_point_num - 1):
                    matrix_entry_index_row = dir_submatrix_index + tmp_pt_index
                    matrix_entry_index_col = negative_dir_in_matrix_index + tmp_pt_index + 1
                    if self.mesh[tmp_pt_index + 1].isInterface():
                        A[matrix_entry_index_row][matrix_entry_index_col] = -self.mesh[tmp_pt_index].scattering_section() * wt[self.discrete_direction_num / 2 + tmp_dir_index] / 2.0
                    else:
                        A[matrix_entry_index_row][matrix_entry_index_col] = -self.mesh[tmp_pt_index + 1].scattering_section() * wt[self.discrete_direction_num / 2 + tmp_dir_index] / 2.0
                B[dir_submatrix_index] += self.mesh[0].material().scattering_section() * wt[self.discrete_direction_num / 2 + tmp_dir_index] * self.grid.right_boundary_condition() / 2.0
        # Return the solution of this linear system, note the result is both undetermined and unchecked, need to put more tests for this solve procedure and refactorizing this since
        # it is such a big block lol
        self.local_cache['complete_solution'] = numpy.linalg.solve(A, B)
        return self.local_cache['complete_solution']

    def solve_required_points(self):
        if 'complete_solution' in self.local_cache:
            complete_solution = self.local_cache['complete_solution']
        else:
            complete_solution = self.solve()
        # Then we should throw away all additional points
        # Plus we need to take boundary condition into account
        # Note Im using solution points class to make the complete required solution easily accessible
        req_solution = [Solution_Point(self.mesh[p_index], p_index) for p_index in range(len(self.mesh)) if self.mesh[p_index].required()]
        req_set = {}
        for p_ind in range(len(req_solution)):
            p = req_solution[p_ind]
            req_set[p.index] = p_ind
        for cur_dir in range(self.discrete_direction_num):
            for cur_p in range(self.n - 1):
                if cur < self.discrete_direction_num / 2:
                    ax_p = cur_p
                    req_solution[-1].add_dir_val(cur_dir, self.grid.right_boundary_condition())
                else:
                    ax_p = cur_p + 1
                    req_solution[0].add_dir_val(cur_dir, self.grid.left_boundary_condition())
                if ax_p in req_set:
                    req_solution[req_set[ax_p]].add_dir_val(cur_dir, complete_solution[cur_dir * (self.n - 1) + cur_p])
        self.local_cache['req_solution'] = req_solution
        return req_solution

    def solve_scalar_flux(self):
        if 'req_solution' in self.local_cache:
            req_solution = self.local_cache['req_solution']
        else:
            req_solution = self.solve_required_points()
        scalar_flux = [0.0] * len(req_solution)
        wt = self.gauss_weight()
        for solution_p_index in range(len(req_solution)):
            solution_p = req_solution[solution_p]
            scalar_flux[solution_p_index] = sum([wt[dir_index] * solution_p.dir_val[dir_index] for dir_index in range(self.discrete_direction_num)])
        return scalar_flux

    # def plot_scalar_flux(self):
    #     if 'scalar_flux' in self.local_cache:
    #         scalar_flux = self.local_cache['scalar_flux']
    #     else:
    #         scalar_flux = self.solve_scalar_flux()
    #     req_solution = self.local_cache['req_solution']
    #     plt.plot([p.spatial_point().x() for p in req_solution], scalar_flux)
    #     plt.show()

    def add_point(self, x):
        for i in range(len(self.mesh)):
            if self.mesh[i] == x:
                return
            if self.mesh[i] < x and x < self.mesh[i+1]:
                self.mesh.insert(i+1, x)
                return

class Model_1D_Periodic_Solver(Model_1D_Numerical_Solver):
        def __init__(self, grid_model, point_num, gauss_discrete_direction_num = 2, start_index = 0):
            for mat in grid_model.material_list():
                assert int(mat.thickness() / self.step_size) == mat.thickness() / self.step_size, "Interfaces must be included in discretization."
            Model_1D_Numerical_Solver.__init__(self, grid = grid_model, gauss_discrete_direction_num = (discrete_direction_num / 2) * 2, total_point_num = point_num, discretization_stepsize = total_len / point_num)
            self.start_index = start_index

        def solve(self):
            mat1, mat2 = self.grid.materials[0], self.grid.materials[1]
            T = self.grid.len
            m1,m2 = mat1.thickness(), mat2.thickness()
            n,N = self.n, self.discrete_direction_num
            Es1,Es2 = mat1.scattering_section(), mat2.scattering_section()
            Et1,Et2 = mat1.cross_section(), mat2.cross_section()
            yo,y_ = self.grid.left_boundary_condition(), self.grid.right_boundary_condition()
            Q1,Q2 = mat1.source(), mat2.source()
            u,wt = self.gauss_u(), self.gauss_weight()
            a = self.start_index
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

class Model_1D_Homogeneous_Solver(Model_1D_Numerical_Solver):
        def __init__(self, grid_model, point_num, gauss_discrete_direction_num = 2, start_index = 0):
            for mat in grid_model.material_list():
                assert int(mat.thickness() / self.step_size) == mat.thickness() / self.step_size, "Interfaces must be included in discretization."
            Model_1D_Numerical_Solver.__init__(self, grid = grid_model, gauss_discrete_direction_num = (discrete_direction_num / 2) * 2, total_point_num = point_num, discretization_stepsize = total_len / point_num)

        def solve(self):
            return
