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

# TODO: 1. one more BUG
#       2. Why involve right cond while dealing with positive direction

import numpy as np

class Spatial_Point():
    def __init__(self, place, material = None, isRequired = True):
        self.material = None
        self.x = place
        self.rek = isRequired

    def isInterface(self):
        return self.material is None

    def material(self):
        return self.material

    def x(self):
        return self.x

    def required():
        return self.rek

class Model_1D_Numerical_Solver():
    cache = {}

    def __init__(self, grid, gauss_discrete_direction_num = 2, total_point_num = -1, discretization_stepsize = -1,):
        """ NOTE: please clean up the discrete_direction_num before call this
            super constructor. Eg, for Gauss Legendre Quadrature, plz ensure it
            is even number.
            TODO: Add slab geometry
        """
        self.n = total_point_num
        self.step_size = discretization_stepsize
        self.discrete_direction_num = discrete_direction_num / 2 * 2
        self.grid = grid
        self.mesh = []
        # TODO: Calculate GAUSS-LEGENDRE QUADRATURE
        if self.discrete_direction_num in cache:
            self.guass_legendre = cache[discrete_direction_num]
        else:
            self.guass_legendre = np.polynomial.legendre.leggauss(self.discrete_direction_num)
            cache[discrete_direction_num] = self.guass_legendre

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
        assert((base_point_num == -1 or base_step_size == -1) and not (base_point_num == -1 and base_step_size == -1), "Please either specify base point number or base step size.")
        Model_1D_Numerical_Solver.__init__(self, grid = grid_model, gauss_discrete_direction_num = (discrete_direction_num / 2) * 2)
        if base_point_num != -1
            self.base_step_size = grid_model.len() / float(base_point_num)
            self.base_point_num = base_point_num
        else:
            self.base_step_size = base_step_size
            self.base_point_num = round(grid_model.len() / base_step_size)
        # Start constructing the mesh, note we will add additional points on
        # interface and inside intercal if the basic points number and step
        # size is not enought to cover all the intervals
        # TODO: BUG may exist, need to debug， Is n * step_size guaranteed to
        #       be the same as total length??
        self.mesh += [Spatial_Point(0.0, None)]
        last_required_point = 0.0
        for interval in grid_model:
            while self.mesh[-1].x() < interval.right():
                next_base_point = self.base_step_size + last_required_point
                if next_base_point < interval.right():
                    # If adding a point does not cause exceeding the interface
                    self.mesh += [Spatial_Point(next_base_point, interval.material())]
                    last_required_point = next_base_point
                elif next_base_point >= interval.right():
                    if interval.left() == self.mesh[-1].x():
                        self.mesh += [Spatial_Point(interval.mid_point(), interval.material()), Spatial_Point(interval.right())]
                    else：
                        self.mesh += [Spatial_Point(interval.right())]

        self.mesh_interval_len = [(self.mesh[i+1].x() - self.mesh[i].x()) for i in range(len(self.mesh) - 1)]
        self.n = len(self.mesh)

    def solve(self):
        """ Build matrix and solve, more docs later
            Currently using three points to simulate. We shall try more points
            later
            # TODO: More flexibility that can use more points (say, 5?)
            #       And add abstraction to simplify the calculation
        """
        mesh_point_num = len(self.mesh) - 1
        matrix_size = discrete_direction_num * mesh_point_numb
        A = [[0.0 for _ in range(matrix_size)] for __ in range(matrix_size)]
        B = [0.0 for _ in range(matrix_size)]
        # Begin constructing the matrix, start by iterating through directions
        u = self.gauss_u()
        wt = self.gauss_weight()
        h = self.mesh_interval_len

        # First half
        for dir_index in range(self.discrete_direction_num / 2):
            dir_submatrix_index = dir_index * mesh_point_num
            # Deal with edge case, in which left part does not exist
            # TODO: Test different points, 2,3,4,5, make this more flexible
            A[dir_submatrix_index][dir_submatrix_index] = -u[dir_index] * (1/h[0]+1/(h[0]+h[1])) + self.mesh[0].material().cross_section() - self.mesh[0].material().scattering_section() * wt[dir_index] / 2.0
            A[dir_submatrix_index][dir_submatrix_index+1] = u[dir_index] * (1/h[0]+1/h[1])
            B[dir_submatrix_index] = self.mesh[0].material().source() / 2.0 - u[dir_index] * h[0]/h[1] * 1.0 / (h[0]+h[1]) * self.grid.left_boundary_condition()
            
            A[dir_submatrix_index + self.n - 1][dir_submatrix_index + self.n - 2] = u[dir_index] / h[-1]
            A[dir_submatrix_index + self.n - 1][dir_submatrix_index + self.n - 1]= -u[dir_index] / h[-1] + self.mesh[-1].material().cross_section() - self.mesh[-1].material().scattering_section() * wt[dir_index] / 2.0
            B[dir_submatrix_index + self.n - 1] = self.mesh[-1].material().source() / 2.0
            # The first and last points are already taken care of
            for spatial_point_index in range(1, len(self.mesh[:-1])):
                cur_index = dir_submatrix_index + spatial_point.index()
                spatial_point = self.mesh[spatial_point_index]
                h_ = h[spatial_point_index]
                _h = h[spatial_point_index - 1]
                if spatial_point.isInterface():
                    next_mat = self.mesh[spatial_point_index + 1].material()
                    h__ = h[spatial_point_index + 1]
                    th = h_ + h__
                    hh = 1/h_ + 1/h__
                    dh = h_/h__
                    A[cur_index][cur_index] = -u[dir_index] * (1.0/h_ + 1.0/th) + next_mat.cross_section() - next_mat.scattering_section() * wt[dir_index] / 2.0
                    A[cur_index][cur_index + 1] = u[dir_index] * hh
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
            # TODO: Better variable naming
            # BUG HERE!!! How to deal with the last point?
            for j in range(1, self.n):
                for i in range(self.discrete_direction_num / 2):
                    if i == dir_index:
                        continue
                    S = i * (self.n - 1)
                    if self.mesh[j].isInterface():
                        A[dir_submatrix_index + j][S + j] = -self.mesh[j+1].scattering_section() * wt[i] / 2.0
                    else:
                        A[dir_submatrix_index + j][S + j] = -self.mesh[j].scattering_section() * wt[i] / 2.0
                        
            for i in range(self.discrete_direction_num / 2):
                S = (self.discrete_direction_num / 2 + i) * (self.n - 1)
                for j in range(1, self.n - 1):
                    if self.mesh[j].isInterface():
                        A[dir_submatrix_index + j][S + j - 1] = -self.mesh[j+1].scattering_section() * wt[self.discrete_direction_num / 2 + i] / 2.0
                    else:
                        A[dir_submatrix_index + j][S + j - 1] = -self.mesh[j].scattering_section() * wt[self.discrete_direction_num / 2 + i] / 2.0
                B[dir_submatrix_index] += self.mesh[0].material().scattering_section() * wt[self.discrete_direction_num / 2 + i] * self.grid.left_boundary_condition() / 2.0

        # Second half, which is basically the same, thus plz refactorizing this part
        # TODO: Recheck formula with Richard, thus leaving the second part unchanged
        for dir_index in range(self.discrete_direction_num / 2, self.discrete_direction_num):
            dir_submatrix_index = dir_index * mesh_point_num
            # Deal with edge case, in which left part does not exist
            # TODO: Test different points, 2,3,4,5, make this more flexible
            A[dir_submatrix_index][dir_submatrix_index] = -u[dir_index] / h[0] + self.mesh[0].material().cross_section() - self.mesh[0].material().scattering_section() * wt[dir_index] / 2.0
            A[dir_submatrix_index][dir_submatrix_index+1] = u[dir_index] / h[0]
            B[dir_submatrix_index] = self.mesh[0].material().source() / 2.0
            
            A[dir_submatrix_index + self.n - 1][dir_submatrix_index + self.n - 2] = u[dir_index] * (1/h[-1]+1/h[-2])
            A[dir_submatrix_index + self.n - 1][dir_submatrix_index + self.n - 1]= -u[dir_index] * (1/h[-1]+1/(h[-1]+h[-2])) + self.mesh[-1].material().cross_section() - self.mesh[-1].material().scattering_section() * wt[dir_index] / 2.0
            B[dir_submatrix_index + self.n - 1] = self.mesh[-1].material().source() / 2.0 - u[dir_index] * h[-1]/h[-2] * 1.0 / (h[-1]+h[-2]) * self.grid.right_boundary_condition()
            # The first and last points are already taken care of
            for spatial_point_index in range(1, len(self.mesh[:-1])):
                cur_index = dir_submatrix_index + spatial_point.index()
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
                    A[cur_index][cur_index + 2] = u[dir_index] * dh * 1.0/th
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
            # TODO: Better variable naming
            # BUG HERE!!! How to deal with the last point?
            for j in range(self.n - 2):
                for i in range(self.discrete_direction_num / 2, self.discrete_direction_num):
                    if i == dir_index:
                        continue
                    S = i * (self.n - 1)
                    if self.mesh[j].isInterface():
                        A[dir_submatrix_index + j][S + j] = -self.mesh[j].scattering_section() * wt[i] / 2.0
                    else:
                        A[dir_submatrix_index + j][S + j] = -self.mesh[j+1].scattering_section() * wt[i] / 2.0
                        
            for i in range(self.discrete_direction_num / 2):
                S = (self.discrete_direction_num / 2 + i) * (self.n - 1)
                for j in range(self.n - 1):
                    if self.mesh[j + 1].isInterface():
                        A[dir_submatrix_index + j][S + j] = -self.mesh[j].scattering_section() * wt[self.discrete_direction_num / 2 + i] / 2.0
                    else:
                        A[dir_submatrix_index + j][S + j]= -self.mesh[j+1].scattering_section() * wt[self.discrete_direction_num / 2 + i] / 2.0
                B[dir_submatrix_index] += self.mesh[0].material().scattering_section() * wt[self.discrete_direction_num / 2 + i] * self.grid.right_boundary_condition() / 2.0
        # Return the solution of this linear system, note the result is both undetermined and unchecked, need to put more tests for this solve procedure and refactorizing this since
        # it is such a big block lol
        return numpy.linalg.solve(A, B)

    def solve_required_points(self):
        complete_solution = self.solve()
        # Then we should throw away all additional points
        


    def add_point(self, x):
        for i in range(len(self.mesh)):
            if self.mesh[i] == x:
                return
            if self.mesh[i] < x and x < self.mesh[i+1]:
                self.mesh.insert(i+1, x)
                return
