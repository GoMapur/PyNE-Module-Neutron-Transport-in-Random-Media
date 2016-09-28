from OneDGridModelBenchMark import *
from OneDGridModel import *
from OneDGridModelSolver import *

materials1 = Model_Material(thickness = 1.0, cross_section = 1.0, scattering_ratio = 1.0, homogeneous_isotropic_source = 0.0)
materials2 = Model_Material(thickness = 1.0, cross_section = 1.0, scattering_ratio = 1.0, homogeneous_isotropic_source = 0.0)
bench = Model_1D_Periodic_Finite_Step_Benchmark(total_len = 10.0, point_num = 100, boundary_cond = [1.0, 0.0], materials = [materials1, materials2], gauss_discrete_direction_num = 2)
bench.benchmark()
