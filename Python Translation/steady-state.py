# %THIS CODE SOLVES THE STEADY-STATE, MONOENERGETIC TRANSPORT EQUATION
# %IN A HOMOGENEOUS MEDIUM WITH ISOTROPIC SCATTERING AND ISOTROPIC SOURCE
# %IN ROD AND SLAB (SN) GEOMETRIES. IT USES 2-POINT CENTRAL DIFFERENCING
# %(ORDER H^2) WITH 3-POINT FORWARD/BACKWARD AT THE BOUNDARIES.

# %MAIN OUTPUTS: transmission 'TR'; reflection 'RL'
# %			  scalar flux 'SCAL'

import numpy as np
import matplotlib.pyplot as plt
import os, sys
import random

def config_file(): #only supports homogenous right now. When refactoring code to clean out repeat functions, make generic
	folder = os.path.join(os.getcwd(), os.path.dirname(__file__))
	config = os.path.join(folder, 'steady-state.config')
	result = []
	#add error checking later
	with open(config, 'rb') as f:
		for line in f:
			line = line.strip()
			if len(line) == 0:
				continue
			if line[0] == '%':
				continue
			line = line[1:len(line)-1].split(',')
			if len(line) == 9:
				if line[0]:
					if line[0] == 'SLAB':
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
				print 'Malformed input data, please check.'
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
			plt.savefig(name +'.png')
		else:
			plt.show()
		plt.clf()
	except Exception as e:
		print e
		print len(x), SCAL.shape
	return
	# not compatible with virtualenv

def save_csv(matrix, name):
	np.savetxt('csv/' + name + '.csv', matrix, delimiter=',')

#extra order does matter
#check h ascending order

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



def main(sysargs):
	# default state is file configuration

	# if not len(sysargs):
	# 	save = False
	# elif sysargs[0].lower() == 'manual':
	# 	var = homogenous_config_manual()
	# elif sysargs[0].lower() == 'save':
	# 	save = True
	# var = config_file()
	# for v in var:
	# 	homogenous(v, save)

	#lots of redudant code, need to make them generic later

	periodic([0, 2, 20.0, 100, 0.0, 0.0, 1.0, 1.0, 1.0, 0.5, 1.0, 0.0, 0.0, 0.0])

if __name__ == '__main__':
	main(sys.argv[1:])
