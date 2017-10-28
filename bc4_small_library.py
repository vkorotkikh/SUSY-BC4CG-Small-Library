# ******************************************************************************
# Name:    Verifying BC4 Space Vij coefficient
# Author:  Vadim Korotkikh
# Email:   va.korotki@gmail.com
# Date:    November 2016
# Version: 1.0A
#
# Description: Function for verifying the BC4 space coefficient library
#
# ******************************************************************************
# Library Imports

import re
import sys
import math
import time
import numpy.matlib
import itertools
import numpy as np
from numpy import array
from numpy.linalg import inv

# ******************************************************************************
# Function Imports
# import matrix_calc_vijmat2
# import matrix_calc_vijmat_nopr
import vij_holoraumy_calc
import matrix_calc_vijmat2

p_switch	= 0
# ******************************************************************************
# Main() function.
def main(pset_str):

	print("# ***********************************************************************")
	print("# Name:    Calculating BC4 CG Small Library")
	print("# Author:  Vadim Korotkikh	")
	print("# Email:   va.korotki@gmail.com")
	print("# Date:    June 2017	")
	print("# Version: N/A")
	print("#	")
	# print("# Description: Function for verifying the BC4 space coefficient library")
	# print("#	")
	print("# ***********************************************************************")
	print("		")

	""" PALL seems to be broken for now due to lack of if elif logic at end
		of string_to_tetrad() function """
	""" Options for BC4 Library:
		PALL - Entire small library, one at a time.
		P1, P2, P3, P4, P5, P6 - Only one section. """

	# pset_str = "PALL"
	# pset_str = "P1"
	bc4cg_holoraumy_calc(pset_str)

# ******************************************************************************
# BC4 Calculation options organizer -  director
def bc4_validation_organizer(pset_arg, *args):
	''' Based on pset_arg and *args figures out which function to pass
		arguments too.
		PALL, ALL, P1 - P6,  'coef' ->  only fermi
		PALL, P1 - P6, 'mats' - fermi / boson
	'''

	args_tuple	= args
	argslow	= [x.lower() for x in args_tuple]
	print(argslow, type(args))
	# sys.exit()
	if all(isinstance(argx, str) for argx in args_tuple):
		print("All *args are strings")
		if len(args_tuple) == 2 and args_tuple[0] == 'coef':
			''' redirect to calculating the coefficients '''
			pass
		elif len(args_tuple) == 2 and 'coef' in args_tuple:
			if 'boson' in args_tuple:
				print("ERROR")
				# sys.exit()5
			elif 'fermi' in args_tuple and not 'boson' in args_tuple:
				# pass args to whatever function will do this
				pass

		elif len(args_tuple) == 2 and 'mats' in args_tuple:
			if 'boson' in args_tuple or args_tuple[1] == 'boson':
				bosind = args_tuple.index('boson')
				matind = args_tuple.index('mats')
				bc4cg_holoraumy_mats(pset_arg, matind, bosind)
			elif 'fermi' in args_tuple or args_tuple[1] == 'fermi':
				ferind = args_tuple.index('fermi')
				matind = args_tuple.index('mats')
			else:
				print("this shouldn't have happened")

		elif len(args_tuple) == 2 and 'Vmats' in args_tuple:
			if 'boson' in args_tuple or args_tuple[1] == 'boson':
				bosind = args_tuple.index('boson')
				matind = args_tuple.index('mats')
				bc4cg_holoraumy_mats(pset_arg, matind, bosind)
			elif 'fermi' in args_tuple or args_tuple[1] == 'fermi':
				ferind = args_tuple.index('fermi')
				matind = args_tuple.index('mats')
				bc4cg_holoraumy_mats(pset_arg, matind, bind)
			else:
				print("this shouldn't have happened")

		else:
			# What goes here? Technically this shuldn't happen. Unless I change
			# the code later
			pass
	# bc4_validation_organizer('PALL', 'boson', 'mats')
	# bc4_validation_organizer('PALL', 'fermi', 'coef')


# ******************************************************************************
# BC4 Validation function process organizer
def bc4cg_holoraumy_mats(pset_arg, *args):
	"""
	ITT, pset is tracked via P1, P2, P3....so pset_arg must be a Pn or
	it must be made one. Otherwise function will not executed
	"""

	psets_list		= ["P1", "P2", "P3", "P4", "P5", "P6"]
	psets_dict		= {}
	holotype		= ""

	''' Figure out *args and implement calc accordingly '''
	args_tuple	= args
	arsglen		= len(args)
	if all(isinstance(xarg, str) for xarg in args_tuple):
		print("all strings")
		# if 'mats' in args_tuple or 'Vmats' in args_tuple:
		for xarg in args_tuple:
			if xarg == 'fermi':
				holotype = 'fermionic'
			elif xarg == 'boson':
				holotype = 'bosonic'
			else:
				print("bc4cg_holoraumy_mats-ERROR")  # configure this better later

	if 'fermi' in args_tuple:
		sys.exit("Fermi borked for now")

	""" Generate individual Library slices	"""
	for ps in psets_list:
		# temp_tetrad	= tetrad_setgen(ps)
		temp_tetrad = tetrad_setgen_detailed(ps)
		if ps not in psets_dict:
			psets_dict['%s' % ps] = temp_tetrad
		else:
			print("Unknown P set -  ERROR")
		'''
		adding slice info	'''
		looplist = tetrad_setgen_detailed(ps)
		print(len(looplist), len(looplist[0]), "", type(looplist[0]))


	if pset_arg in psets_list:
	# elif pset_arg != "PALL":
		temp_plist = psets_dict[pset_arg]
		print("		")
		print("Execute Bosonic Holoraumy Calc for", pset_arg)
		print("		")
		vij_holoraumy_calc.calc_holoraumy_mats(temp_plist, pset_arg, holotype)
		print("Holoraumy calc. for:", pset_arg,"finished")
		print("")
	elif pset_arg == "PALL":
		for pset, plist in sorted(psets_dict.items()):
			pint = 0
			pint = int(pset.lstrip("P"))
			print("		")
			print("Execute Bosonic Holoraumy Calc for", pset)
			print("		")
			vij_holoraumy_calc.calc_holoraumy_mats(plist, pset, holotype)
			print("Holoraumy calc. for:", pset,"finished")
			print("")

# ******************************************************************************
def tetrad_setgen(pset):
	"""	Generates a {P#} set of tetrads via execution of pset_string_format()
		and	string_to_tetrad() function. Creates a set of python numpy tetrads
		from string representation	"""

	pset_boold 	= []
	""" Transform P# slices into list index by getting the # """
	pint 	= 0
	if pset[0].lower() == 'p':
		pint = int(pset.lstrip("P")) - 1
	else:
		pint = int(pset) - 1
	# pint 	= int(pset.lstrip("P")) - 1
	pslice, pdef 	= bc4cg_libsets(pint)
	if p_switch:
		print("# ********************************")
		print("Starting conversion process for P slice:", pset)
		# print("Length of", pset, "tetrad set:", len(run_group))
		print("")
	""" Perform boolean calculations """
	pbools	= flips_org_lib(pset)

	for ind, booleans in enumerate(pbools):
		# print("")
		# print(booleans)
		bool_list = []
		for bins in booleans:
			bins_list 	= binaries(bins)
			# temp		= np.array(bins_list)
			bool_mat	= np.asmatrix(np.diag(bins_list))
			bool_list.append(bool_mat)
		# print(bool_list)
		# booled_adinkra 	= [(np.dot(pslice[x], bool_list[x])) for x in range(0, len(pslice))]
		booled_adinkra	= [(np.dot(bool_list[x], pslice[x])) for x in range(0, len(pslice))]
		if p_switch:
			print("Boolean Factor for: ", booleans)
			for i in booled_adinkra:
				print(i)
		pset_boold.append(booled_adinkra)

		# temp_adinkra	= [(np.dot(bool_list[x],pslice[x])) for x in range(0,len(pslice))]
	return pset_boold

# ******************************************************************************
def tetrad_setgen_detailed(pset):
	"""	Generates a {P#} set of tetrads via execution of pset_string_format()
		and	string_to_tetrad() function. Creates a set of python numpy tetrads
		from string representation
		Returns a list containing tuples ( adinkra, (booleans, pslice))
	"""

	pset_boold 	= []
	""" Transform P# slices into list index by getting the #
		subtract 1 for `index`ing
	"""
	pint 	= 0
	pint 	= int(pset.lstrip("P")) - 1
	pslice, pdef 	= bc4cg_libsets(pint)
	if p_switch:
		print("# ********************************")
		print("Starting conversion process for P slice:", pset)
		# print("Length of", pset, "tetrad set:", len(run_group))
		print("")
	""" Perform boolean calculations """
	pbools	= flips_org_lib(pset)

	for ind, booleans in enumerate(pbools):
		bool_list = []
		for bins in booleans:
			bins_list 	= binaries(bins)
			# temp		= np.array(bins_list)
			bool_mat	= np.asmatrix(np.diag(bins_list))
			bool_list.append(bool_mat)
		# booled_adinkra 	= [(np.dot(pslice[x], bool_list[x])) for x in range(0, len(pslice))]
		booled_adinkra	= [(np.dot(bool_list[x], pslice[x])) for x in range(0, len(pslice))]
		if p_switch:
			print("Boolean Factor for: ", booleans)
			for i in booled_adinkra:
				print(i)
		# pset_boold.append(booled_adinkra)
		pset_boold.append((booled_adinkra, (booleans, pdef)))

	return pset_boold

##************************************
# Defining the P slices of original BC4 CG library
def bc4cg_libsets(p_index):
# def lib_pslices(p_index):

	""" {P1} = { (243), (123), (134), (142) }	"""
	p1	= 	[np.matrix([[1, 0, 0, 0], [0, 0, 0, 1], [0, 1, 0, 0], [0, 0, 1, 0]]),
			np.matrix([[0, 1, 0, 0], [0, 0, 1, 0], [1, 0, 0, 0], [0, 0, 0, 1]]),
			np.matrix([[0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 0, 1], [1, 0, 0, 0]]),
			np.matrix([[0, 0, 0, 1], [1, 0, 0, 0], [0, 0, 1, 0], [0, 1, 0, 0]])
			]

	""" {P2} = { (234), (124), (132), (143) }	"""
	p2	=	[np.matrix([[1, 0, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1], [0, 1, 0, 0]]),
			# np.matrix([[0, 0, 0, 1], [1, 0, 0, 0], [0, 0, 1, 0], [0, 1, 0, 0]]),
			np.matrix([[0, 1, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0], [1, 0, 0, 0]]),
			np.matrix([[0, 0, 1, 0], [1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 1]]),
			np.matrix([[0, 0, 0, 1], [0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 1, 0]])
			]
	""" {P3} = { (1243), (23), (14), (1342) } 	"""
	p3	=	[np.matrix([[0, 1, 0, 0], [0, 0, 0, 1], [1, 0, 0, 0], [0, 0, 1, 0]]),
			np.matrix([[1, 0, 0, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 0, 1]]),
			np.matrix([[0, 0, 0, 1], [0, 1, 0, 0], [0, 0, 1, 0], [1, 0, 0, 0]]),
			np.matrix([[0, 0, 1, 0], [1, 0, 0, 0], [0, 0, 0, 1], [0, 1, 0, 0]])
			]
	""" {P4} = { (24), (1234), (13), (1432) } 	"""
	p4	= 	[np.matrix([[1, 0, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0], [0, 1, 0, 0]]),
			 np.matrix([[0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1], [1, 0, 0, 0]]),
			 np.matrix([[0, 0, 1, 0], [0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 0, 1]]),
			 np.matrix([[0, 0, 0, 1], [1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0]])
			]
	""" {P5} = { (34), (12), (1324), (1423) }	"""
	p5	=	[np.matrix([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]]),
			 np.matrix([[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]),
		 	 np.matrix([[0, 0, 1, 0], [0, 0, 0, 1], [0, 1, 0, 0], [1, 0, 0, 0]]),
			 np.matrix([[0, 0, 0, 1], [0, 0, 1, 0], [1, 0, 0, 0], [0, 1, 0, 0]])
			]
	""" {P6} = { (), (12)(34), (13)(24), (14)(23) }	"""
	p6	=	[np.matrix([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]),
			np.matrix([[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]]),
			np.matrix([[0, 0, 1, 0], [0, 0, 0, 1], [1, 0, 0, 0], [0, 1, 0, 0]]),
			np.matrix([[0, 0, 0, 1], [0, 0, 1, 0], [0, 1, 0, 0], [1, 0, 0, 0]])
			]

	p_strings = [ "{(243), (123), (134), (142)}", "{(234), (124), (132), (143)}",
			 	"{(1243), (23), (14), (1342)}", "{(24), (1234), (13), (1432)}",
				"{(34), (12), (1324), (1423)}", "{(), (12)(34), (13)(24), (14)(23)}"]

	p_slices = [ p1, p2, p3, p4, p5, p6 ]

	return p_slices[p_index], p_strings[p_index]

##************************************
# Defining the Pizza slices
""" Definitions used in Adinkra Condense paper """
def pieslices(pie_index):

	p1	= 	[np.matrix([[1, 0, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1], [0, 1, 0, 0]]),
			np.matrix([[0, 1, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0], [1, 0, 0, 0]]),
			np.matrix([[0, 0, 1, 0], [1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 1]]),
			np.matrix([[0, 0, 0, 1], [0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 1, 0]])
			]

	p2	= 	[np.matrix([[1, 0, 0, 0], [0, 0, 0, 1], [0, 1, 0, 0], [0, 0, 1, 0]]),
			np.matrix([[0, 1, 0, 0], [0, 0, 1, 0], [1, 0, 0, 0], [0, 0, 0, 1]]),
			np.matrix([[0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 0, 1], [1, 0, 0, 0]]),
			np.matrix([[0, 0, 0, 1], [1, 0, 0, 0], [0, 0, 1, 0], [0, 1, 0, 0]])
			]

	p3	=	[np.matrix([[1, 0, 0, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 0, 1]]),
			np.matrix([[0, 1, 0, 0], [0, 0, 0, 1], [1, 0, 0, 0], [0, 0, 1, 0]]),
			np.matrix([[0, 0, 1, 0], [1, 0, 0, 0], [0, 0, 0, 1], [0, 1, 0, 0]]),
			np.matrix([[0, 0, 0, 1], [0, 1, 0, 0], [0, 0, 1, 0], [1, 0, 0, 0]])
			]

	p4	=	[np.matrix([[1, 0, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0], [0, 1, 0, 0]]),
			np.matrix([[0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1], [1, 0, 0, 0]]),
			np.matrix([[0, 0, 1, 0], [0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 0, 1]]),
			np.matrix([[0, 0, 0, 1], [1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0]])
			]

	p5	=	[np.matrix([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]]),
			np.matrix([[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]),
			np.matrix([[0, 0, 1, 0], [0, 0, 0, 1], [0, 1, 0, 0], [1, 0, 0, 0]]),
			np.matrix([[0, 0, 0, 1], [0, 0, 1, 0], [1, 0, 0, 0], [0, 1, 0, 0]])
			]

	p6	=	[np.matrix([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]),
			np.matrix([[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]]),
			np.matrix([[0, 0, 1, 0], [0, 0, 0, 1], [1, 0, 0, 0], [0, 1, 0, 0]]),
			np.matrix([[0, 0, 0, 1], [0, 0, 1, 0], [0, 1, 0, 0], [1, 0, 0, 0]])
			]

	pie_complete = [ p1, p2, p3, p4, p5, p6 ]

	return pie_complete[pie_index]


##************************************
# Perform flop operation over Adinkra color space
def colorspace_flop(adinkra, flop_op):
	# print("Executing colorspace_flip", flop_op)
	""" Moving around the Lmatrices in the Adinkra	"""
	new_adinkra	= [ adinkra[(ind - 1)] for ind in flop_op]

	return new_adinkra


##************************************
# Perform flip operation over Adinkra color space
def colorspace_flip(adinkra, flip_op):

	# print("Executing colorspace_flip", flip_op)
	""" Weird bug here if you do the algorith this way """
	# new_adinkra	= []
	# for i in range(0, len(adinkra)):
	# 	tmat 		= adinkra[i]
	# 	tmat[1:]	= tmat[1:] * flip_op[i]
	# 	new_adinkra.append(tmat)
	""" Normal algorithm """
	new_adinkra	= []
	for i in range(0, len(flip_op)):
		# new_adinkra[i][1:]	= new_adinkra[i][1:] * flip_op[i]
		temp_mat	= adinkra[i] * flip_op[i]
		new_adinkra.append(temp_mat)
	return new_adinkra


##************************************
# Defining the binary multiplication arrays
def binaries(bin_code):
	bin_code = int(bin_code)
	binaries_lt	= [(0, [1,1,1,1]), (2, [1,-1,1,1]), (4, [1,1,-1,1]),
					(6, [1,-1,-1,1]), (8, [1,1,1,-1]), (10, [1,-1,1,-1]),
					(12, [1,1,-1,-1]), (14, [1,-1,-1,-1])]

	for btuple in binaries_lt:
		if bin_code == btuple[0]:
			return btuple[1]

##************************************
# Defining the flips used in original CG BC4 Library
def flips_org_lib(flip_set):

	p_slice 	= {}

	# p_slice['P1']	= [[0,12,10,6], [0,6,12,10], [2,4,14,8], [2,14,8,4],
	# 					[8,4,2,14], [8,14,4,2], [10,12,6,0], [10,6,0,12]]
	p_slice['P1']	= [[0,6,12,10], [0,12,10,6], [2,4,14,8], [2,14,8,4],
						[4,2,8,14], [4,8,14,2], [6,0,10,12], [6,10,12,0],
						[8,4,2,14], [8,14,4,2], [10,6,0,12], [10,12,6,0],
						[12,0,6,10], [12,10,0,6], [14,2,4,8], [14,8,2,4]]

	p_slice['P2']	= [[0,10,6,12], [0,12,10,6], [2,8,4,14], [2,14,8,4],
						[4,8,14,2], [4,14,2,8], [6,10,12,0], [6,12,0,10],
						[8,2,14,4], [8,4,2,14], [10,0,12,6], [10,6,0,12],
						[12,0,6,10], [12,6,10,0], [14,2,4,8], [14,4,8,2]]

	p_slice['P3']	= [[0,6,10,12], [0,12,6,10], [2,4,8,14], [2,14,4,8],
						[4,2,14,8], [4,8,2,14], [6,0,12,10], [6,10,0,12],
						[8,4,14,2], [8,14,2,4], [10,6,12,0], [10,12,0,6],
						[12,0,10,6], [12,10,6,0], [14,2,8,4], [14,8,4,2]]

	p_slice['P4']	= [[0,10,12,6], [0,12,6,10], [2,8,14,4], [2,14,4,8],
						[4,8,2,14], [4,14,8,2], [6,10,0,12], [6,12,10,0],
						[8,2,4,14], [8,4,14,2], [10,0,6,12], [10,6,12,0],
						[12,0,10,6], [12,6,0,10], [14,2,8,4], [14,4,2,8]]

	p_slice['P5']	= [[0,6,10,12], [0,10,12,6], [2,4,8,14], [2,8,14,4],
						[4,2,14,8], [4,14,8,2], [6,0,12,10], [6,12,10,0],
						[8,2,4,14], [8,14,2,4], [10,0,6,12], [10,12,0,6],
						[12,6,0,10], [12,10,6,0], [14,4,2,8], [14,8,4,2]]

	p_slice['P6']	= [[0,6,12,10], [0,10,6,12], [2,4,14,8], [2,8,4,14],
						[4,2,8,14], [4,14,2,8], [6,0,10,12], [6,12,0,10],
						[8,2,14,4], [8,14,4,2], [10,0,12,6], [10,12,6,0],
						[12,6,10,0], [12,10,0,6], [14,4,8,2], [14,8,2,4]]

	return p_slice[flip_set]


##************************************
# Defining the elle binary representations for the Vierergruppe
def flip_ellebin(flip_set):

	vgrp_elle			= {}

	vgrp_elle['()']		= [[14,8,2,4], [2,4,14,8], [4,2,8,14],[8,14,4,2],
							[6,0,10,12], [10,12,6,0], [12,10,0,6], [0,6,12,10]]

	vgrp_elle['(12)'] 	= [[14,4,2,8], [2,8,14,4], [4,14,8,2], [8,2,4,14],
							[6,12,10,0], [10,0,6,12], [12,6,0,10], [0,10,12,6]]

	vgrp_elle['(13)'] 	= [[14,2,8,4], [2,14,4,8], [4,8,2,14], [8,4,14,2],
							[6,10,0,12], [10,6,12,0], [12,0,10,6], [0,12,6,10]]

	vgrp_elle['(23)'] 	= [[2,4,8,14], [14,8,4,2], [8,14,2,4], [4,2,14,8],
							[10,12,0,6], [6,0,12,10], [0,6,10,12], [12,10,6,0]]

	vgrp_elle['(123)']	= [[14,4,8,2], [2,8,4,14], [4,14,2,8], [8,2,14,4],
							[6,12,0,10], [10,0,12,6], [12,6,10,0], [0,10,6,12]]

	vgrp_elle['(132)']	= [[14,2,4,8], [2,14,8,4], [4,8,14,2], [8,4,2,14],
							[6,10,12,0], [10,6,0,12], [12,0,6,10], [0,12,10,6]]

	return vgrp_elle[flip_set]


##************************************
# Defining the six Vierergruppe representations
def vierergruppe_sets():

	""" Defining V and the rest. These are elements for flopping bosonic fields
	"""
	vp1 	= np.matrix([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
	vp2 	= np.matrix([[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]])
	vp3 	= np.matrix([[0, 0, 1, 0], [0, 0, 0, 1], [1, 0, 0, 0], [0, 1, 0, 0]])
	vp4 	= np.matrix([[0, 0, 0, 1], [0, 0, 1, 0], [0, 1, 0, 0], [1, 0, 0, 0]])
	vprime 	= [vp1, vp2, vp3, vp4]

	b12 	= np.matrix([[0,1,0,0], [1,0,0,0], [0,0,1,0], [0,0,0,1]])
	b13 	= np.matrix([[0,0,1,0], [0,1,0,0], [1,0,0,0], [0,0,0,1]])
	b23 	= np.matrix([[1,0,0,0], [0,0,1,0], [0,1,0,0], [0,0,0,1]])
	b123	= np.matrix([[0,0,1,0], [1,0,0,0], [0,1,0,0], [0,0,0,1]])
	b132	= np.matrix([[0,1,0,0], [0,0,1,0], [1,0,0,0], [0,0,0,1]])

	""" Vierergrupe via S3 left cosets """
	vgruppe	= 	{'()': vprime,
				'(12)': [np.dot(b12, i) for i in vprime],
				'(13)': [np.dot(b13, i) for i in vprime],
				'(23)': [np.dot(b23, i) for i in vprime],
				'(123)':[np.dot(b123, i) for i in vprime],
				'(132)':[np.dot(b132, i) for i in vprime]
				}

	return vgruppe


##************************************
# Defining the Vierergruppe representations for flop operations
def vierergruppe_flops():

	vprime 	= [ "()", "(12)(34)", "(13)(24)", "(14)(23)" ]
	# vprime	= [ "1234", "2143", "3412", "4321"]
	vprime	= [ [1,2,3,4], [2,1,4,3], [3,4,1,2], [4,3,2,1] ]


	vgrpv 		= [ ("()", [1,2,3,4]), ("(12)(34)", [2,1,4,3]),
					("13)(24)", [3,4,1,2]), ("(14)(23)", [4,3,2,1])]

	vgrp12v		= [ ("(12)", [2,1,3,4]), ("(34)", [1,2,4,3]),
					# ("(1324)", [3,4,2,1]), ("(1423)", [4,3,1,2]) ]
					("(1324)", [4,3,1,2]), ("(1423)", [3,4,2,1]) ]

	vgrp13v		= [ ("(13)", [3,2,1,4]), ("(1234)", [4,1,2,3]),
					("(24)", [1,4,3,2]), ("(1432)", [2,3,4,1])	]

	vgrp23v		= [ ("(23)", [1,3,2,4]), ("(1342)", [2,4,1,3]),
					("(1243)", [3,1,4,2]), ("(14)", [4,2,3,1])	]

	vgrp123v	= [ ("(132)", [3,1,2,4]), ("(132)", [4,2,1,3]),
					("(243)", [1,3,4,2]), ("(142)", [2,4,3,1])	]

	vgrp132v	= [ ("(132)", [2,3,1,4]), ("(234)", [1,4,2,3]),
					("(124)", [4,1,3,2]), ("(143)", [3,2,4,1])	]


	vgruppe	   = [ ('()', vgrpv), ('(12)',vgrp12v), ('(13)', vgrp13v),
				('(23)', vgrp23v), ('(123)', vgrp123v), ('(132)', vgrp132v) ]

	return vgruppe


# **************************************************************************
# Use the binary representation info to perform flips on L mats in each tetrad
def lmat_flipping(vbasis, binaries_list):

	lmat_list		= []

	for xbin in binaries_list:
		print("Flip:", xbin)
		binmats = [binaries(b) for b in xbin]
		temp	= [np.dot(binmats[i], vbasis[i]) for i in range(0, len(binmats))]
		lmat_list.append(temp)

	return lmat_list


# **************************************************************************
# Check user input for issues
def verify_input(userstr):

	'''
		Check for PALL and for P1 - P6
	'''
	# restr = re.compile('[pP1-6]', re.IGNORECASE)
	# str2 = restr.findall(userstr)
	loc_userstr	= userstr.lower()
	if userstr.isdigit():
		if int(userstr) in list(range(1,7)):
			return "P" + userstr
		else:
			return 0
	elif userstr[0].isdigit():
		if int(userstr[0]) in list(range(1,7)):
			return "P" + userstr
		else:
			return 0

	if (loc_userstr).startswith('p'):
		recomp = re.compile('[pP1-6]', re.IGNORECASE)
		relist = recomp.findall(userstr)
		if recomp.match(relist[0]) is not None:
			if len(relist) > 1:
				if relist[1].isdigit():
					print('Sucesfull BC4 CG Library matching input')
					return 'P' + relist[1]
		else:
			sys.exit('STRING ISSUE')
	elif not (loc_userstr).startswith('p'):
		ustr_len = len(userstr)
		recomp = re.compile('[1-6]')
		if userstr.isdigit():
			if recomp.match(userstr) is not None:
				if ustr_len == 1:
					print(userstr)
					return userstr
				elif ustr_len != 1:
					print("Digit error")
					print("")
					sys.exit('Add more crap here')
			else:
				print("NONE ISSUE")


# **************************************************************************
# Adding calculation options
def calc_options():

	base_options_print()
	uinput	= input("Choose wisely: ")

	if int(uinput) == 1:
		print("")
		print("Would you like to calculate entire BC4 CG Small Library?")
		print("")
		userlib = input(" yes/no ")
		if userlib == 'yes':
			pass
		else:
			userchk = input("exit?  yes/no?")
			if userchk.lower() == 'yes':
				sys.exit("EXITING")
			else:
				base_options_print()
				uinput	= input("Choose wisely: ")
				''' *** ** * ** *** **** *** ** * ** ***
					FIGURE OUT HOW TO MAKE THESE LOOP
					*** ** * ** *** **** *** ** * ** *** '''
				pass

	elif int(uinput) == 2:
		print("")
		print("Pick P-set to calculate")
		print("P1 or 1 - {P1} ")
		print("P2 or 1 - {P2} ")
		print("P3 or 3 - {P3} ")
		print("P4 or 4 - {P4} ")
		print("P5 or 5 - {P5} ")
		print("P6 or 6 - {P6} ")
		print("PALL - All P sets")
		print("")
		userpset = input(":")
		checkstr = verify_input(userpset)

	elif int(uinput) == 3:
		print("")
		print("Pick P-set to calculate")
		print("P1 or 1 - {P1} ")
		print("P2 or 1 - {P2} ")
		print("P3 or 3 - {P3} ")
		print("P4 or 4 - {P4} ")
		print("P5 or 5 - {P5} ")
		print("P6 or 6 - {P6} ")
		print("PALL - All P sets")
		print("")
		userpset = input(":")
		checkstr = verify_input(userpset)
	elif int(uinput) == 4:
		print("")
		print("Please enter file naming convention")
		fileconv = input("")
		pass
	elif int(uinput) == 5 or uinput.lower() == 'exit':
		print("")
		sys.exit("EXITING")
	else:
		print("Please re-enter selection or type 'exit'")
		uinput = input("")


# **************************************************************************
# Default standard P set options
def pset_options_std():
	''' Making a func of printing the P-set options. Tired of rewriting it '''
	userinput = ""
	print("")
	print("Pick P-set to calculate. Use 1-6 or P1 - P6")
	print(" < 1 >  -  {P1} ")
	print(" < 2 >  -  {P2} ")
	print(" < 3 >  -  {P3} ")
	print(" < 4 >  -  {P4} ")
	print(" < 5 >  -  {P5} ")
	print(" < 6 >  -  {P6} ")
	print(" < 7 / exit >  -  Go back")
	userinput = input(": ")
	checkstr  = verify_input(userinput)
	if not checkstr:
		print("INVALID INPUT")
		print("TRY AGAIN")
		pset_options_std()
		pass
	elif checkstr.strip() == '7' or checkstr.lower() == 'exit':
		user_options()
	elif checkstr:
		return userinput
	else:
		pass


# **************************************************************************
# User options
def user_options():

	print("#***********************************************************************")
	print("# Name: BC4 Coxeter Group Small Library Toolset ")
	print("# Author:  Vadim Korotkikh	")
	print("# Date:    June 2017	")
	print("# Description: Function for verifying the BC4 space coefficient library")
	print("#	")
	# print("# Description: Function for verifying the BC4 space coefficient library")
	# print("#	")
	print("#***********************************************************************")
	print("	")

	# **************************************************************************
	def core_options():
		print("Choose from one of the following calculation options:")
		print("")
		print(" < 1 >  -  Calculate P-set Bosonic Matrices")
		print(" < 2 >  -  Calculate P-set Fermionic Matrices")
		print(" < 3 >  -  Display/Print BC4 CG Small Library P-sets")
		print(" < 4 >  -  Verify BC4 CG Small Library ~V coefficient values")
		print(" < 5 >  -  Exit")
		# print(" < 1 >  -  Display BC4 CG Small Library P sets")
		# print(" < 2 >  -  Verify BC4 CG Small Library ~V coefficient values")
		# print(" < 3 >  -  Calculate P-set Fermionic Matrices")
		# print(" < 4 >  -  Calculate P-set Bosonic Matrices")
		# print(" < 5 >  -  Set output file string")
		# print(" < 6 >  =  Nevermind. Get me outa here! Exit")
		print("")
		return input(": ")

	# Set loopcount = 0 of no arg is supplied for first time
	def option_one(loopcount=0):
		print("")
		print(" < 1 >  -  Calculate all P-sets")
		print(" < 2 >  -  Calculate select P-set from the Small Library")
		print(" < 3 >  -  Back to main menu")
		ninput = input(": ")
		if ninput.strip() == '1':
			# bc4_validation_organizer('PALL', 'Vmats', 'fermi')
			bc4_validation_organizer('PALL', 'mats', 'boson')
			pass
		elif ninput.strip() == '2':
			usr_pset = pset_options_std()
			print(usr_pset)
			bc4_validation_organizer(usr_pset, 'mats', 'boson')
		elif ninput.strip() == '3':
			option_activator('core')
		else:
			loopcount += 1
			print("Unrecognized option")
			if loopcount <= 5:
				option_one(loopcount)
			else:
				print("Returning to core options...")
				option_activator('core')

	# Set loopcount = 0 of no arg is supplied for first time
	def option_two(loopcount=0):
		print("")
		print(" < 1 >  -  Calculate all P-sets")
		print(" < 2 >  -  Calculate select P-set from the Small Library")
		print(" < 3 >  -  Back to main menu")
		ninput = input(": ")
		if ninput.strip() == '1':
			# bc4_validation_organizer('PALL', 'Vmats', 'fermi')
			bc4_validation_organizer('PALL', 'mats', 'fermi')
			pass
		elif ninput.strip() == '2':
			usr_pset = pset_options_std()
			print(usr_pset)
			bc4_validation_organizer(usr_pset, 'mats', 'fermi')
		elif ninput.strip() == '3':
			option_activator('core')
		else:
			loopcount += 1
			print("Unrecognized option")
			if loopcount <= 5:
				option_two(loopcount)
			else:
				print("Returning to core options...")
				option_activator('core')

	def option_three(loopcount=0):
		counter	= 0
		print("")
		print(" < 1 >  -  Display entire BC4 CG Small Library")
		print(" < 2 >  -  Display select P-set from the Small Library")
		print(" < 3 >  -  Back to main menu")
		ninput = input(": ")
		if ninput.strip() == '1':
			pass
		elif ninput.strip() == '2':
			usr_pset = pset_options_std()
			# Lets make the P-Set options a function
		elif ninput.strip() == '3':
			option_activator('core')
		else:
			loopcount += 1
			print("Unrecognized option")
			if loopcount <= 5:
				option_three(loopcount)
			else:
				print("Returning to core options...")
				option_activator('core')

	def option_four(loopcount=0):
		print("")
		print(" < 1 >  -  Verify entire BC4 CG Library")
		print(" < 2 >  -  Verify select P-set from the Small Library")
		print(" < 3 >  -  Back to main menu")
		ninput = input(": ")
		if ninput.strip() == '1':
			# for now
			pass
		elif ninput.strip() == '2':
			usr_pset = pset_options_std()
			print(usr_pset)
		elif ninput.strip() == '3':
			option_activator('core')
		else:
			loopcount += 1
			print("Unrecognized option")
			if loopcount <= 5:
				option_four(loopcount)
			else:
				print("Returning to core options...")
				option_activator('core')

	# Preset loopcount = 0 if no arg supplied for first time.
	# def option_four(loopcount=0):
	# 	print("")
	# 	print(" < 1 >  -  Calculate all P-sets")
	# 	print(" < 2 >  -  Calculate select P-set from the Small Library")
	# 	print(" < 3 >  -  Back to main menu")
	# 	opt_str = input(": ")
	#
	# 	if opt_str.strip() == '1':
	# 		# for now
	# 		bc4_validation_organizer('PALL', 'mats', 'boson')
	# 		pass
	# 	elif opt_str.strip() == '2':
	# 		usr_pset = pset_options_std()
	# 		bc4_validation_organizer(usr_pset, 'mats', 'boson')
	# 	elif opt_str.strip() == '3':
	# 		option_activator('core')
	# 	else:
	# 		loopcount += 1
	# 		print("Unrecognized option")
	# 		if loopcount <= 5:
	# 			option_four(loopcount)
	# 		else:
	# 			print("Returning to core options...")
	# 			option_activator('core')

	def option_five():
		print("NOT ACTIVATED (code not finished)")
		print("Going back to main menu")
		option_activator('core')

	def option_six():
		print("")
		print("Quiting script. Are you sure (yes/no)?")
		ninput = input(": ")
		if ninput.lower() == 'yes':
			sys.exit("EXITING BC4 CG Library Utility")
		elif ninput.lower() == 'no':
			print("Going back to main menu")
			option_activator('core')
		else:
			pass

	# **************************************************************************
	# Executes options inner functions based on input_str
	def option_activator(input_str, mcounter=0):
		''' option_activator - Inner function of user_options that depending on
			input_str either executes other inner functions or rexecutes itself
			with core_options() func. providing input string.
		'''
		if input_str.strip() == '1':
			option_one()
		elif input_str.strip() == '2':
			option_two()
		elif input_str.strip() == '3':
			option_three()
		elif input_str.strip() == '4':
			option_four()
		elif input_str.strip() == '5':
			option_five()
		elif input_str.strip() == '6' or input_str.lower() == 'exit':
			option_six()
		elif input_str.strip() == 'core' or input_str.lower() == 'core':
			option_activator(core_options())
		else:
			print("UNRECOGNIZED SELECTION: ", input_str)
			print("Try again...")
			if mcounter <= 6:
				mcounter+=1
				option_activator(core_options(), mcounter)
			elif mcounter >= 7:
				print("Too many attempts. Try again later.")
				# print("EXITING BC4 CG Library Utility")

				sys.exit("EXITING BC4 CG Library Utility")

	uinput = core_options()
	option_activator(uinput)


# **************************************************************************
# Execute main()
if __name__ == "__main__":
	start_time = time.time()

	pset_def = "P1"
	try:
		main(sys.argv[1])
	except IndexError:
		# calc_options()
		user_options()

	print("-- Execution time --")
	print("---- %s seconds ----" % (time.time() - start_time))
