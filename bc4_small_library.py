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
	bc4_validation_seq(pset_str)

# ******************************************************************************
# BC4 Validation function process organizer
def bc4_validation_seq(pset_arg):

	psets_list		= ["P1", "P2", "P3", "P4", "P5", "P6"]
	psets_dict		= {}
	'''
		Define Holoraumy matrix type
	'''
	holotype	= "bosonic"

	""" Generate individual Library slices	"""
	for ps in psets_list:
		temp_tetrad	= tetrad_setgen(ps)
		if ps not in psets_dict:
			psets_dict['%s' % ps] = temp_tetrad
		else:
			print("Unknown P set -  ERROR")

	# """ Changing PALL to ALL to calculate for one big chunk	"""
	if pset_arg == "ALL":
		temp_plist	= []
		psl 		= psets_list
		temp_plist	= tetrad_setgen(pset_arg)
	# 	print("# ********************************")
	# 	print("Execute Holoraumy calc for P sets:", psl[0], psl[1], psl[2], psl[3], psl[4], psl[5])
	# 	print("		")
	# 	vij_holoraumy_calc.calc_holoraumy_mats(temp_plist, pset_arg, holotype)

	elif pset_arg in psets_list:
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
	pint 	= int(pset.lstrip("P")) - 1
	pslice 	= lib_pslices(pint)
	if p_switch:
		print("# ********************************")
		print("Starting conversion process for P slice:", pset)
		# print("Length of", pset, "tetrad set:", len(run_group))
		print("")
	""" Perform boolean calculations """
	pbools	= flips_org_lib(pset)
	# print(pbools)
	# print("pbools type", type(pbools))

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

##************************************
# Defining the P slices of original BC4 CG library
def lib_pslices(pie_index):

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

	p_slices = [ p1, p2, p3, p4, p5, p6 ]

	return p_slices[pie_index]

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


##************************************
# Function for calling calculate_vij_matrices
def calculate_vgruppe_sets(gruppe_adinkras, gruppe_binaries):
	"""
	Function for printing out the details of each Adinkra - Vij matrix
	sixset calculation, including the binary representation and the
	corresponding resulting Vij matrices and their elles/tilde elles
	Coefficients.
	"""

	for i in range(0, len(gruppe_adinkras)):

		print("")
		print("Calculating for binary flip:", gruppe_binaries[i])

		vijset = vij_holoraumy_4x4.calculate_vijmatset(gruppe_adinkras[i])
		for i in vijset:
			print(i)
	# for vgrp, binaries_list in vierergruppe_tilde.items():
	# 	vbasis	= vgruppe_sets[vgrp]
	# 	temp 	= lmat_flipping(vbasis, binaries_list)
	# 	# for i, tet in enumerate(temp):
	# 	# 	print("Length of tet:", len(tet), "Type", type(tet))
	# 	# print("Length lmat_flipping", len(temp), vgrp, binaries_list)
	# 	main_tetrad.extend(temp)


##************************************
# Use the binary representation info to perform flips on L mats in each tetrad
def lmat_flipping(vbasis, binaries_list):

	lmat_list		= []

	for xbin in binaries_list:
		print("Flip:", xbin)
		binmats = [binaries(b) for b in xbin]
		temp	= [np.dot(binmats[i], vbasis[i]) for i in range(0, len(binmats))]
		lmat_list.append(temp)

	return lmat_list


# ********************************
# Alpha and Beta matrices hardcoded
def alphas_betas():

	""" These are the alpha and beta matrices multiplied by 2i
		alpha and beta are tensor products of Pauli spin matrices
 		Identity matrix They are defined in equtions (4.5)
		in Isaac Chappell II, S. James Gates, Jr - 2012
	"""

	alpha1i = np.matrix([[0, 0, 0, 2], [0, 0, 2, 0], [0, -2, 0, 0], [-2, 0, 0, 0]])
	alpha2i = np.matrix([[0, 2, 0, 0], [-2, 0, 0, 0], [0, 0, 0, 2], [0, 0, -2, 0]])
	alpha3i = np.matrix([[0, 0, 2, 0], [0, 0, 0, -2], [-2, 0, 0, 0], [0, 2, 0, 0]])

	# Betas
	beta1i = np.matrix([[0, 0, 0, 2], [0, 0, -2, 0], [0, 2, 0, 0], [-2, 0, 0, 0]])
	beta2i = np.matrix([[0, 0, 2, 0], [0, 0, 0, 2], [-2, 0, 0, 0], [0, -2, 0, 0]])
	beta3i = np.matrix([[0, 2, 0, 0], [-2, 0, 0, 0], [0, 0, 0, -2], [0, 0, 2, 0]])

	return [alpha1i, alpha2i, alpha3i, beta1i, beta2i, beta3i]


# **************************************************************************
# Execute main()
if __name__ == "__main__":
	start_time = time.time()

	pset_def = "P1"

	try:
		main(sys.argv[1])
	except IndexError:
		print("Using hardcoded default P-set", pset_def)
		main(pset_def)
	print("-- Execution time --")
	print("---- %s seconds ----" % (time.time() - start_time))
