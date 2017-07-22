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
# import matrix_outerprod_calc
import matrix_calc_vijmat2
import matrix_calc_vijmat_nopr
import vij_holoraumy_calc

# ******************************************************************************
# Main() function.
def main():

	print("# ***********************************************************************")
	print("# Name:    Verifying BC4 Space Vij coefficient")
	print("# Author:  Vadim Korotkikh	")
	print("# Email:   va.korotki@gmail.com")
	print("# Date:    November 2016		")
	print("# Version: 1.0A				")
	print("#							")
	print("# Description: Function for verifying the BC4 space coefficient library")
	print("#	")
	print("# ***********************************************************************")
	print("		")

	""" PALL seems to be broken for now due to lack of if elif logic at end
		of string_to_tetrad() function """
	""" Options for BC4 Library:
		ALL - Entire small library, one big chunk.
		PALL - Entire small library, one at a time.
		P1, P2, P3, P4, P5, P6 - Only one section. 	"""
	# pset_str = "PALL"
	pset_str = "PALL"
	bc4_validation_seq(pset_str)

# ******************************************************************************
# BC4 Validation function process organizer
def bc4_validation_seq(pset_arg):

	psets_list		= ["P1", "P2", "P3", "P4", "P5", "P6"]
	psets_dict		= {}

	""" Generate individual Library slices	"""
	for ps in psets_list:
		temp_tetrad	= tetrad_setgen(ps)
		if ps not in psets_dict:
			psets_dict['%s' % ps] = temp_tetrad
		else:
			print("Unknown P set -  ERROR")

	# temp_l			= []
	""" Changing PALL to ALL to calculate for one big chunk	"""
	if pset_arg == "ALL":
		temp_plist	= []
		psl 		= psets_list
		temp_plist	= tetrad_setgen(pset_arg)
		# for p, plist in psets_dict.items():
		# 	temp_plist.extend(plist)

		print("# ********************************")
		print("Execute Gadget calc for P sets:", psl[0], psl[1], psl[2], psl[3], psl[4], psl[5])
		print("		")
		matrix_calc_vijmat2.calculate_vij_matrices(temp_plist)

	elif pset_arg in psets_list:
	# elif pset_arg != "PALL":
		temp_plist = psets_dict[pset_arg]
		print("		")
		print("Execute Gadget calc for P-set:", pset_arg)
		print("		")
		matrix_calc_vijmat2.calculate_vij_matrices(temp_plist)
		print("Gadget calc. for:", pset_arg,"finished")
		print("")
	elif pset_arg == "PALL":
		for pset, plist in sorted(psets_dict.items()):
			pint = 0
			pint = int(pset.lstrip("P"))
			print("		")
			print("Execute Gadget calc for P-set:", pset)
			print("		")
			matrix_calc_vijmat_nopr.calculate_vij_matrices(plist)
			print("Gadget calc. for:", pset,"finished")
			print("")
		""" Code for doing Gadget value calc for P - Pairs 	"""
			# for pxset, txlist in psets_dict.items():
			# 	pxint = 0
			# 	pxint = int(pxset.lstrip("P"))
			#
			# 	combonum = [pint, pxint]
			# 	# combonum.sort()
			# 	if combonum not in temp_l:
			# 	# if combonum == [2,1]:
			# 		temp_l.append(combonum)
			# 		# print(combonum)
			# 		print("		")
			# 		print("Execute Gadget calc for pair:", pset, pxset)
			# 		print("		")
			# 		temp_pspx = []
			# 		temp_pspx.extend(tlist)
			# 		temp_pspx.extend(txlist)
			# 		print("Print length of Pi-Pj set", len(temp_pspx))
			# 		# klein_check(tlist, txlist)
			# 		# print("Klein check finished for:", pset,"-",pxset)
			# 		matrix_calc_vijmat2.calculate_vij_matrices(temp_pspx)
			# 		print("Gadget calc. for:", pset,"-",pxset, "finished")
			# 		print("")
			# 	else:
			# 		pass
	# temp_l.sort()
	#
	# print("Printing length of P lists")
	# print(len(psets_dict))
	# print(temp_l)
	# matrix_calc_vijmat2.calculate_vij_matrices(temp_tetrad)


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
# Compiling the tetrads from predfined Adinkras
def assemble_tetrads():

	"""Vierergrupe dictionaries with binary quadsets for ells and tilde-ells
	"""

	main_tetrad			= []

	vgruppe_sets		= vierergruppe_sets()
	vierergruppe_elle	= flip_ellebin()
	vierergruppe_tilde	= flip_tildebin()

	for vgrp, binaries_list in vierergruppe_elle.items():
		vbasis	= vgruppe_sets[vgrp]
		print("")
		print("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx ")
		print("Calculating Vij elle coefficients")
		print("							")
		print("Vierergruppe flop: ",vgrp)
		temp 	= lmat_flipping(vbasis, binaries_list)
		# print("Flip sets:", binaries_list)
		# vij_holoraumy_prime.calculate_vij_matrices(temp)
		calculate_vgruppe_sets(temp, binaries_list)

		print("<<<>>>")
		main_tetrad.extend(temp)

##************************************
# Defining the six Vierergruppe representations
def vierergruppe_sets():

	vp1 	= np.matrix([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
	vp2 	= np.matrix([[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]])
	vp3 	= np.matrix([[0, 0, 1, 0], [0, 0, 0, 1], [1, 0, 0, 0], [0, 1, 0, 0]])
	vp4 	= np.matrix([[0, 0, 0, 1], [0, 0, 1, 0], [0, 1, 0, 0], [1, 0, 0, 0]])
	vprime 	= [vp1, vp2, vp3, vp4]

	"""Elements for flopping bosonic fields"""
	b12 	= np.matrix([[0,1,0,0], [1,0,0,0], [0,0,1,0], [0,0,0,1]])
	b13 	= np.matrix([[0,0,1,0], [0,1,0,0], [1,0,0,0], [0,0,0,1]])
	b23 	= np.matrix([[1,0,0,0], [0,0,1,0], [0,1,0,0], [0,0,0,1]])
	b123	= np.matrix([[0,0,1,0], [1,0,0,0], [0,1,0,0], [0,0,0,1]])
	b132	= np.matrix([[0,1,0,0], [0,0,1,0], [1,0,0,0], [0,0,0,1]])

	vgruppe	= 	{'()': vprime,
				'(12)': [np.dot(b12, i) for i in vprime],
				'(13)': [np.dot(b13, i) for i in vprime],
				'(23)': [np.dot(b23, i) for i in vprime],
				'(123)':[np.dot(b123, i) for i in vprime],
				'(132)':[np.dot(b132, i) for i in vprime]
				}

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

# ******************************************************************************
def tetrad_setgen(pset):
	"""	Generates a {P#} set of tetrads via execution of pset_string_format()
		and	string_to_tetrad() function. Creates a set of python numpy tetrads
		from string representation	"""

	p_switch = 0
	run_group = pset_string_format(pset)
	if pset == "PALL":
		pset = "P1, P2, P3, P4, P5, P6"
	if p_switch:
		print("# ********************************")
		print("Starting conversion process for", pset)
		print("Length of", pset, "tetrad set:", len(run_group))
		print("")
		print("Printing", pset, "tetrads string representation before conversion")
		for x in run_group:
			x = re.sub(r"\[", "<", x)
			x = re.sub(r"\]", ">", x)
			print(x)
		print("")

	tetrad_sig_perm = gen_sign_perm(4)
	tsp = tetrad_sig_perm
	tetrad_bc4chk =[]

	for ind, itet in enumerate(run_group):
		# print("IND DEBUG:", ind)
		tempt = string_to_tetrad(pset, ind, itet)
		tetrad_bc4chk.append(tempt)

	return tetrad_bc4chk


# ******************************************************************************
# Transform a tetrad string representation into tetrad of matrices
def string_to_tetrad(p_str,indx_num,tet_strrep):

	""" Turn debug_pr = 1 or True to turn on print output of
	string to tetrad conversion process. Simple shortcut for now """
	debug_pr 	= 0

	qt_temp		= []
	dec_indx = indx_num * 4
	tet_rep = tet_strrep.split()
	t1 = np.array((1,0,0,0))
	t2 = np.array((0,1,0,0))
	t3 = np.array((0,0,1,0))
	t4 = np.array((0,0,0,1))
	vtl = [t1,t2,t3,t4]
	# onesl	= np.ones(4, int)
	# vtl		= np.diag(onesl)
	tetint_list	= []
	if debug_pr:
		tet_nicerep = re.sub(r"\[", "<", tet_strrep)
		tet_nicerep2 = re.sub(r"\]", ">", tet_nicerep)
		print("# ********************************")
		print("Converting tetrad #:", indx_num)
		print("Str Representation: ", tet_nicerep2)
		print("")
	elif debug_pr == 0:
		pass
	for i, m in enumerate(tet_rep):
		xstr = re.sub(r"\[", "<", m)
		xstr = re.sub(r"\]", ">", xstr)
		xstr = xstr.strip(",")
		mtemp = m[m.index("[") + 1:m.rindex("]")]
		# print(mtemp)
		mint_list = [int(s) for s in mtemp.split(',')]
		tetint_list.append(mint_list)
		# tf = any(not isinstance(x, int) for x in mint_list)
		sign_ind = [ x / abs(x) for x in mint_list]
		sgi = [ x / abs(x) for x in mint_list]
		if debug_pr:
			print(mint_list)
			print(sign_ind)
		else:
			pass
		mi = [abs(x) - 1 for x in mint_list]
		if debug_pr:
			print(mi)
			print("")
		else:
			pass
		# tempm = numpy.column_stack((vtl[mi[0]]*sgi[0],vtl[mi[1]]*sgi[1],vtl[mi[2]]*sgi[2],vtl[mi[3]]*sgi[3]))
		tempm = numpy.vstack((vtl[mi[0]]*sgi[0],vtl[mi[1]]*sgi[1],vtl[mi[2]]*sgi[2],vtl[mi[3]]*sgi[3]))
		matint = tempm.astype(int)
		matint1 = np.asmatrix(matint)
		if debug_pr:
			print("String:", xstr)
			print(matint1)
		else:
			pass
		qt_temp.append((i+dec_indx, matint1))
	# qt_temp[0][1] =  np.multiply(qt_temp[0][1], -1)
	f_temp 	= []
	if p_str == "P1":
		f_temp = qt_temp
	elif p_str == "P2":
		f_temp = qt_temp
	elif p_str == "P3":
		f_temp = [qt_temp[1], qt_temp[0], qt_temp[3], qt_temp[2]]
	elif p_str == "P4":
		# f_temp = [qt_temp[3], qt_temp[0], qt_temp[1], qt_temp[2]]
		f_temp = qt_temp
	elif p_str == "P5":
		# f_temp = [qt_temp[2], qt_temp[3], qt_temp[1], qt_temp[0]]
		f_temp = qt_temp
	elif p_str == "P6":
		# f_temp = [qt_temp[2], qt_temp[3], qt_temp[0], qt_temp[1]]
		# f_temp = [qt_temp[2], qt_temp[3], qt_temp[1], qt_temp[0]]
		f_temp = qt_temp
	return f_temp
	# return qt_temp

# ******************************************************************************
# Find all patterns within the list of tetrads.
def klein_check(tet_lista, tet_listb):

	self_kein	= []
	kein_flip	= []

	pair_match	= []
	print("Printing Tetrads plugged into klein_check")
	print("")

	for ind, itet in enumerate(tet_lista):
		for i in itet:
			it = np.transpose(i[1])
		ivt = [np.transpose(xm[1]) for xm in itet]

		for jnd, jtet in enumerate(tet_listb):

			print("Tetrad #: ", ind, jnd)
			jtet_is = [jm[1] for jm in jtet]
			# if ind != jnd:
			# print(ivt[0], jtet[0][1])
			print("")
			if np.array_equal(ivt[0], jtet[0][1]) and np.array_equal(ivt[1], jtet[1][1]):
				if np.array_equal(ivt[2], jtet[2][1]) and np.array_equal(ivt[3], jtet[3][1]):
					print("Klein Flip found")
					print("", ind, jnd)
					demp = [ind, jnd]
					demp.sort()
					if demp not in kein_flip:
						kein_flip.append(demp)
					else:
						print("Duplicate Klein")
						pass
			elif any( mj for im in ivt for mj in jtet if np.array_equal(im, mj[1])):
				print("Self-Klein Flip found")
				print(itet[0], jtet[0])
				print("", ind, jnd)
				demp = [ind, jnd]
				demp.sort()
				if demp not in self_kein:
					self_kein.append(demp)
				else:
					print("Duplicate Self-Klein")
					pass
			else:
				print("Non-Klein flip")
				# print(ivt)
				print("")
				# print(jtet)
	print("")
	print("Length of Kein Flip list:", len(kein_flip))
	print("")
	print("Length of Self Kein Flip:", len(self_kein))
	print("")
	for isk in self_kein:
		print(isk)

# ******************************************************************************
# Generate all sign permutations of an nxn Identity Matrix
def gen_sign_perm(n):

	n 			= int(n)
	items		= [1] * n
	ptemp = []

	for signs in itertools.product([-1,1], repeat=len(items)):
		temp = np.array([a*sign for a,sign in zip(items,signs)],dtype=int)
		ptemp.append(temp)
	return ptemp

# ******************************************************************************
# Hardcoded bracket-overbar representation of tetrad sets.
def pset_string_format(req_pgroup):

	pq_set 		= {}

	""" {CM} L-matrix tetrad set - P2 Gadget values """
	pq_set['P1']	= [
		"[1,4,2,3], [2,3,-1,-4], [3,-2,4,-1], [4,-1,-3,2]",
		"[1,4,2,3], [2,-3,-1,4], [3,2,-4,-1], [4,-1,3,-2]",
		"[1,-4,2,3], [2,3,-1,4], [3,-2,-4,-1], [4,1,3,-2]",
		"[1,-4,2,3], [2,-3,-1,-4], [3,2,4,-1], [4,1,-3,2]",
		"[1,4,2,-3], [2,3,-1,4], [3,-2,4,1], [4,-1,-3,-2]",
		"[1,4,2,-3], [2,-3,-1,-4], [3,2,-4,1], [4,-1,3,2]",
		"[1,-4,2,-3], [2,3,-1,-4], [3,-2,-4,1], [4,1,3,2]",
		"[1,-4,2,-3], [2,-3,-1,4], [3,2,4,1], [4,1,-3,-2]",
		"[1,4,-2,3], [2,3,1,-4], [3,-2,-4,-1], [4,-1,3,2]",
		"[1,4,-2,3], [2,-3,1,4], [3,2,4,-1], [4,-1,-3,-2]",
		"[1,-4,-2,3], [2,3,1,4], [3,-2,4,-1], [4,1,-3,-2]",
		"[1,-4,-2,3], [2,-3,1,-4], [3,2,-4,-1], [4,1,3,2]",
		"[1,4,-2,-3], [2,3,1,4], [3,-2,-4,1], [4,-1,3,-2]",
		"[1,4,-2,-3], [2,-3,1,-4], [3,2,4,1], [4,-1,-3,2]",
		"[1,-4,-2,-3], [2,3,1,-4], [3,-2,4,1], [4,1,-3,2]",
		"[1,-4,-2,-3], [2,-3,1,4], [3,2,-4,1], [4,1,3,-2]"	# 16 count
          ]

	""" {TM} L-matrix tetrad set - P2 Gadget values """
	pq_set['P2'] 	= [
		"[1,3,4,2], [2,-4,3,-1], [3,-1,-2,4], [4,2,-1,-3]",
		"[1,3,4,2], [2,4,-3,-1], [3,-1,2,-4], [4,-2,-1,3]",
		"[1,3,-4,2], [2,4,3,-1], [3,-1,-2,-4], [4,-2,1,3]",
		"[1,3,-4,2], [2,-4,-3,-1], [3,-1,2,4], [4,2,1,-3]",
		"[1,-3,4,2], [2,4,3,-1], [3,1,-2,4], [4,-2,-1,-3]",
		"[1,-3,4,2], [2,-4,-3,-1], [3,1,2,-4], [4,2,-1,3]",
		"[1,-3,-4,2], [2,-4,3,-1], [3,1,-2,-4], [4,2,1,3]",
		"[1,-3,-4,2], [2,4,-3,-1], [3,1,2,4], [4,-2,1,-3]",
		"[1,3,4,-2], [2,-4,3,1], [3,-1,-2,-4], [4,2,-1,3]",
		"[1,3,4,-2], [2,4,-3,1], [3,-1,2,4], [4,-2,-1,-3]",
		"[1,3,-4,-2], [2,4,3,1], [3,-1,-2,4], [4,-2,1,-3]",
		"[1,3,-4,-2], [2,-4,-3,1], [3,-1,2,-4], [4,2,1,3]",
		"[1,-3,4,-2], [2,4,3,1], [3,1,-2,-4], [4,-2,-1,3]",
		"[1,-3,4,-2], [2,-4,-3,1], [3,1,2,4], [4,2,-1,-3]",
		"[1,-3,-4,-2], [2,-4,3,1], [3,1,-2,4], [4,2,1,-3]",
		"[1,-3,-4,-2], [2,4,-3,1], [3,1,2,-4], [4,-2,1,3]"
		]
	# pq_set['P2'] 	= [
	# 	"[1,3,2,4], [2,-4,-1,3], [3,-1,4,-2], [4,2,-3,-1]",
	# 	"[1,3,2,4], [2,4,-1,-3], [3,-1,-4,2], [4,-2,3,-1]",
	# 	"[1,3,2,-4], [2,4,-1,3], [3,-1,-4,-2], [4,-2,3,1]",
	# 	"[1,3,2,-4], [2,-4,-1,-3], [3,-1,4,2], [4,2,-3,1]",
	# 	"[1,-3,2,4], [2,4,-1,3], [3,1,4,-2], [4,-2,-3,-1]",
	# 	"[1,-3,2,4], [2,-4,-1,-3], [3,1,-4,2], [4,2,3,-1]",
	# 	"[1,-3,2,-4], [2,-4,-1,3], [3,1,-4,-2], [4,2,3,1]",
	# 	"[1,-3,2,-4], [2,4,-1,-3], [3,1,4,2], [4,-2,-3,1]",
	# 	"[1,3,-2,4], [2,-4,1,3], [3,-1,-4,-2], [4,2,3,-1]",
	# 	"[1,3,-2,4], [2,4,1,-3], [3,-1,4,2], [4,-2,-3,-1]",
	# 	"[1,3,-2,-4], [2,4,1,3], [3,-1,4,-2], [4,-2,-3,1]",
	# 	"[1,3,-2,-4], [2,-4,1,-3], [3,-1,-4,2], [4,2,3,1]",
	# 	"[1,-3,-2,4], [2,4,1,3], [3,1,-4,-2], [4,-2,3,-1]",
	# 	"[1,-3,-2,4], [2,-4,1,-3], [3,1,4,2], [4,2,-3,-1]",
	# 	"[1,-3,-2,-4], [2,-4,1,3], [3,1,4,-2], [4,2,-3,1]",
	# 	"[1,-3,-2,-4], [2,4,1,-3], [3,1,-4,2], [4,-2,3,1]"
	# 	]

	""" {VM} L-matrix tetrad set - P3  Gadget values """
	pq_set['P3']	= [
		"[1,3,2,4], [2,-4,-1,3], [3,-1,4,-2], [4,2,-3,-1]",
		"[1,3,2,4], [2,4,-1,-3], [3,-1,-4,2], [4,-2,3,-1]",
		"[1,3,2,-4], [2,4,-1,3], [3,-1,-4,-2], [4,-2,3,1]",
		"[1,3,2,-4], [2,-4,-1,-3], [3,-1,4,2], [4,2,-3,1]",
		"[1,-3,2,4], [2,4,-1,3], [3,1,4,-2], [4,-2,-3,-1]",
		"[1,-3,2,4], [2,-4,-1,-3], [3,1,-4,2], [4,2,3,-1]",
		"[1,-3,2,-4], [2,-4,-1,3], [3,1,-4,-2], [4,2,3,1]",
		"[1,-3,2,-4], [2,4,-1,-3], [3,1,4,2], [4,-2,-3,1]",
		"[1,3,-2,4], [2,-4,1,3], [3,-1,-4,-2], [4,2,3,-1]",
		"[1,3,-2,4], [2,4,1,-3], [3,-1,4,2], [4,-2,-3,-1]",
		"[1,3,-2,-4], [2,4,1,3], [3,-1,4,-2], [4,-2,-3,1]",
		"[1,3,-2,-4], [2,-4,1,-3], [3,-1,-4,2], [4,2,3,1]",
		"[1,-3,-2,4], [2,4,1,3], [3,1,-4,-2], [4,-2,3,-1]",
		"[1,-3,-2,4], [2,-4,1,-3], [3,1,4,2], [4,2,-3,-1]",
		"[1,-3,-2,-4], [2,-4,1,3], [3,1,4,-2], [4,2,-3,1]",
		"[1,-3,-2,-4], [2,4,1,-3], [3,1,-4,2], [4,-2,3,1]"
		]

	""" {VM1} L-matrix tetrad set - P4 Gadget values """
	pq_set['P4']	 = [
		"[1,4,3,2], [2,3,-4,-1], [3,-2,-1,4], [4,-1,2,-3]",
		"[1,4,3,2], [2,-3,4,-1], [3,2,-1,-4], [4,-1,-2,3]",
		"[1,-4,3,2], [2,3,4,-1], [3,-2,-1,-4], [4,1,-2,3]",
		"[1,-4,3,2], [2,-3,-4,-1], [3,2,-1,4], [4,1,2,-3]",
		"[1,4,-3,2], [2,3,4,-1], [3,-2,1,4], [4,-1,-2,-3]",
		"[1,4,-3,2], [2,-3,-4,-1], [3,2,1,-4], [4,-1,2,3]",
		"[1,-4,-3,2], [2,3,-4,-1], [3,-2,1,-4], [4,1,2,3]",
		"[1,-4,-3,2], [2,-3,4,-1], [3,2,1,4], [4,1,-2,-3]",
		"[1,4,3,-2], [2,3,-4,1], [3,-2,-1,-4], [4,-1,2,3]",
		"[1,4,3,-2], [2,-3,4,1], [3,2,-1,4], [4,-1,-2,-3]",
		"[1,-4,3,-2], [2,3,4,1], [3,-2,-1,4], [4,1,-2,-3]",
		"[1,-4,3,-2], [2,-3,-4,1], [3,2,-1,-4], [4,1,2,3]",
		"[1,4,-3,-2], [2,3,4,1], [3,-2,1,-4], [4,-1,-2,3]",
		"[1,4,-3,-2], [2,-3,-4,1], [3,2,1,4], [4,-1,2,-3]",
		"[1,-4,-3,-2], [2,3,-4,1], [3,-2,1,4], [4,1,2,-3]",
		"[1,-4,-3,-2], [2,-3,4,1], [3,2,1,-4], [4,1,-2,3]"
		]

	""" {VM2} L-matrix tetrad set - P5 Gadget values """
	pq_set['P5'] 	= [
		"[1,2,4,3]. [2,-1,3,-4], [3,4,-2,-1], [4,-3,-1,2]",
		"[1,2,4,3], [2,-1,-3,4], [3,-4,2,-1], [4,3,-1,-2]",
		"[1,2,-4,3], [2,-1,3,4], [3,-4,-2,-1], [4,3,1,-2]",
		"[1,2,-4,3], [2,-1,-3,-4], [3,4,2,-1], [4,-3,1,2]",
		"[1,2,4,-3], [2,-1,3,4], [3,4,-2,1], [4,-3,-1,-2]",
		"[1,2,4,-3], [2,-1,-3,-4], [3,-4,2,1], [4,3,-1,2]",
		"[1,2,-4,-3], [2,-1,3,-4], [3,-4,-2,1], [4,3,1,2]",
		"[1,2,-4,-3], [2,-1,-3,4], [3,4,2,1], [4,-3,1,-2]",
		"[1,-2,4,3], [2,1,3,-4], [3,-4,-2,-1], [4,3,-1,2]",
		"[1,-2,4,3], [2,1,-3,4], [3,4,2,-1], [4,-3,-1,-2]",
		"[1,-2,-4,3], [2,1,3,4], [3,4,-2,-1], [4,-3,1,-2]",
		"[1,-2,-4,3], [2,1,-3,-4], [3,-4,2,-1], [4,3,1,2]",
		"[1,-2,4,-3], [2,1,3,4], [3,-4,-2,1], [4,3,-1,-2]",
		"[1,-2,4,-3], [2,1,-3,-4], [3,4,2,1], [4,-3,-1,2]",
		"[1,-2,-4,-3], [2,1,3,-4], [3,4,-2,1], [4,-3,1,2]",
		"[1,-2,-4,-3], [2,1,-3,4], [3,-4,2,1], [4,3,1,-2]"
		]

	""" {VM3} L-matrix tetrad set - P6 Gadget values """
	pq_set['P6']	 = [
		"[1,2,3,4], [2,-1,-4,3], [3,4,-1,-2], [4,-3,2,-1]",
		"[1,2,3,4], [2,-1,4,-3], [3,-4,-1,2], [4,3,-2,-1]",
		"[1,2,3,-4], [2,-1,4,3], [3,-4,-1,-2], [4,3,-2,1]",
		"[1,2,3,-4], [2,-1,-4,-3], [3,4,-1,2], [4,-3,2,1]",
		"[1,2,-3,4], [2,-1,4,3], [3,4,1,-2], [4,-3,-2,-1]",
		"[1,2,-3,4], [2,-1,-4,-3], [3,-4,1,2], [4,3,2,-1]",
		"[1,2,-3,-4], [2,-1,-4,3], [3,-4,1,-2], [4,3,2,1]",
		"[1,2,-3,-4], [2,-1,4,-3], [3,4,1,2], [4,-3,-2,1]",
		"[1,-2,3,4], [2,1,-4,3], [3,-4,-1,-2], [4,3,2,-1]",
		"[1,-2,3,4], [2,1,4,-3], [3,4,-1,2], [4,-3,-2,-1]",
		"[1,-2,3,-4], [2,1,4,3], [3,4,-1,-2], [4,-3,-2,1]",
		"[1,-2,3,-4], [2,1,-4,-3], [3,-4,-1,2], [4,3,2,1]",
		"[1,-2,-3,4], [2,1,4,3], [3,-4,1,-2], [4,3,-2,-1]",
		"[1,-2,-3,4], [2,1,-4,-3], [3,4,1,2], [4,-3,2,-1]",
		"[1,-2,-3,-4], [2,1,-4,3], [3,4,1,-2], [4,-3,2,1]",
		"[1,-2,-3,-4], [2,1,4,-3], [3,-4,1,2], [4,3,-2,1]"
		]
	""" {VM3} L-matrix tetrad set - P6 Gadget values """
	# pq_set['P6']	 = [
	# 	"[3,4,1,2], [4,-3,-2,1], [1,2,-3,-4], [2,-1,4,-3]",
	# 	"[3,4,1,2], [4,-3,2,-1], [1,-2,-3,4], [2,1,-4,-3]",
	# 	"[3,-4,1,2], [4,3,-2,1], [1,-2,-3,-4], [2,1,4,-3]",
	# 	"[3,-4,1,2], [4,3,2,-1], [1,2,-3,4], [2,-1,-4,-3]",
	# 	"[3,4,-1,2], [4,-3,2,1], [1,2,3,-4], [2,-1,-4,-3]",
	# 	"[3,4,-1,2], [4,-3,-2,-1], [1,-2,3,4], [2,1,4,-3]",
	# 	"[3,-4,-1,2], [4,3,2,1], [1,-2,3,-4], [2,1,-4,-3]",
	# 	"[3,-4,-1,2], [4,3,-2,-1], [1,2,3,4], [2,-1,4,-3]",
	# 	"[3,4,1,-2], [4,-3,2,1], [1,-2,-3,-4], [2,1,-4,3]",
	# 	"[3,4,1,-2], [4,-3,-2,-1], [1,2,-3,4], [2,-1,4,3]",
	# 	"[3,-4,1,-2], [4,3,2,1], [1,2,-3,-4], [2,-1,-4,3]",
	# 	"[3,-4,1,-2], [4,3,-2,-1], [1,-2,-3,4], [2,1,4,3]",
	# 	"[3,4,-1,-2], [4,-3,-2,1], [1,-2,3,-4], [2,1,4,3]",
	# 	"[3,4,-1,-2], [4,-3,2,-1], [1,2,3,4], [2,-1,-4,3]",
	# 	"[3,-4,-1,-2], [4,3,-2,1], [1,2,3,-4], [2,-1,4,3]",
	# 	"[3,-4,-1,-2], [4,3,2,-1], [1,-2,3,4], [2,1,-4,3]"
	# 	]

	pq_set['ALL']	 = [
		"[1,4,2,3], [2,3,-1,-4], [3,-2,4,-1], [4,-1,-3,2]",
		"[1,4,2,3], [2,-3,-1,4], [3,2,-4,-1], [4,-1,3,-2]",
		"[1,-4,2,3], [2,3,-1,4], [3,-2,-4,-1], [4,1,3,-2]",
		"[1,-4,2,3], [2,-3,-1,-4], [3,2,4,-1], [4,1,-3,2]",
		"[1,4,2,-3], [2,3,-1,4], [3,-2,4,1], [4,-1,-3,-2]",
		"[1,4,2,-3], [2,-3,-1,-4], [3,2,-4,1], [4,-1,3,2]",
		"[1,-4,2,-3], [2,3,-1,-4], [3,-2,-4,1], [4,1,3,2]",
		"[1,-4,2,-3], [2,-3,-1,4], [3,2,4,1], [4,1,-3,-2]",
		"[1,4,-2,3], [2,3,1,-4], [3,-2,-4,-1], [4,-1,3,2]",
		"[1,4,-2,3], [2,-3,1,4], [3,2,4,-1], [4,-1,-3,-2]",
		"[1,-4,-2,3], [2,3,1,4], [3,-2,4,-1], [4,1,-3,-2]",
		"[1,-4,-2,3], [2,-3,1,-4], [3,2,-4,-1], [4,1,3,2]",
		"[1,4,-2,-3], [2,3,1,4], [3,-2,-4,1], [4,-1,3,-2]",
		"[1,4,-2,-3], [2,-3,1,-4], [3,2,4,1], [4,-1,-3,2]",
		"[1,-4,-2,-3], [2,3,1,-4], [3,-2,4,1], [4,1,-3,2]",
		"[1,-4,-2,-3], [2,-3,1,4], [3,2,-4,1], [4,1,3,-2]", # 16 end of pl1

		"[1,3,2,4], [2,-4,-1,3], [3,-1,4,-2], [4,2,-3,-1]",
		"[1,3,2,4], [2,4,-1,-3], [3,-1,-4,2], [4,-2,3,-1]",
		"[1,3,2,-4], [2,4,-1,3], [3,-1,-4,-2], [4,-2,3,1]",
		"[1,3,2,-4], [2,-4,-1,-3], [3,-1,4,2], [4,2,-3,1]",
		"[1,-3,2,4], [2,4,-1,3], [3,1,4,-2], [4,-2,-3,-1]",
		"[1,-3,2,4], [2,-4,-1,-3], [3,1,-4,2], [4,2,3,-1]",
		"[1,-3,2,-4], [2,-4,-1,3], [3,1,-4,-2], [4,2,3,1]",
		"[1,-3,2,-4], [2,4,-1,-3], [3,1,4,2], [4,-2,-3,1]",
		"[1,3,-2,4], [2,-4,1,3], [3,-1,-4,-2], [4,2,3,-1]",
		"[1,3,-2,4], [2,4,1,-3], [3,-1,4,2], [4,-2,-3,-1]",
		"[1,3,-2,-4], [2,4,1,3], [3,-1,4,-2], [4,-2,-3,1]",
		"[1,3,-2,-4], [2,-4,1,-3], [3,-1,-4,2], [4,2,3,1]",
		"[1,-3,-2,4], [2,4,1,3], [3,1,-4,-2], [4,-2,3,-1]",
		"[1,-3,-2,4], [2,-4,1,-3], [3,1,4,2], [4,2,-3,-1]",
		"[1,-3,-2,-4], [2,-4,1,3], [3,1,4,-2], [4,2,-3,1]",
		"[1,-3,-2,-4], [2,4,1,-3], [3,1,-4,2], [4,-2,3,1]", # end of pl2

		"[1,3,4,2], [2,-4,3,-1], [3,-1,-2,4], [4,2,-1,-3]",
		"[1,3,4,2], [2,4,-3,-1], [3,-1,2,-4], [4,-2,-1,3]",
		"[1,3,-4,2], [2,4,3,-1], [3,-1,-2,-4], [4,-2,1,3]",
		"[1,3,-4,2], [2,-4,-3,-1], [3,-1,2,4], [4,2,1,-3]",
		"[1,-3,4,2], [2,4,3,-1], [3,1,-2,4], [4,-2,-1,-3]",
		"[1,-3,4,2], [2,-4,-3,-1], [3,1,2,-4], [4,2,-1,3]",
		"[1,-3,-4,2], [2,-4,3,-1], [3,1,-2,-4], [4,2,1,3]",
		"[1,-3,-4,2], [2,4,-3,-1], [3,1,2,4], [4,-2,1,-3]",
		"[1,3,4,-2], [2,-4,3,1], [3,-1,-2,-4], [4,2,-1,3]",
		"[1,3,4,-2], [2,4,-3,1], [3,-1,2,4], [4,-2,-1,-3]",
		"[1,3,-4,-2], [2,4,3,1], [3,-1,-2,4], [4,-2,1,-3]",
		"[1,3,-4,-2], [2,-4,-3,1], [3,-1,2,-4], [4,2,1,3]",
		"[1,-3,4,-2], [2,4,3,1], [3,1,-2,-4], [4,-2,-1,3]",
		"[1,-3,4,-2], [2,-4,-3,1], [3,1,2,4], [4,2,-1,-3]",
		"[1,-3,-4,-2], [2,-4,3,1], [3,1,-2,4], [4,2,1,-3]",
		"[1,-3,-4,-2], [2,4,-3,1], [3,1,2,-4], [4,-2,1,3]", # end of pl3

		"[1,4,3,2], [2,3,-4,-1], [3,-2,-1,4], [4,-1,2,-3]",
		"[1,4,3,2], [2,-3,4,-1], [3,2,-1,-4], [4,-1,-2,3]",
		"[1,-4,3,2], [2,3,4,-1], [3,-2,-1,-4], [4,1,-2,3]",
		"[1,-4,3,2], [2,-3,-4,-1], [3,2,-1,4], [4,1,2,-3]",
		"[1,4,-3,2], [2,3,4,-1], [3,-2,1,4], [4,-1,-2,-3]",
		"[1,4,-3,2], [2,-3,-4,-1], [3,2,1,-4], [4,-1,2,3]",
		"[1,-4,-3,2], [2,3,-4,-1], [3,-2,1,-4], [4,1,2,3]",
		"[1,-4,-3,2], [2,-3,4,-1], [3,2,1,4], [4,1,-2,-3]",
		"[1,4,3,-2], [2,3,-4,1], [3,-2,-1,-4], [4,-1,2,3]",
		"[1,4,3,-2], [2,-3,4,1], [3,2,-1,4], [4,-1,-2,-3]",
		"[1,-4,3,-2], [2,3,4,1], [3,-2,-1,4], [4,1,-2,-3]",
		"[1,-4,3,-2], [2,-3,-4,1], [3,2,-1,-4], [4,1,2,3]",
		"[1,4,-3,-2], [2,3,4,1], [3,-2,1,-4], [4,-1,-2,3]",
		"[1,4,-3,-2], [2,-3,-4,1], [3,2,1,4], [4,-1,2,-3]",
		"[1,-4,-3,-2], [2,3,-4,1], [3,-2,1,4], [4,1,2,-3]",
		"[1,-4,-3,-2], [2,-3,4,1], [3,2,1,-4], [4,1,-2,3]", # end of pl4

		"[1,2,4,3]. [2,-1,3,-4], [3,4,-2,-1], [4,-3,-1,2]",
		"[1,2,4,3], [2,-1,-3,4], [3,-4,2,-1], [4,3,-1,-2]",
		"[1,2,-4,3], [2,-1,3,4], [3,-4,-2,-1], [4,3,1,-2]",
		"[1,2,-4,3], [2,-1,-3,-4], [3,4,2,-1], [4,-3,1,2]",
		"[1,2,4,-3], [2,-1,3,4], [3,4,-2,1], [4,-3,-1,-2]",
		"[1,2,4,-3], [2,-1,-3,-4], [3,-4,2,1], [4,3,-1,2]",
		"[1,2,-4,-3], [2,-1,3,-4], [3,-4,-2,1], [4,3,1,2]",
		"[1,2,-4,-3], [2,-1,-3,4], [3,4,2,1], [4,-3,1,-2]",
		"[1,-2,4,3], [2,1,3,-4], [3,-4,-2,-1], [4,3,-1,2]",
		"[1,-2,4,3], [2,1,-3,4], [3,4,2,-1], [4,-3,-1,-2]",
		"[1,-2,-4,3], [2,1,3,4], [3,4,-2,-1], [4,-3,1,-2]",
		"[1,-2,-4,3], [2,1,-3,-4], [3,-4,2,-1], [4,3,1,2]",
		"[1,-2,4,-3], [2,1,3,4], [3,-4,-2,1], [4,3,-1,-2]",
		"[1,-2,4,-3], [2,1,-3,-4], [3,4,2,1], [4,-3,-1,2]",
		"[1,-2,-4,-3], [2,1,3,-4], [3,4,-2,1], [4,-3,1,2]",
		"[1,-2,-4,-3], [2,1,-3,4], [3,-4,2,1], [4,3,1,-2]", # end of pl5

		"[1,2,3,4], [2,-1,-4,3], [3,4,-1,-2], [4,-3,2,-1]",
		"[1,2,3,4], [2,-1,4,-3], [3,-4,-1,2], [4,3,-2,-1]",
		"[1,2,3,-4], [2,-1,4,3], [3,-4,-1,-2], [4,3,-2,1]",
		"[1,2,3,-4], [2,-1,-4,-3], [3,4,-1,2], [4,-3,2,1]",
		"[1,2,-3,4], [2,-1,4,3], [3,4,1,-2], [4,-3,-2,-1]",
		"[1,2,-3,4], [2,-1,-4,-3], [3,-4,1,2], [4,3,2,-1]",
		"[1,2,-3,-4], [2,-1,-4,3], [3,-4,1,-2], [4,3,2,1]",
		"[1,2,-3,-4], [2,-1,4,-3], [3,4,1,2], [4,-3,-2,1]",
		"[1,-2,3,4], [2,1,-4,3], [3,-4,-1,-2], [4,3,2,-1]",
		"[1,-2,3,4], [2,1,4,-3], [3,4,-1,2], [4,-3,-2,-1]",
		"[1,-2,3,-4], [2,1,4,3], [3,4,-1,-2], [4,-3,-2,1]",
		"[1,-2,3,-4], [2,1,-4,-3], [3,-4,-1,2], [4,3,2,1]",
		"[1,-2,-3,4], [2,1,4,3], [3,-4,1,-2], [4,3,-2,-1]",
		"[1,-2,-3,4], [2,1,-4,-3], [3,4,1,2], [4,-3,2,-1]",
		"[1,-2,-3,-4], [2,1,-4,3], [3,4,1,-2], [4,-3,2,1]",
		"[1,-2,-3,-4], [2,1,4,-3], [3,-4,1,2], [4,3,-2,1]"
		]

	pl13 = [
		"[1,4,2,3], [2,3,-1,-4], [3,-2,4,-1], [4,-1,-3,2]",
		"[1,4,2,3], [2,-3,-1,4], [3,2,-4,-1], [4,-1,3,-2]",
		"[1,-4,2,3], [2,3,-1,4], [3,-2,-4,-1], [4,1,3,-2]",
		"[1,-4,2,3], [2,-3,-1,-4], [3,2,4,-1], [4,1,-3,2]",
		"[1,4,2,-3], [2,3,-1,4], [3,-2,4,1], [4,-1,-3,-2]",
		"[1,4,2,-3], [2,-3,-1,-4], [3,2,-4,1], [4,-1,3,2]",
		"[1,-4,2,-3], [2,3,-1,-4], [3,-2,-4,1], [4,1,3,2]",
		"[1,-4,2,-3], [2,-3,-1,4], [3,2,4,1], [4,1,-3,-2]",
		"[1,4,-2,3], [2,3,1,-4], [3,-2,-4,-1], [4,-1,3,2]",
		"[1,4,-2,3], [2,-3,1,4], [3,2,4,-1], [4,-1,-3,-2]",
		"[1,-4,-2,3], [2,3,1,4], [3,-2,4,-1], [4,1,-3,-2]",
		"[1,-4,-2,3], [2,-3,1,-4], [3,2,-4,-1], [4,1,3,2]",
		"[1,4,-2,-3], [2,3,1,4], [3,-2,-4,1], [4,-1,3,-2]",
		"[1,4,-2,-3], [2,-3,1,-4], [3,2,4,1], [4,-1,-3,2]",
		"[1,-4,-2,-3], [2,3,1,-4], [3,-2,4,1], [4,1,-3,2]",
		"[1,-4,-2,-3], [2,-3,1,4], [3,2,-4,1], [4,1,3,-2]",
		"[1,3,4,2], [2,-4,3,-1], [3,-1,-2,4], [4,2,-1,-3]",
		"[1,3,4,2], [2,4,-3,-1], [3,-1,2,-4], [4,-2,-1,3]",
		"[1,3,-4,2], [2,4,3,-1], [3,-1,-2,-4], [4,-2,1,3]",
		"[1,3,-4,2], [2,-4,-3,-1], [3,-1,2,4], [4,2,1,-3]",
		"[1,-3,4,2], [2,4,3,-1], [3,1,-2,4], [4,-2,-1,-3]",
		"[1,-3,4,2], [2,-4,-3,-1], [3,1,2,-4], [4,2,-1,3]",
		"[1,-3,-4,2], [2,-4,3,-1], [3,1,-2,-4], [4,2,1,3]",
		"[1,-3,-4,2], [2,4,-3,-1], [3,1,2,4], [4,-2,1,-3]",
		"[1,3,4,-2], [2,-4,3,1], [3,-1,-2,-4], [4,2,-1,3]",
		"[1,3,4,-2], [2,4,-3,1], [3,-1,2,4], [4,-2,-1,-3]",
		"[1,3,-4,-2], [2,4,3,1], [3,-1,-2,4], [4,-2,1,-3]",
		"[1,3,-4,-2], [2,-4,-3,1], [3,-1,2,-4], [4,2,1,3]",
		"[1,-3,4,-2], [2,4,3,1], [3,1,-2,-4], [4,-2,-1,3]",
		"[1,-3,4,-2], [2,-4,-3,1], [3,1,2,4], [4,2,-1,-3]",
		"[1,-3,-4,-2], [2,-4,3,1], [3,1,-2,4], [4,2,1,-3]",
		"[1,-3,-4,-2], [2,4,-3,1], [3,1,2,-4], [4,-2,1,3]"
		]


	if req_pgroup in pq_set:
		return pq_set[req_pgroup]
	# if req_pgroup == 'pl1':
	# 	return pl1
	# elif req_pgroup == 'pl2':
	# 	return pl2
	# elif req_pgroup == 'pl3':
	# 	return pl3
	# elif req_pgroup == 'pl4':
	# 	return pl4
	# elif req_pgroup == "pl5":
	# 	return pl5
	# elif req_pgroup == "pl6":
	# 	return pl6
	# elif req_pgroup == "plall":
	# 	return plall
	# elif req_pgroup == "pl13":
	# 	return pl13

# **************************************************************************
# Execute main()
if __name__ == "__main__":
	start_time = time.time()

	main()
	print("-- Execution time --")
	print("---- %s seconds ----" % (time.time() - start_time))
