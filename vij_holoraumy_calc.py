# ******************************************************************************
# Name:    Calculate Vij matrices and the elle & tilde~elle Coefficients
# Author:  Vadim Korotkikh
# Email:   va.korotki@gmail.com
# Date:    November 2016
# Version: 1.3
#
# Description: Calculates the corresponding Vij Holoraumy matrices for Adinkras
# Using the Holoraumy matrices, finds the corresponding elle or tilde`elle
# coefficients for ALpha or Beta matrices
#
# ******************************************************************************


# ******************************************************************************
# Library Imports

import sys
import math
import time
import itertools
import numpy as np
from numpy import array
from numpy.linalg import inv

# ******************************************************************************
# Main() function.
def main():
	# gen_signm(4)
	pass

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


vij_possibilities 	= alphas_betas()

# ******************************************************************************
# Do the final Vij calculation
def calculate_vij_matrices(main_tetrad_list):

	""" Remember that the main_tetrad_ark is a list of lists,
		with each list containing four tuples, with tuples being
		matrix number and the matrices itself. """


	# vij_possibilities 	= alphas_betas()

	anomaly_switch  = 0
	debug			= 0

	for ti, teti in enumerate(main_tetrad_list):
		if debug:
			print("# ********************************")
			print("								     ")
			print("Tetrad i: ", ti)
			# calculate_vijmatset(teti)
			fermionic_holomats(teti)

# ******************************************************************************
# Calculating Fermionic holoraumy matrices for given Adinkra
def fermionic_holomats(adinkra):

	# vij_possibilities 	= alphas_betas()
	debug			= 0

	""" Store Fermionic Holoraumy matrices """
	vij_fermi	= []
	# r_matrices	= []
	r_matrices	= [np.transpose(mat) for mat in adinkra]
	""" Store 6 Vij matrices in temp_vijmat"""
	# temp_vijmat		= []
	ij_indices	= list(itertools.combinations([0,1,2,3], 2))

	for ijtup in ij_indices:
		limat 		= adinkra[ijtup[0]]
		ljmat 		= adinkra[ijtup[1]]
		ij_temp		= str(ijtup[0] + 1) + str(ijtup[1] + 1)
		""" Enhance appearance of ijstr - V_{ ij } """
		ijstr		= "~V_{" + ij_temp + "}"
		# tr_limat	= np.transpose(limat)
		# tr_ljmat	= np.transpose(ljmat)
		rimat		= np.transpose(limat)
		rjmat		= np.transpose(ljmat)
		""" Vij eq from 1601.00 (3.2) """
		holo_mat	= np.dot(rimat, ljmat) - np.dot(rjmat, limat)
		# temp_mat	= np.dot(tr_limat, ljmat) - np.dot(tr_ljmat, limat)
		vij_fermi.append(holomat)

	return vij_fermi, r_matrices

""" This needs work. Probably later	"""
def gadgetizing(holomats):
		# """ Compare against the 6 possible matrix solutions """

		tf_bool = 0
		for xi, ijx in enumerate(vij_possibilities):
			ijx_neg = np.multiply(ijx, -1)
			# print(xi)
			if np.array_equal(temp_mat, ijx):
				tf_bool = 1
				if debug:
					print("*************$$$$$$$$$$$$$$$$$$ ")
					print("l-solution found:")
					print(ijx)
				tmint = np.int(1)
				if xi < 3:
					tmp_str = "alpha^" + str((xi + 1))
					# print(tmp_str)
					# vij_tempset.append([tmp_str, ijstr, tmint])
					res_str	= ijstr + " = " + "1 * " + tmp_str
					# vij_tempset.append([tmint, ijstr, tmp_str])
					vij_tempset.append(res_str)
				elif xi >= 3:
					tmp_str = "beta^" + str((xi - 2))
					res_str	= ijstr + " = " + "1 * " + tmp_str
					# vij_tempset.append([tmint, ijstr, tmp_str])
					vij_tempset.append(res_str)
			elif np.array_equal(temp_mat, ijx_neg):
				tf_bool = 1
				if debug:
					print("*************$$$$$$$$$$$$$$$$$$ ")
					print("l-solution found:")
					print(ijx_neg)
				# xint = (xi + 1) * ( -1)
				tmint = np.int(-1)
				if xi < 3:
					tmp_str = "alpha^" + str((xi + 1))
					res_str	= ijstr + " = " + "-1 * " + tmp_str
					# print(tmp_str)
					# vij_tempset.append([tmint, ijstr, tmp_str])
					vij_tempset.append(res_str)
				elif xi >= 3:
					tmp_str = "beta^" + str((xi - 2))
					res_str	= ijstr + " = " + "-1 * " + tmp_str
					# vij_tempset.append([tmint, ijstr, tmp_str])
					vij_tempset.append(res_str)
			else:
				if tf_bool == 0 and xi >= 5:
					if not(np.array_equal(temp_mat, ijx)) or not np.array_equal(temp_mat, ijx_neg):
						print("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx ")
						print("Anomaly found:",ijstr)
						print(temp_mat)
						anomaly_switch = 1
		tf_bool = 0

	return vij_tempset

	# print("*************$$$$$$$$$$$$$$$$$$ ")
	# print("Vij Matrix Coefficients Results:")
	# print("")
	# for mvals in calc_check:
	# 	if any(x for x in mvals if x[0].startswith('alpha')) and any(x for x in mvals if x[0].startswith('beta')):
	# 		print("MIXED ALPHA_BETA ERROR")
	# 		print(mvals)
	# 	else:
	# 		print(mvals)
	#
	# print("Length Vij alphas adinkras:", len(vij_alphas))
	# print("Length Vij beta adikras:", len(vij_betas))

# ******************************************************************************
# Creating an abstract Adinkra handler
def calculate_vijmatset_nicely(one_adinkra):

	vij_possibilities 	= alphas_betas()

	debug			= 0

	vij_detail	=	[]
	vij_tempset = 	[]
	""" Store 6 Vij matrices in temp_vijmat"""
	# temp_vijmat		= []
	ij_indices			= list(itertools.combinations([0,1,2,3], 2))

	for ijtup in ij_indices:
		limat 		= one_adinkra[ijtup[0]]
		ljmat 		= one_adinkra[ijtup[1]]
		ij_temp		= str(ijtup[0] + 1) + str(ijtup[1] + 1)
		""" Enhance appearance of ijstr - V_{ ij } """
		ijstr		= "~V_{" + ij_temp + "}"
		tr_limat	= np.transpose(limat)
		tr_ljmat	= np.transpose(ljmat)
		""" Vij eq from 1601.00 (3.2) """
		temp_mat	= np.dot(tr_limat, ljmat) - np.dot(tr_ljmat, limat)
		""" Compare against the 6 possible matrix solutions """
		tf_bool = 0

		for xi, ijx in enumerate(vij_possibilities):
			ijx_neg = np.multiply(ijx, -1)
			# print(xi)
			if np.array_equal(temp_mat, ijx):
				tf_bool = 1
				if debug:
					print("*************$$$$$$$$$$$$$$$$$$ ")
					print("l-solution found:")
					print(ijx)
				tmint = np.int(1)
				if xi < 3:
					tmp_str = "alpha^" + str((xi + 1))
					# print(tmp_str)
					# vij_tempset.append([tmp_str, ijstr, tmint])
					res_str	= ijstr + " = " + "1 * " + tmp_str
					# vij_tempset.append([tmint, ijstr, tmp_str])
					vij_detail.append([ij_temp, tmint, tmp_str])
					vij_tempset.append(res_str)
				elif xi >= 3:
					tmp_str = "beta^" + str((xi - 2))
					res_str	= ijstr + " = " + "1 * " + tmp_str
					# vij_tempset.append([tmint, ijstr, tmp_str])
					vij_detail.append([ij_temp, tmint, tmp_str])
					vij_tempset.append(res_str)
			elif np.array_equal(temp_mat, ijx_neg):
				tf_bool = 1
				if debug:
					print("*************$$$$$$$$$$$$$$$$$$ ")
					print("l-solution found:")
					print(ijx_neg)
				# xint = (xi + 1) * ( -1)
				tmint = np.int(-1)
				if xi < 3:
					tmp_str = "alpha^" + str((xi + 1))
					res_str	= ijstr + " = " + "-1 * " + tmp_str
					# print(tmp_str)
					# vij_tempset.append([tmint, ijstr, tmp_str])
					vij_detail.append([ij_temp, tmint, tmp_str])
					vij_tempset.append(res_str)
				elif xi >= 3:
					tmp_str = "beta^" + str((xi - 2))
					res_str	= ijstr + " = " + "-1 * " + tmp_str
					# vij_tempset.append([tmint, ijstr, tmp_str])
					vij_detail.append([ij_temp, tmint, tmp_str])
					vij_tempset.append(res_str)
			else:
				if tf_bool == 0 and xi >= 5:
					if not(np.array_equal(temp_mat, ijx)) or not np.array_equal(temp_mat, ijx_neg):
						print("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx ")
						print("Anomaly found:",ijstr)
						print(temp_mat)
						anomaly_switch = 1
		tf_bool = 0

	repackaging	=	[(vij_detail[0], vij_detail[5]), (vij_detail[1], vij_detail[4]), (vij_detail[2], vij_detail[3])]
	end_pack	=	[ "", "", ""]

	for itup, vijtup in enumerate(repackaging):
		p1 		= vijtup[0]
		p2 		= vijtup[1]
		ind_pos	= int(p1[2][-1:]) - 1
		if itup == 0:
			if p1[1] == p2[1]:
				if p1[1] < 0:
					end_pack[ind_pos] = "alpha^2"
				elif p1[1] > 0:
					end_pack[ind_pos] = "-alpha^2"
			elif p1[1] != p2[1]:
				if p1[1] > 0:
					end_pack[ind_pos] = "-beta^3"
				elif p1[1] < 0:
					end_pack[ind_pos] = "beta^3"
		if itup == 1:
			if p1[1] != p2[1]:
				if p1[1] < 0 and p2[1] > 0:
					end_pack[ind_pos] = "alpha^3"
				elif p1[1] > 0 and p2[1] < 0:
					end_pack[ind_pos] = "-alpha^3"
			elif p1[1] == p2[1]:
				if p1[1] < 0 and p2[1] < 0:
					end_pack[ind_pos] = "beta^2"
				elif p1[1] > 0 and p2[1] > 0:
					end_pack[ind_pos] = "-beta^2"
		if itup == 2:
			if p1[1] == p2[1]:
				if p1[1] < 0 and p2[1] < 0:
					end_pack[ind_pos] = "alpha^1"
				elif p1[1] > 0 and p2[1] > 0:
					end_pack[ind_pos] = "-alpha^1"
			if p1[1] != p2[1]:
				if p1[1] < 0 and p2[1] > 0:
					end_pack[ind_pos] = "beta^1"
				else:
					end_pack[ind_pos] = "-beta^1"

	return vij_tempset, end_pack


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
