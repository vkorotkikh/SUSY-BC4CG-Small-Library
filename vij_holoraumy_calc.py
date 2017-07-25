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
from pprint import pprint

pr_sw	= 0

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


# ******************************************************************************
# Do the final Vij calculation
def calc_holoraumy_mats(main_tetrad_list, pset_arg):

	""" Remember that the main_tetrad_ark is a list of lists,
		with each list containing four tuples, with tuples being
		matrix number and the matrices itself. """

	holo_mats	= []
	r_matrices	= []

	for ti, teti in enumerate(main_tetrad_list):
		if pr_sw:
			print("# ********************************")
			print("								     ")
			print("Tetrad i: ", ti)
			# calculate_vijmatset(teti)
		# fermionic_holomats(teti)
		holomat, rmat = bosonic_holomats(teti)
		holo_mats.append(holomat)
		r_matrices.append(rmat)

	nicely_print(holo_mats, r_matrices, pset_arg)


# ******************************************************************************
# Calculating Fermionic holoraumy matrices for given Adinkra
def fermionic_holomats(adinkra):

	# vij_possibilities 	= alphas_betas()

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
		""" Probably needs 1/2	"""
		holo_mat	= np.dot(rimat, ljmat) - np.dot(rjmat, limat)
		# temp_mat	= np.dot(tr_limat, ljmat) - np.dot(tr_ljmat, limat)
		vij_fermi.append(holomat)

	return vij_fermi, r_matrices

# ******************************************************************************
# Calculating Bosonic holoraumy matrices for given Adinkra
def bosonic_holomats(adinkra):

	""" Store n Vij bosonic matrices in vij_bosonic	"""
	vij_bosonic	= []
	r_matrices	= [np.asarray(np.transpose(mat)) for mat in adinkra]

	ij_indices = list(itertools.combinations([0,1,2,3], 2))

	for ijtup in ij_indices:

		limat	= adinkra[ijtup[0]]
		ljmat	= adinkra[ijtup[1]]
		rimat	= np.transpose(limat)
		rjmat	= np.transpose(ljmat)

		""" Vij bosonic eq	1501.00101 3.5	"""
		holo_mat = np.dot(limat, rjmat) - np.dot(ljmat, rimat)
		vij_bosonic.append(holo_mat)

	return vij_bosonic, r_matrices


# ******************************************************************************
# Calculating Bosonic holoraumy matrices for given Adinkra
def nicely_print(holo_mats, rmats, pset_arg):

	holos	= []
	holos = [np.asarray(x) for x in holo_mats]

	print("# ********************************")
	print("Calculated matrices for: ", pset_arg)
	# print("Bosonic or Fermionic....")
	print("")
	lenh, lenr = len(holos), len(rmats)
	# print("Length holo_mats: ", len(holo_mats))
	print("Length holo_mats: ", lenh)
	# print("Length r_matrices: ", len(r_matrices))
	print("Length r_matrices: ", lenr)
	print("")

	np.set_printoptions(precision=2, suppress=True, linewidth=100)
	ij_ind	= list(itertools.combinations([0,1,2,3], 2))
	setlen = 0
	if len(holos) == len(rmats):
		setlen = len(holos)
	else:
		print("LENGTH MISMATCH ERROR")

	for zi in range(0, lenh):
		temph = holos[zi]
		tempr = rmats[zi]
		print("#********************************")
		print("Adinkra #", zi, "Bosonic Holoraumy Matrices")
		vij_strings	= []
		for ijtup in ij_ind:
			ij_temp		= str(ijtup[0] + 1) + str(ijtup[1] + 1)
			ijstr		= "~V_{" + ij_temp + "}"
			vij_strings.append(ijstr)
		v13strings = " \t" + vij_strings[0] + " \t \t \t" + vij_strings[1] + \
		" \t \t \t" + vij_strings[2]
		v46strings = " \t" + vij_strings[3] + " \t \t \t" + vij_strings[4] + \
		" \t \t \t" + vij_strings[5]

		""" Convoluted way of printing out numpy matrices	"""
		mat13 		= temph[0:3]
		mat13str 	= [np.array_str(y)[1:-1] for y in mat13]
		tm13		= []
		for matstr in mat13str:
			# onemat = [ix.lstrip() for ix in matstr.split('\n')]
			tm13.append([ix.lstrip() for ix in matstr.split('\n')])
		print(v13strings)
		for ix in range(0,4):
			pstr = tm13[0][ix] + "\t\t" + tm13[1][ix] + "\t\t" + tm13[2][ix]
			print(pstr)

		mat46		= temph[3:6]
		mat46str	= [np.array_str(y)[1:-1] for y in mat46]
		tm46		= []
		for matstr in mat46str:
			tm46.append([ix.lstrip() for ix in matstr.split('\n')])
		print(v46strings)
		for ix in range(0,4):
			pstr = tm46[0][ix] + "\t\t" + tm46[1][ix] + "\t\t" + tm46[2][ix]
			print(pstr)
		print("")

		# tempr = rmats[zi]
		print("Adinkra #", zi, "R matrices")
		for rl in [ tempr[:2], tempr[2:]]:
			rltostr = [np.array_str(y)[1:-1] for y in rl]
			rtm	= []
			for matstr in rltostr:
				rtm.append([ix.lstrip() for ix in matstr.split('\n')])
			for ix in range(0,4):
				pstr = rtm[0][ix] + "\t\t" + rtm[1][ix]
				print(pstr)
			print("")
""" This needs work. Probably later	"""
# def gadgetizing(holomats):
# 		# """ Compare against the 6 possible matrix solutions """
#
# 		tf_bool = 0
# 		for xi, ijx in enumerate(vij_possibilities):
# 			ijx_neg = np.multiply(ijx, -1)
# 			# print(xi)
# 			if np.array_equal(temp_mat, ijx):
# 				tf_bool = 1
# 				if debug:
# 					print("*************$$$$$$$$$$$$$$$$$$ ")
# 					print("l-solution found:")
# 					print(ijx)
# 				tmint = np.int(1)
# 				if xi < 3:
# 					tmp_str = "alpha^" + str((xi + 1))
# 					# print(tmp_str)
# 					# vij_tempset.append([tmp_str, ijstr, tmint])
# 					res_str	= ijstr + " = " + "1 * " + tmp_str
# 					# vij_tempset.append([tmint, ijstr, tmp_str])
# 					vij_tempset.append(res_str)
# 				elif xi >= 3:
# 					tmp_str = "beta^" + str((xi - 2))
# 					res_str	= ijstr + " = " + "1 * " + tmp_str
# 					# vij_tempset.append([tmint, ijstr, tmp_str])
# 					vij_tempset.append(res_str)
# 			elif np.array_equal(temp_mat, ijx_neg):
# 				tf_bool = 1
# 				if debug:
# 					print("*************$$$$$$$$$$$$$$$$$$ ")
# 					print("l-solution found:")
# 					print(ijx_neg)
# 				# xint = (xi + 1) * ( -1)
# 				tmint = np.int(-1)
# 				if xi < 3:
# 					tmp_str = "alpha^" + str((xi + 1))
# 					res_str	= ijstr + " = " + "-1 * " + tmp_str
# 					# print(tmp_str)
# 					# vij_tempset.append([tmint, ijstr, tmp_str])
# 					vij_tempset.append(res_str)
# 				elif xi >= 3:
# 					tmp_str = "beta^" + str((xi - 2))
# 					res_str	= ijstr + " = " + "-1 * " + tmp_str
# 					# vij_tempset.append([tmint, ijstr, tmp_str])
# 					vij_tempset.append(res_str)
# 			else:
# 				if tf_bool == 0 and xi >= 5:
# 					if not(np.array_equal(temp_mat, ijx)) or not np.array_equal(temp_mat, ijx_neg):
# 						print("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx ")
# 						print("Anomaly found:",ijstr)
# 						print(temp_mat)
# 						anomaly_switch = 1
# 		tf_bool = 0
# 	return vij_tempset

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
