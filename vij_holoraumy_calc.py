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
# Calc. Bosonic or Fermionic Holoraumy mats for Adinkra.
def calc_holoraumy_mats(main_tetrad_list, pset_arg, holotype):
	""" Remember that the main_tetrad_ark is a list of lists,
		with each list four numpy matrices, or a 4 L matrix Adinkra
		> This is hardcoded for a 4 Matrix Adinkra
	"""
	lholotype	= holotype.lower()
	holo_mats	= []
	r_matrices	= []
	adink_def	= []

	if lholotype.startswith('boson'):
		for ti, teti in enumerate(main_tetrad_list):
			if pr_sw:
				print("# ********************************")
				print("								     ")
				print("Tetrad i: ", ti)
				# calculate_vijmatset(teti)
			# fermionic_holomats(teti)
			if len(teti) > 1 and isinstance(teti, tuple) is True:
				print("YES", teti, "", len(teti))
				if isinstance(teti[1], tuple) is True:
					print("TUPLES!!!")
					holomat, rmat = bosonic_holomats(teti[0])
					holo_mats.append(holomat)
					r_matrices.append(rmat)
					adink_def.append(teti[1])
			else:
				print("NO", teti, "", len(teti))
				holomat, rmat = bosonic_holomats(teti)
				holo_mats.append(holomat)
				r_matrices.append(rmat)

	elif lholotype.startswith('fermi'):
		for ti, teti in enumerate(main_tetrad_list):
			if pr_sw:
				print("# ********************************")
				print("								     ")
				print("Tetrad i: ", ti)
			if len(teti) > 1 and isinstance(teti, tuple) is True:
				# print("YES, Tup!", teti, "", len(teti))
				if isinstance(teti[1], tuple) is True:
					print("Tuple inside too")
					holomat, rmat = fermionic_holomats(teti[0])
					holo_mats.append(holomat)
					r_matrices.append(rmat)
					adink_def.append(teti[1])
			else:
				print("Not Tup", teti, "", len(teti))
				holomat, rmat = fermionic_holomats(teti)
				holo_mats.append(holomat)
				r_matrices.append(rmat)

	if lholotype.startswith('boson'):
		# nicely_print_boson(holo_mats, r_matrices, pset_arg, adink_def)
		full_nprint_boson(main_tetrad_list, holo_mats, r_matrices, pset_arg, adink_def)
	elif lholotype.startswith('fermi'):
		nicely_print_fermi(holo_mats, r_matrices, pset_arg, adink_def)


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
		rimat		= np.transpose(limat)
		rjmat		= np.transpose(ljmat)
		""" Vij eq from 1601.00 (3.2) """
		""" Probably needs 1/2	"""
		holo_mat	= np.divide((np.dot(rimat, ljmat) - np.dot(rjmat, limat)),2)
		vij_fermi.append(holo_mat)

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
		holo_mat = np.divide((np.dot(limat, rjmat) - np.dot(ljmat, rimat)),2)
		# holo_mat = np.dot(limat, rjmat) - np.dot(ljmat, rimat)
		vij_bosonic.append(holo_mat)

	return vij_bosonic, r_matrices

# ******************************************************************************
def full_nprint_boson(lmat_list, holo_mats, rmats, pset_arg, adink_def):
	"""
		lmat_list - List of lists w/ each cont. Adinkra 4 L matrices
		holo_mats - List of lists w/ each containg 6 Bosnic Holoraumy matrices
		rmats 	  -	List w/ lists, each containing 4 R matrices
		pset_arg  - String specifying P slices of library, ie P1 or P6
		adink_def - List contains tuples of boolean factors, matrix slice
	"""

	# T/F whether to add boolean factor information to the output text file
	adinkdef_yn	= 0
	if len(adink_def) > 0:
		adinkdef_yn = 1
	holos = [np.asarray(x) for x in holo_mats]
	text_list	= []

	# print("# ********************************")
	text_list.append("#********************************")
	text_list.append("BC4 Coxeter Group Small Library")
	if adinkdef_yn:
		text_list.append(pset_arg + " set: " + adink_def[0][1])
	else:
		text_list.append(pset_arg + " set")
	text_list.append("Calculated Bosonic Holoraumy matrices and R matrices")
	text_list.append("")
	lenh, lenr = len(holos), len(rmats)
	print("Length holo_mats: ", lenh)
	print("Length r_matrices: ", lenr)
	print("")

	np.set_printoptions(precision=2, suppress=True, linewidth=100)
	ij_ind	= list(itertools.combinations([0,1,2,3], 2))
	setlen = 0
	# if len(holos) == len(rmats):
	if lenh == lenr:
		setlen = lenh
	else:
		print("LENGTH MISMATCH ERROR")

	for zi in range(0, lenh):
		temph = holos[zi]
		tempr = rmats[zi]
		templ = lmat_list[zi][0]
		print(templ)
		text_list.append("#********************************")
		text_list.append("Adinkra # " + str(zi))
		if adinkdef_yn:
			print("Def ", adink_def[zi], " Bool Fct", adink_def[zi][0])
			boolstr = (",").join(['(' + str(ix) + ')' for ix in adink_def[zi][0]])
			boolstr = "{" + boolstr + "}"
			print("Bool String:", boolstr)
			text_list.append("Boolean Factor: " + boolstr)
			# text_list.append("Boolean Factor:" + adink_def[zi][0] + " P-set " + adink_def[zi][1])
		text_list.append("Bosonic Holoraumy Matrices")
		# print("Adinkra #", zi)
		# print("Bosonic Holoraumy Matrices")
		vij_strings	= []
		for ijtup in ij_ind:
			ij_temp		= str(ijtup[0] + 1) + str(ijtup[1] + 1)
			ijstr		= "V_{" + ij_temp + "}"
			vij_strings.append(ijstr)
		""" List comprehension for the above 3line forloop, 4line if counting
			vij_strings list declaration	"""
		'''
		vij_strings = [ "V_{" + (str(ijt[0]+1) + str(ijt[1]+1)) + "}" for ijt in ij_ind]
		'''
		# v13strings = " \t" + vij_strings[0] + " \t \t \t" + vij_strings[1] + \
		# " \t \t \t" + vij_strings[2]
		v13strings = "\t" + vij_strings[0] + "\t\t   " + vij_strings[1] + \
		"\t\t   " + vij_strings[2]
		# v46strings = " \t" + vij_strings[3] + " \t \t \t" + vij_strings[4] + \
		# " \t \t \t" + vij_strings[5]
		v46strings = "\t" + vij_strings[3] + "\t\t   " + vij_strings[4] + \
		"\t\t   " + vij_strings[5]

		""" Convoluted way of printing out numpy matrices	"""
		mat13 		= temph[0:3]
		mat13str 	= [np.array_str(y)[1:-1] for y in mat13]
		tm13		= []
		for matstr in mat13str:
			tm13.append([ix.lstrip() for ix in matstr.split('\n')])
		# print(v13strings)
		text_list.append(v13strings)
		for ix in range(0,4):
			pstr = tm13[0][ix] + " \t" + tm13[1][ix] + " \t" + tm13[2][ix]
			# print(pstr)
			text_list.append(pstr)

		mat46		= temph[3:6]
		mat46str	= [np.array_str(y)[1:-1] for y in mat46]
		tm46		= []
		for matstr in mat46str:
			tm46.append([ix.lstrip() for ix in matstr.split('\n')])
		# print(v46strings)
		text_list.append(v46strings)
		for ix in range(0,4):
			pstr = tm46[0][ix] + " \t" + tm46[1][ix] + " \t" + tm46[2][ix]
			# print(pstr)
			text_list.append(pstr)
		text_list.append("")
		# print("")
		""" Printing L Matrices """
		text_list.append("L matrices")
		for ind, lmi in enumerate([templ[:2], templ[2:]]):
			lm	= [np.asarray(x) for x in lmi]
			lmtostr = [np.array_str(y)[1:-1] for y in lm]
			ltm	= []
			for matstr in lmtostr:
				ltm.append([ix.lstrip() for ix in matstr.split('\n')])
			if ind == 0:
				len1, len2 	= len(ltm[0][0]), len(ltm[1][0])
				lbl_str	= ""
				if len2 > len1:
					lbl_str	= (len1//3)*" " + "L1" + len2*" " + "L2"
				elif len1 > len2:
					lbl_str	= (len1//2)*" " + "L1" + len2*" " + " L2"
				elif len1 == len2:
					lbl_str = (len1//2)*" " + "L1" + len2*" " + "L2"
				# print(lbl_str)
				text_list.append(lbl_str)
			elif ind == 1:
				len1, len2 	= len(ltm[0][0]), len(ltm[1][0])
				lbl_str		= "  L3" + len2*" " + " L4"
				# print(lbl_str)
				text_list.append(lbl_str)
			for ix in range(0,4):
				pstr = ltm[0][ix] + "\t" + ltm[1][ix]
				text_list.append(pstr)
		text_list.append("")

		"""	Printing R matrices """
		text_list.append("R matrices")
		for ind, rl in enumerate([ tempr[:2], tempr[2:]]):
			rltostr = [np.array_str(y)[1:-1] for y in rl]
			rtm	= []
			for matstr in rltostr:
				rtm.append([ix.lstrip() for ix in matstr.split('\n')])

			if ind == 0:
				# print("Length: ", len(rtm[0
				len1, len2 	= len(rtm[0][0]), len(rtm[1][0])
				lbl_str	= ""
				if len2 > len1:
					lbl_str	= (len1//3)*" " + "R1" + len2*" " + "R2"
				elif len1 > len2:
					lbl_str	= (len1//2)*" " + "R1" + len2*" " + " R2"
				elif len1 == len2:
					lbl_str = (len1//2)*" " + "R1" + len2*" " + "R2"
				# print(lbl_str)
				text_list.append(lbl_str)
			elif ind == 1:
				len1, len2 	= len(rtm[0][0]), len(rtm[1][0])
				# print(len1, len2)
				lbl_str		= "  R3" + len2*" " + " R4"
				# print(lbl_str)
				text_list.append(lbl_str)
			for ix in range(0,4):
				pstr = rtm[0][ix] + "\t" + rtm[1][ix]
				# print(pstr)
				text_list.append(pstr)
		text_list.append("")

	pmatsfile = "BC4-CoxeterGroup " + pset_arg + " withLmats.txt"
	# pmatsfile = "BC4-CG " + pset_arg + "-BosonicHolos.txt"
	with open(pmatsfile, "w") as wfile:
		for item in text_list:
			wfile.write("%s \n" % item)


# ******************************************************************************
# Calculating Bosonic holoraumy matrices for given Adinkra
def nicely_print_boson(holo_mats, rmats, pset_arg, adink_def):
	""" holo_mats - List of lists w/ each containg 6 Bosnic Holoraumy matrices
		rmats 	  -	List w/ lists, each containing 4 R matrices
		pset_arg  - String specifying P slices of library, ie P1 or P6
		adink_def - List contains tuples of boolean factors, matrix slice
	"""

	adinkdef_yn	= 0
	if len(adink_def) > 0:
		adinkdef_yn = 1

	holos = [np.asarray(x) for x in holo_mats]
	text_list	= []

	# print("# ********************************")
	text_list.append("#********************************")
	text_list.append("BC4 Coxeter Group Small Library")
	if adinkdef_yn:
		text_list.append(pset_arg + " set: " + adink_def[0][1])
	else:
		text_list.append(pset_arg + " set")
	text_list.append("Calculated Bosonic Holoraumy matrices and R matrices")
	text_list.append("")
	# print("Bosonic Holoraumy matrices for: ", pset_arg)
	# print("")
	lenh, lenr = len(holos), len(rmats)
	print("Length holo_mats: ", lenh)
	print("Length r_matrices: ", lenr)
	print("")

	np.set_printoptions(precision=2, suppress=True, linewidth=100)
	ij_ind	= list(itertools.combinations([0,1,2,3], 2))
	setlen = 0
	# if len(holos) == len(rmats):
	if lenh == lenr:
		setlen = lenh
	else:
		print("LENGTH MISMATCH ERROR")

	for zi in range(0, lenh):
		temph = holos[zi]
		tempr = rmats[zi]
		text_list.append("#********************************")
		text_list.append("Adinkra # " + str(zi))
		if adinkdef_yn:
			print("Def ", adink_def[zi], " Bool Fct", adink_def[zi][0])
			boolstr = (",").join(['(' + str(ix) + ')' for ix in adink_def[zi][0]])
			boolstr = "{" + boolstr + "}"
			print("Bool String:", boolstr)
			text_list.append("Boolean Factor: " + boolstr)
			# text_list.append("Boolean Factor:" + adink_def[zi][0] + " P-set " + adink_def[zi][1])
		text_list.append("Bosonic Holoraumy Matrices")
		# print("Adinkra #", zi)
		vij_strings	= []
		for ijtup in ij_ind:
			ij_temp		= str(ijtup[0] + 1) + str(ijtup[1] + 1)
			ijstr		= "V_{" + ij_temp + "}"
			vij_strings.append(ijstr)
		""" List comprehension for the above 3line forloop, 4line if counting
			vij_strings list declaration	"""
		'''
		vij_strings = [ "V_{" + (str(ijt[0]+1) + str(ijt[1]+1)) + "}" for ijt in ij_ind]
		'''
		# v13strings = " \t" + vij_strings[0] + " \t \t \t" + vij_strings[1] + \
		# " \t \t \t" + vij_strings[2]
		v13strings = "\t" + vij_strings[0] + "\t\t   " + vij_strings[1] + \
		"\t\t   " + vij_strings[2]
		# v46strings = " \t" + vij_strings[3] + " \t \t \t" + vij_strings[4] + \
		# " \t \t \t" + vij_strings[5]
		v46strings = "\t" + vij_strings[3] + "\t\t   " + vij_strings[4] + \
		"\t\t   " + vij_strings[5]

		""" Convoluted way of printing out numpy matrices	"""
		mat13 		= temph[0:3]
		mat13str 	= [np.array_str(y)[1:-1] for y in mat13]
		tm13		= []
		for matstr in mat13str:
			tm13.append([ix.lstrip() for ix in matstr.split('\n')])
		# print(v13strings)
		text_list.append(v13strings)
		for ix in range(0,4):
			pstr = tm13[0][ix] + " \t" + tm13[1][ix] + " \t" + tm13[2][ix]
			text_list.append(pstr)

		mat46		= temph[3:6]
		mat46str	= [np.array_str(y)[1:-1] for y in mat46]
		tm46		= []
		for matstr in mat46str:
			tm46.append([ix.lstrip() for ix in matstr.split('\n')])
		# print(v46strings)
		text_list.append(v46strings)

		for ix in range(0,4):
			pstr = tm46[0][ix] + " \t" + tm46[1][ix] + " \t" + tm46[2][ix]
			# print(pstr)
			text_list.append(pstr)
		text_list.append("")

		"""	Printing R matrices """
		# print("R matrices")
		text_list.append("R matrices")
		for ind, rl in enumerate([ tempr[:2], tempr[2:]]):
			rltostr = [np.array_str(y)[1:-1] for y in rl]
			rtm	= []
			for matstr in rltostr:
				rtm.append([ix.lstrip() for ix in matstr.split('\n')])

			if ind == 0:
				# print("Length: ", len(rtm[0
				len1, len2 	= len(rtm[0][0]), len(rtm[1][0])
				lbl_str	= ""
				if len2 > len1:
					lbl_str	= (len1//3)*" " + "R1" + len2*" " + "R2"
				elif len1 > len2:
					lbl_str	= (len1//2)*" " + "R1" + len2*" " + " R2"
				elif len1 == len2:
					lbl_str = (len1//2)*" " + "R1" + len2*" " + "R2"
				text_list.append(lbl_str)
			elif ind == 1:
				len1, len2 	= len(rtm[0][0]), len(rtm[1][0])
				# print(len1, len2)
				lbl_str		= "  R3" + len2*" " + " R4"
				# print(lbl_str)
				text_list.append(lbl_str)
			for ix in range(0,4):
				pstr = rtm[0][ix] + "\t" + rtm[1][ix]
				# print(pstr)
				text_list.append(pstr)
		text_list.append("")

	pmatsfile = "BC4-CoxeterGroup " + pset_arg + ".txt"
	# pmatsfile = "BC4-CG-BosonicH " + pset_arg + ".txt"
	# pmatsfile = "BC4-CG " + pset_arg + "-BosonicHolos.txt"
	with open(pmatsfile, "w") as wfile:
		for item in text_list:
			wfile.write("%s \n" % item)

# ******************************************************************************
# Calculating Fermionic Holoraumy matrices for given Adinkra
def nicely_print_fermi(fermi_mats, rmats, pset_arg, adink_def):

	""" ***CHECK THIS*** I'm not sure this is even correct to do,
		since every x in fermi_mats is a list of 6 numpy.matrices
	"""
	fermis	= [np.asarray(x) for x in fermi_mats]
	text_list = []

	adinkdef_yn	= 0
	if len(adink_def) > 0:
		adinkdef_yn = 1

	text_list.append("#********************************")
	text_list.append("BC4 Coxeter Group Small Library")
	if adinkdef_yn:
		text_list.append(pset_arg + " set: " + adink_def[0][1])
	else:
		text_list.append(pset_arg + " set")
	text_list.append("Calculated Fermionic Holoraumy matrices and R matrices")
	text_list.append("")
	lenh, lenr = len(fermis), len(rmats)
	print("Length holo_mats: ", lenh)
	print("Length r_matrices: ", lenr)
	print("")

	np.set_printoptions(precision=2, suppress=True, linewidth=100)
	ij_ind	= list(itertools.combinations([0,1,2,3], 2))
	setlen = 0
	if lenh == lenr:
		setlen = lenh
	else:
		print("LENGTH MISMATCH ERROR")

	for zi in range(0, lenh):
		temph = fermis[zi]
		tempr = rmats[zi]
		text_list.append("#********************************")
		text_list.append("Adinkra # " + str(zi))
		if adinkdef_yn:
			print("Def ", adink_def[zi], " Bool Fct", adink_def[zi][0])
			boolstr = (",").join(['(' + str(ix) + ')' for ix in adink_def[zi][0]])
			boolstr = "{" + boolstr + "}"
			print("Bool String:", boolstr)
			text_list.append("Boolean Factor: " + boolstr)
			# text_list.append("Boolean Factor:" + adink_def[zi][0] + " P-set " + adink_def[zi][1])
		text_list.append("Fermionic Holoraumy Matrices")

		vij_strings	= []
		for ijtup in ij_ind:
			ij_temp		= str(ijtup[0] + 1) + str(ijtup[1] + 1)
			ijstr		= "~V_{" + ij_temp + "}"
			vij_strings.append(ijstr)
		'''
		vij_strings = [ "V_{" + (str(ijt[0]+1) + str(ijt[1]+1)) + "}" for ijt in ij_ind]
		'''
		# v13strings = " \t" + vij_strings[0] + " \t \t \t" + vij_strings[1] + \
		# " \t \t \t" + vij_strings[2]
		v13strings = "\t" + vij_strings[0] + "\t\t   " + vij_strings[1] + \
		"\t\t   " + vij_strings[2]
		# v46strings = " \t" + vij_strings[3] + " \t \t \t" + vij_strings[4] + \
		# " \t \t \t" + vij_strings[5]
		v46strings = "\t" + vij_strings[3] + "\t\t   " + vij_strings[4] + \
		"\t\t   " + vij_strings[5]

		""" Obtuse way of nicely printing out n>1 numpy matrices per row	"""
		mat13 		= temph[0:3]
		mat13str 	= [np.array_str(y)[1:-1] for y in mat13]
		tm13		= []
		for matstr in mat13str:
			# onemat = [ix.lstrip() for ix in matstr.split('\n')]
			tm13.append([ix.lstrip() for ix in matstr.split('\n')])
		text_list.append(v13strings)
		# print(v13strings)
		for ix in range(0,4):
			pstr = tm13[0][ix] + " \t" + tm13[1][ix] + " \t" + tm13[2][ix]
			text_list.append(pstr)

		mat46		= temph[3:6]
		mat46str	= [np.array_str(y)[1:-1] for y in mat46]
		tm46		= []
		for matstr in mat46str:
			tm46.append([ix.lstrip() for ix in matstr.split('\n')])
		# print(v46strings)
		text_list.append(v46strings)

		for ix in range(0,4):
			pstr = tm46[0][ix] + " \t" + tm46[1][ix] + " \t" + tm46[2][ix]
			text_list.append(pstr)
		text_list.append("")
		'''*** End upgrade here***'''

		"""	Printing R matrices """
		print("R matrices")
		for ind, rl in enumerate([ tempr[:2], tempr[2:]]):
			rltostr = [np.array_str(y)[1:-1] for y in rl]
			rtm	= []
			for matstr in rltostr:
				rtm.append([ix.lstrip() for ix in matstr.split('\n')])

			if ind == 0:
				# print("Length: ", len(rtm[0
				len1, len2 	= len(rtm[0][0]), len(rtm[1][0])
				lbl_str	= ""
				if len2 > len1:
					lbl_str	= (len1//3)*" " + "R1" + len2*" " + "R2"
				elif len1 > len2:
					lbl_str	= (len1//2)*" " + "R1" + len2*" " + " R2"
				elif len1 == len2:
					lbl_str = (len1//2)*" " + "R1" + len2*" " + "R2"
				text_list.append(lbl_str)
			elif ind == 1:
				len1, len2 	= len(rtm[0][0]), len(rtm[1][0])
				lbl_str		= "  R3" + len2*" " + " R4"
				text_list.append(lbl_str)
			for ix in range(0,4):
				pstr = rtm[0][ix] + "\t" + rtm[1][ix]
				text_list.append(pstr)
		text_list.append("")

	pmatsfile = "BC4-CG-FermionicH " + pset_arg + ".txt"
	with open(pmatsfile, "w") as wfile:
		for item in text_list:
			wfile.write("%s \n" % item)

# ******************************************************************************
# Write results to text file
def write_results(result_file, write_list):
	with open(result_file, "w") as wfile:
		for item in write_list:
			wfile.write("%s \n" % item)

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
