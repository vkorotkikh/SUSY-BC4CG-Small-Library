# ******************************************************************************
# Name:    Calculating BC4 CG Small Library
# Author:  Vadim Korotkikh
# Email:   va.korotki@gmail.com
# Date:    November 2017
# Version: 2.0
#
# Description: Code for verifying the BC4 space coefficient library
#
# ******************************************************************************
# Library Imports

import re, sys, math, time
import itertools, logging
import numpy as np

#>******************************************************************************
# Function Imports
import vij_holoraumy_calc

p_switch = 0
py_ver = ''
if sys.version_info >= (3, 5):
	py_ver = '3.5'
elif sys.version_info >= (2, 6) and sys.version_info < (2, 7):
	py_ver = '2.6'
	print("Python version :", sys.version_info[0:3])
elif sys.version_info >= (2, 7):
	py_ver = '2.7'
	print("Python version :", sys.version_info[0:3])
else:
	raise Exception("Minimum Python 2.6 or Python 3.5 are required for this code.")
logging.basicConfig(level=logging.DEBUG)
bc4logger = logging.getLogger(__file__)

#>******************************************************************************
def main(pset_str):

	print("# ***********************************************************************")
	print("# Name:    Calculating BC4 CG Small Library")
	print("# Author:  Vadim Korotkikh	")
	print("# Email:   va.korotki@gmail.com")
	print("# Date:    June 2017	")
	print("# Version: N/A")
	print("#	")
	print("# ***********************************************************************")
	print("		")
	user_options()
	# bc4cg_holoraumy_calc(pset_str)

#>******************************************************************************
def bc4_validation_organizer(pset_arg, *args):
	''' Based on pset_arg and *args figures out which function to pass
		arguments too.
		PALL, ALL, P1 - P6,  'coef' ->  only fermi
		PALL, P1 - P6, 'mats' - fermi / boson
	'''
	args_tuple	= args
	argslow	= [x.lower() for x in args_tuple]
	# print(argslow, type(args))
	bc4logger.info("Input args: %s" % pset_arg)
	if all(isinstance(argx, str) for argx in args_tuple):
		# print("All *args are strings")
		if len(args_tuple) == 2 and args_tuple[0] == 'coef':
			''' redirect to calculating the coefficients '''
			pass
		elif len(args_tuple) == 2 and 'coef' in args_tuple:
			if 'boson' in args_tuple:
				print("ERROR")
			elif 'fermi' in args_tuple and not 'boson' in args_tuple:
				# pass args to whatever function will do this
				pass

		elif len(args_tuple) == 2 and 'mats' in args_tuple:
			if 'boson' in args_tuple or args_tuple[1] == 'boson':
				bosi = args_tuple.index('boson')
				mati = args_tuple.index('mats')
				# Send input to bc4cg_holo via matching the indices to content

				bc4cg_holoraumy_mats(pset_arg, args_tuple[bosi], args_tuple[mati])
			elif 'fermi' in args_tuple or args_tuple[1] == 'fermi':
				feri = args_tuple.index('fermi')
				mati = args_tuple.index('mats')
				bc4cg_holoraumy_mats(pset_arg, args_tuple[feri], args_tuple[mati])
			else:
				bc4logger.debug("ERROR %s %s" % (pset_arg, args_tuple[1]))
				print("this shouldn't have happened")

		elif len(args_tuple) == 2 and 'Vmats' in args_tuple:
			if 'boson' in args_tuple or args_tuple[1] == 'boson':
				bosi = args_tuple.index('boson')
				mati = args_tuple.index('mats')
				bc4cg_holoraumy_mats(pset_arg, args_tuple[bosi], args_tuple[mati])
			elif 'fermi' in args_tuple or args_tuple[1] == 'fermi':
				ferind = args_tuple.index('fermi')
				matind = args_tuple.index('mats')
				# bc4cg_holoraumy_mats(pset_arg, matind, bind)
			else:
				bc4logger.debug("ERROR %s %s" % (pset_arg, args_tuple[1]))
				print("this shouldn't have happened")
		else:
			# What goes here? Technically this shuldn't happen. Unless I change
			# the code later
			pass


#>******************************************************************************
def bc4cg_holoraumy_mats(pset_arg, *args):
	'''
	Process organizer
	ITT, pset is tracked via P1, P2, P3....so pset_arg must be a Pn or
	it must be made one. Otherwise function will not executed
	'''
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
				pass
				# print("bc4cg_holoraumy_mats-ERROR")  # configure this better later
	""" Generate individual Library slices	"""
	for ps in psets_list:
		# temp_tetrad	= tetrad_setgen(ps)
		temp_tetrad = tetrad_setgen_detailed(ps)
		if ps not in psets_dict:
			psets_dict['%s' % ps] = temp_tetrad
		else:
			print("Unknown P set -  ERROR")

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
			# print("		")
			vij_holoraumy_calc.calc_holoraumy_mats(plist, pset, holotype)
			print("Holoraumy calc. for:", pset,"finished")
			print("")

#>******************************************************************************
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
	return pset_boold

#>******************************************************************************
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

#>******************************************************************************
def bc4cg_libsets(p_index):
	""" Defining the P set slices of the original BC4 CG Adinkra library	"""

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
def binaries(bin_code):
	""" Defining the booleans as n=4 vector of 1s, corresponding to boolean
	value. 0 - all 1, no flipping done. 14 - [ 1, -1, -1, -1 ]	"""
	bin_code = int(bin_code)
	binaries_lt	= [(0, [1,1,1,1]), (2, [1,-1,1,1]), (4, [1,1,-1,1]),
					(6, [1,-1,-1,1]), (8, [1,1,1,-1]), (10, [1,-1,1,-1]),
					(12, [1,1,-1,-1]), (14, [1,-1,-1,-1])]
	for btuple in binaries_lt:
		if bin_code == btuple[0]:
			return btuple[1]

##************************************
def flips_org_lib(flip_set):
	""" Hardcoded definitions of the boolean 'factors' used for each P
		slice of the original BC4 Coxeter Group Small Library	"""
	p_slice 	= {}

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

#>******************************************************************************
def lmat_flipping(vbasis, binaries_list):
	""" Use the binary representation values to perform flips on L matrices in
	each Adinkra	"""

	lmat_list		= []

	for xbin in binaries_list:
		print("Flip:", xbin)
		binmats = [binaries(b) for b in xbin]
		temp	= [np.dot(binmats[i], vbasis[i]) for i in range(0, len(binmats))]
		lmat_list.append(temp)
	return lmat_list

#>******************************************************************************
def verify_input(userstr):
	'''
	Check whether string is PALL or one of P1 - P6 or 7/exit option
	'''
	loc_str = userstr.strip().lower()
	if loc_str.isdigit():
		if int(loc_str) in list(range(1,7)):
			return "P" + loc_str
		elif int(loc_str) == 7:
			return loc_str
		else:
			return 0
	elif loc_str[0].isdigit():
		if int(loc_str[0]) in list(range(1,7)):
			return "P" + loc_str[0]
		elif int(loc_str[0]) == 7:
			return loc_str
		else:
			return 0

	if (loc_str).startswith('p'):
		recomp = re.compile('[pP1-6]', re.IGNORECASE)
		relist = recomp.findall(loc_str)
		if recomp.match(relist[0]) is not None:
			if len(relist) > 1:
				if relist[1].isdigit():
					if int(relist[1]) in list(range(1,7)):
						print('Sucesfull BC4 CG Library matching input')
						return 'P' + relist[1]
		else:
			sys.exit('STRING ISSUE')
	elif not (loc_str).startswith('p'):
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
				sys.exit("NONE ISSUE")


#>******************************************************************************
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
	print(" < 7/Back >  -  Go back")
	userinput = pyver_uinput(": ")
	checkstr  = verify_input(userinput)
	if not checkstr:
		print("INVALID INPUT")
		print("TRY AGAIN")
		pset_options_std()
		pass
	elif checkstr.strip() == '7' or checkstr.lower() == 'back':
		user_options()
	elif checkstr:
		return checkstr
	else:
		pass


# **************************************************************************
def pyver_uinput(instr):
	if int(py_ver[0]) >= 3:
		return `input`(instr)
	elif int(py_ver[0]) < 3:
		return raw_input(instr)

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
	bc4logger.info("Executing user options")
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
		# print(" < 6 >  -  Set output file string")(
		return pyver_uinput(": ")

	# Set loopcount = 0 of no arg is supplied for first time
	def option_one(loopcount=0):
		print("")
		print(" < 1 >  -  Calculate all P-sets")
		print(" < 2 >  -  Calculate select P-set from the Small Library")
		print(" < 3 >  -  Back to main menu")
		ninput = pyver_uinput(": ")
		if ninput.strip() == '1':
			# bc4_validation_organizer('PALL', 'Vmats', 'fermi')
			bc4logger.info(" Calculating all P-set %s" % ('Bosonic Holo. matrices'))
			bc4_validation_organizer('PALL', 'mats', 'boson')
		elif ninput.strip() == '2':
			usr_pset = pset_options_std()
			bc4logger.info(" Calculating %s %s" % (usr_pset, 'Bosonic Holo. matrices'))
			bc4_validation_organizer(usr_pset, 'mats', 'boson')
		elif ninput.strip() == '3':
			option_activator('core')
		else:
			loopcount += 1
			print("Unrecognized option")
			bc4logger.debug("Unrecognized input %s" % ninput)
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
		ninput = pyver_uinput(": ")
		if ninput.strip() == '1':
			# bc4_validation_organizer('PALL', 'Vmats', 'fermi')
			bc4logger.info("Calculating all P-set %s" % ('Fermionic Holo. matrices'))
			bc4_validation_organizer('PALL', 'mats', 'fermi')
		elif ninput.strip() == '2':
			usr_pset = pset_options_std()
			bc4logger.info("Calculating %s %s" % (usr_pset, 'Fermionic Holo. matrices'))
			exit()
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
		print(" -- Work in Progress -- ")
		print(" < 1 >  -  Display entire BC4 CG Small Library")
		print(" < 2 >  -  Display select P-set from the Small Library")
		print(" < 3 >  -  Back to main menu")
		ninput = pyver_uinput(": ")
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
		print(" -- Work in Progress -- ")
		print(" < 1 >  -  Verify entire BC4 CG Library")
		print(" < 2 >  -  Verify select P-set from the Small Library")
		print(" < 3 >  -  Back to main menu")
		ninput = pyver_uinput(": ")
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

	def option_five():
		print("")
		print("Quiting script. Are you sure (yes/no)?")
		ninput = pyver_uinput(": ")
		if ninput.lower() == 'yes':
			sys.exit("EXITING BC4 CG Library Utility")
		elif ninput.lower() == 'no':
			print("Going back to main menu")
			option_activator('core')
		else:
			pass
	# Preset loopcount = 0 if no arg supplied for first time.
	# def option_four(loopcount=0):
	# 	print("")
	# 	print(" < 1 >  -  Calculate all P-sets")
	# 	print(" < 2 >  -  Calculate select P-set from the Small Library")
	# 	print(" < 3 >  -  Back to main menu")
	# 	opt_str = input(": ")

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

	# def option_six():
	# 	print("NOT ACTIVATED (code not finished)")
	# 	print("Going back to main menu")
	# 	option_activator('core')

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
			bc4logger.debug("Unknown input  %s", input_str)
			print("UNRECOGNIZED SELECTION: ", input_str)
			print("Try again...")
			if mcounter <= 6:
				mcounter+=1
				option_activator(core_options(), mcounter)
			elif mcounter == 7:
				print("Too many attempts. Try again later.")
				sys.exit("EXITING BC4 CG Library Utility")

	uinput = core_options()
	option_activator(uinput)


# **************************************************************************
# Execute main()
if __name__ == "__main__":
	start_time = time.time()
	try:
		main(sys.argv[1])
	except IndexError:
		user_options()
	print("-- Execution time --")
	print("---- %s seconds ----" % (time.time() - start_time))
