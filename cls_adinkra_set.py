# ******************************************************************************
# Name:    Class Adinkra Set
# Author:  Vadim Korotkikh
# Email:   v
# Date:    September 2017
# Version: 1.0A
#
# Description: Intended to be a class for storing sets containing Adinkras
#
# ******************************************************************************
# Library Imports

import numpy as np



# ******************************************************************************
def main():
	# nothing here for now
	pass


class AdinkraSet():

	adinkra_count = 0

	@classmethod
	def return_info(cls):
		return cls.adinkra_count

	def __init__(self, args, *args):

		self.arg = args
