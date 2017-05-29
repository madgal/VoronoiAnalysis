# -*- coding: utf-8 -*-
#########################################
##### @author: Madeline Galbraith #####
#####   Last Modified: Oct/16/2015   #####
#########################################
##   This class returns the selection  ##
###  based upon method the user calls ###
###       Each method is for a        ###
###        particular molcule         ###
#########################################
from MDAnalysis import *
class Selections:

	""" 
	def SampleSelection(self):
		### this is an example ###
		return "your selection"
	"""

        def getThisSelection(self):
               ### added by EAP ###
               return "resnum 1:20"
	
	def getPNIPAMSelection(self):
	    #### for PNIPAM ####
		return "bynum 0:575"
		
	def getTBUTYL_SubstituentSelection(self):
		#### for tbutyl substituent ####
		return "bynum 0:665"

	def getPropyl_SubstituentSelection(self):
		#### for tbutyl substituent ####
		return "bynum 0:485"
		
	def getMethyl_SubstituentSelection(self):
		#### for methyl substituent ####
		return "bynum 0:395"
		
	def getWaterSelection(self, SelectedNonWaters):
		####          for water selection          ####
		####    returns all other oxygen atoms     ####
		#### in system not part of other selection ####
		return "name O* and not(" + SelectedNonWaters+")"
		
	def getRGroupSelection(self, universe):
		RSelection=""
		selection = universe.bonds.selectBonds(('N','C1'))
		selection = selection.atomgroup_intersection(universe)
		selection = selection.dump_contents()
		for atom in range(0,len(selection)):
			if atom !=len(selection)-1:
				RSelection += "(bynum  " + str(selection[atom][1]+1)+") or "
			else:
				RSelection += "(bynum  " + str(selection[atom][1]+1)+")"
				
		return RSelection
		
	def getAmideSelection(self, universe):
		return "name N*"

	def getCarbonylSelection(self, universe,selection):
		return "(name O*) and " + selection
