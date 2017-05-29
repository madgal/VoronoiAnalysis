# -*- coding: utf-8 -*-
##################################################################
##### @author: Madeline Galbraith ################################
#####   Last Modified: Oct/16/2015   ##############################
###################################################################
###################################################################
###################################################################
#####   Updated 1/12/2017 :  fixed how the weights were calculated
###################################################################
###################################################################
from MDAnalysis.analysis import *
class Weights:

	def LJ_sigma_Weighting(self, allSelections):
		### parameters take are Rmin/2 from the OPLS-AA forcefield
		par = {"C135":1.9643086,"C136":1.9643086,"C137":1.9643086,"C145":1.9923701,"C149":1.9643086,"C157":1.9643086,"C158":1.9643086,"C166":1.9923701,"C206":1.9643086,"C209":1.9643086,"C210":1.9643086,"C214":1.9643086,"C223":1.9643086,"C224":1.9643086,"C235":2.1046163,"C242":1.9643086,"C245":1.9643086,"C246":1.9643086,"C267":2.1046163,"C271":2.1046163,"C274":1.9643086,"C283":1.9643086,"C284":1.9643086,"C285":1.9643086,"C292":1.9643086,"C293":1.9643086,"C295":1.9643086,"C296":1.9643086,"C302":1.2627698,"C307":1.9643086,"C308":1.9643086,"C500":1.9923701,"C501":1.9923701,"C502":1.9923701,"C505":1.9643086,"C506":1.9923701,"C507":1.9923701,"C508":1.9923701,"C509":1.9923701,"C510":1.9923701,"C514":1.9923701,"H140":1.4030776,"H146":1.3581791,"H155":0.0000000,"H168":0.0000000,"H204":0.0000000,"H240":0.0000000,"H241":0.0000000,"H270":0.0000000,"H290":0.0000000,"H301":0.0000000,"H304":0.0000000,"H310":0.0000000,"H504":0.0000000,"H513":0.0000000,"N237":1.8240008,"N238":1.8240008,"N239":1.8240008,"N287":1.8240008,"N300":1.8240008,"N303":1.8240008,"N309":1.8240008,"N503":1.8240008,"N511":1.8240008,"N512":1.8240008,"O154":1.7510408,"O167":1.7229792,"O236":1.6612438,"O268":1.6836931,"O269":1.6612438,"O272":1.6612438,"S200":2.0204317,"S202":2.0204317,"S203":2.0204317,"OT":1.768200,"HT":0.224500,"SOD":2.9969737,"CLA":2.2561487}

		list=[]
		for atom in allSelections:
			### add the sigma parameter to the list as its weight
			if atom.name in par:
				list.append(par[atom.name*(2**(5/6))])
		
			     ## should never get here
		 	     ## append a weight of 1 to maintain the correct weight positions
			else:
				list.append(1.00)
			
		return list

	################
	################
	################
	################
	################
	def VanDerWaalsWeighting(self, allSelections):
	
		list=[]
		for atom in allSelections:
		
			## if the name of the atom starts with O
			## assume it is an oxygen
			## add a relative weight of its radii in angstroms
		
			if atom.name[0]=="O":
				list.append(1.52)#angstroms
			elif atom.name[0]=="C":
				list.append(1.70)#angstroms
			elif atom.name[0]=="N":
				list.append(1.55)#angstroms
			elif atom.name[0]=="H":
				list.append(1.20)#angstroms
			elif atom.name[0]=="K":
				list.append(2.75)#angstroms
			elif atom.name[0]=="S":
				list.append(2.27)#angstroms

			     ## should never get here
		 	     ## append a weight of 1 to maintain the correct weight positions
			else:
				list.append(1.00)
			
		return list
