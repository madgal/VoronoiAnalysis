# -*- coding: utf-8 -*-
#########################################
##### @author: Madeline Galbraith #####
#####   Last Modified: Oct/16/2015   #####
#########################################
##### Purpose: 
#####   Analyze molecular dynamics. 
##### 	Calculate the hydration properties of a polymer or 
#####   protein solvated in water.
#####   A psf or pdb file and a trajectory (*.dcd) file are needed. 
#####   Voronoi polyhedra and 
#####   tetrahedrality parameters
#####   are used. 
#############################################
###   Fix on 1/12/2017:
###     added a wrapping  so that tesellations
###     are calculated in a periodic box
############################################
### Addition on 4/6/2017
### determine the number of waters associated with each side chain
############################################               

#### import the following to analyze trajectory files
###from MDAnalysis.analysis import * 
import numpy
from MDAnalysis import *
## need to import so files can be read to and from
import os 

#### importing classes that are needed for analysis ####
import myFiles
import Selections
import Weights
import pyvoroTesselations
###import Waters
#######################################
### initialize Output files #####
myFiles.myFiles().initializeOutputFiles()

######## Main definition for analyzing  ##################
### creating objects for the myFiles class
topology = myFiles.myFiles()
trajectory = myFiles.myFiles()

### getting the filenames
topFilename = topology.getfilename()
trajFilename = trajectory.getfilelist()

### creating the universe from the *psf or *pdb file 
### and a list of trajectory files (*dcd)
universe = Universe(topFilename, trajFilename)

#### get the selections to be used (in this case the polymer)###
selectionOfinterest = Selections.Selections().getMethyl_SubstituentSelection()
atomsOfInterest = universe.selectAtoms(selectionOfinterest)

### get the water selection to be used ####
selectionOfWater = Selections.Selections().getWaterSelection(selectionOfinterest)
moleculesOfWater = universe.selectAtoms(selectionOfWater)

### create a single selection with all the water and the selection of interest ####
waterAndAtomsOfInterest = atomsOfInterest + moleculesOfWater 
allIndices = waterAndAtomsOfInterest.indices

##### get the weights ####
#### weighting based on VanDerWaalsWeighting ####
### these will not change throughout the simulation ####
#weightingValues = Weights.Weights().VanDerWaalsWeighting(waterAndAtomsOfInterest)
weightingValues = Weights.Weights().LJ_sigma_Weighting(waterAndAtomsOfInterest)
### Uncomment to make sure that all atoms have been given a weight
#if len(weightingValues)!= len(allIndices):
#	print "Error in creating weights"
	
#### get the selection to calculate number of waters around each functional group
###Rgroupselection = Selections.Selections().getRGroupSelection(universe)
###carbonylOxygenSelection = Selections.Selections().getCarbonylSelection(universe,selectionOfinterest)
###amideNitrogenSelection = Selections.Selections().getAmideSelection(universe)

### select atoms from the selection above
###Rgroup= universe.selectAtoms(Rgroupselection)
###carbonylOxygen= universe.selectAtoms(carbonylOxygenSelection)
###amideNitrogen= universe.selectAtoms(amideNitrogenSelection)
 
#### create object for the Water class to get the number of waters
###RgroupWaters = Waters.Waters(Rgroup, universe)
###carbonylOxygenWaters = Waters.Waters(carbonylOxygen, universe)
###amideNitrogenWaters = Waters.Waters( amideNitrogen, universe)

#### initialize counts and values to zero for the hydration shell visualization
firstShellCount=[]
pyvoroclass=[]
index=0
rg = myFiles.myFiles()
for ts in universe.trajectory[::50]:
	print(ts)
	### printing out the radius of gyration
	rg.writeRadius(atomsOfInterest.radiusOfGyration())
	### wrap the system 
	waterAndAtomsOfInterest.pack_into_box()
	## get the coordinates for the system
	allCoordinates = waterAndAtomsOfInterest.coordinates()
        minx = min(allCoordinates[:,0])
        miny = min(allCoordinates[:,1])
        minz = min(allCoordinates[:,2])
        maxx = max(allCoordinates[:,0])
        maxy = max(allCoordinates[:,1])
        maxz = max(allCoordinates[:,2])
	### creating box boundaries using the  max and min coordinates in the xyz planes
	### needed for pyvoro tesselations
        box = [[minx,miny,minz],[maxx,maxy,maxz]]

	### put the box into the correct configuration for calculating tetrahedralities
        BoxForWaters=numpy.array([box[0][0]-box[1][0],box[0][1]-box[1][1],box[0][2]-box[1][2]])

	######### calculating the number of waters and tetrahedrality for the functional groups
###	RgroupWaters.getnumWaters(selectionOfWater,4.5)
###	RgroupWaters.printToTable("RgroupNumberOfNearWaters.dat",ts.frame, RgroupWaters.numOFWater)
###	RgroupWaters.tetrahedrality(moleculesOfWater, ts, BoxForWaters)
###	RgroupWaters.printToTable("Rgroup_Tetrahedrality.dat",ts.frame, RgroupWaters.SgParameter)

###	carbonylOxygenWaters.getnumWaters(selectionOfWater,3.5)
###	carbonylOxygenWaters.printToTable("CarbonylNumberOfNearWaters.dat",ts.frame, carbonylOxygenWaters.numOFWater)
###	carbonylOxygenWaters.tetrahedrality(moleculesOfWater, ts, BoxForWaters)
###	carbonylOxygenWaters.printToTable("Carbonyl_Tetrahedrality.dat",ts.frame, carbonylOxygenWaters.SgParameter)

###	amideNitrogenWaters.getnumWaters(selectionOfWater,3.5)
###	amideNitrogenWaters.printToTable("AmideNumberOfNearWaters.dat", ts.frame, amideNitrogenWaters.numOFWater)
###	amideNitrogenWaters.tetrahedrality(moleculesOfWater, ts, BoxForWaters)
###	amideNitrogenWaters.printToTable("Amide_Tetrahedrality.dat", ts.frame, amideNitrogenWaters.SgParameter)

	#### putting the box  dimensions into the correct configuration for calculating tesselations
	## [[xmin,xmax],[ymin,ymax],[zmin,zmax]]
        boxsize = [[box[0][0],box[1][0]],[box[0][1],box[1][1]],[box[0][2],box[1][2]]]
	
	### determine number of atoms in each selection
	Np = len(atomsOfInterest)# number of atoms in polymer
	Nw = len(moleculesOfWater)	## number of waters selected
	Na = len(waterAndAtomsOfInterest) ## Na + Np ( total number of atoms in system)
	
	########### The data must first be passed to the constructor in the class #####
	###### the tesselations must be created next before any other methods can be called 
	###### then the other methods can be called ####
	pyvoroclass.append(pyvoroTesselations.pyvoroTesselations(ts.frame,Np, Nw, Na, boxsize, allCoordinates, weightingValues, allIndices))
	pyvoroclass[index].voronoiTesselations()
	### make sure the tesselations were created before calling the methods on them
	tesWasCreated = pyvoroclass[index].created
	if( tesWasCreated):
		#####   calculate the polymer volume
		polymerVolume= pyvoroclass[index].polymerVolume()
		
		#### find first hydration shell(calculate volume and number of waters)
		firstHydrationShell = pyvoroclass[index].firstHydrationShellVolumeCount()
		### add the number of waters in the first shell to the list
		### this will be used to find the hydration shell most similar to the average for visualization
		firstShellCount.append(firstHydrationShell[1])	
		
		#### find second hydraiton shell(calculate volume and number of waters)
		secondHydrationShell = pyvoroclass[index].secondHydrationShellVolumeCount()
		
		#### calculate the partial Molar Volume
		pmv= pyvoroclass[index].partialMolarVolume()

		#### calculate the number of Waters per index to be used in water per monomer calculations 
		wc= pyvoroclass[index].watersPerMonomerCount()



	#### we reduce the amount of data stored to save memory 
	#### by finding the average of this group of fifty 
	if index ==50:
		### first we find the avg value  ###
                firstShellCount = numpy.array(firstShellCount)
                avgFirstShellCount = numpy.mean(firstShellCount)
                firstShellCount = abs(firstShellCount - avgFirstShellCount)
                minfirst = min(firstShellCount)
                indexFirst = numpy.where( firstShellCount == minfirst)
                indexForAvg = indexFirst[0][len(indexFirst)-1]
                ### then we remove all of the old values ###
                ### and reinitiate it to a list of avg value ###
                firstShellCount = [firstShellCount[indexForAvg]]
                pyvoroclass = [pyvoroclass[indexForAvg]]
                index =0
	###### increment the index for accessing the object's list
	index+=1

#### find the average number of waters in the first shell ###
firstShellCount = numpy.array(firstShellCount)
avgFirstShellCount = numpy.mean(firstShellCount)
firstShellCount = abs(firstShellCount - avgFirstShellCount)
minfirst = min(firstShellCount)
### find the index pertaining to the average number of waters ###
indexFirst = numpy.where( firstShellCount == minfirst)
indexForAvg = indexFirst[0][len(indexFirst)-1]
myFiles.myFiles().initializeIndividualVolumeFiles(pyvoroclass[indexForAvg].ts)
### output the indices for visualization and the volumes for each index
pyvoroclass[indexForAvg].getFirstShellIndividualVolumes()
pyvoroclass[indexForAvg].getSecondShellIndividualVolumes()
pyvoroclass[indexForAvg].getBulkWaterIndividualVolumes()
pyvoroclass[indexForAvg].getPolymerIndividualVolumes()

