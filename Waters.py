# -*- coding: utf-8 -*-
#########################################
##### @author: Madeline Galbraith #####
#####   Last Modified: Oct/16/2015   #####
################################r#########
from MDAnalysis import *
from MDAnalysis.analysis import *
import MDAnalysis.KDTree.NeighborSearch as NS
import numpy
import math
class Waters(object):
	def __init__(self, group, universe):
		self.group = group
		self.universe = universe
	def getnumWaters(self,waterSelection,distSel):
	
		numOFWater = []
		for index in self.group.indices():
			proximalWaterSelection = waterSelection + " and (around "+ str(distSel)+ " bynum " +str(index +1)+")"
			numOFWater.append(len(self.universe.selectAtoms(proximalWaterSelection)))
		self.numOFWater = numOFWater
	
	def printToTable(self,filename,timestep, ListOfInterest):
		with open(filename, 'a') as file:
			#file.write("%d \t \t " %timestep)
			for element in ListOfInterest:
				file.write("%f  " %element)
			file.write("\n")	

			
	def tetrahedrality(self, Waters, ts, box):
		# box is created using the unmodified box dimensions	
		neighborWater = NS.AtomNeighborSearch(Waters)# where Waters is the water atom selectoin
		SgParameter =[]
		for index in self.group.indices():
			radius =3
			selection = self.universe.selectAtoms("bynum " + str(index +1))
			selectionNeighbors = neighborWater.search_list(selection, radius)
			while len(selectionNeighbors) <4:
				radius +=1 
				selectionNeighbors = neighborWater.search_list(selection, radius)
			distanceToSelection = distances.distance_array(selectionNeighbors.coordinates(ts), selection.coordinates(ts), box)
			AllSelectionindices = selectionNeighbors.indices()
			
			distance = []
			indices = []
			for x in range(0,  len(distanceToSelection)):
				distance.append(min(distanceToSelection[x]))
				indices.append(AllSelectionindices[x])
			fourClosest = self.__FindClosest(distance, indices)
			distance = distance[0:4]
			indices = indices[0:4]
			
			central = index 
			identity =0
			cosineSum =0
			for j in range(1,4):
				a = indices[identity]
				ident =identity+1
				for k in range(j+1, 5):
					b = indices[ident]
					sel = "(bynum " + str(a+1) + " ) or (bynum " + str(central +1) +") or (bynum " + str(b +1)+")"
				
					angle = self.universe.selectAtoms(sel)
					angle = angle.angle()
					
					cosine = math.cos(angle)
					cosine = (cosine+ (1.0/3.0) )**2
					
					cosineSum+=cosine
				
					ident+=1
				identity+=1
			SgParameter.append(1.0 - ((3.0/32.0)*cosineSum))
		self.SgParameter = SgParameter
					
	def __FindClosest(self, distance, indices):
		if len(distance)== 4:
			return [distance, indices]
			
		size = len(distance)
		
		for j in range( 1, size):
			temp = distance[j]
			i=j
			while i>0 and distance[i-1] > temp:
				distance[i] = distance[i-1]
				i-=1
			distance[i] = temp
		return [ [distance[0],distance[1],distance[2],distance[3]],[indices[0],indices[1], indices[2],indices[3]]]
			
