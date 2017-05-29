# -*- coding: utf-8 -*-
#########################################
##### @author: Madeline Galbraith #####
#####   Last Modified: Oct/16/2015   #####
#########################################
#####  This class creates voronoi 
##### tesselations using the 
##### pyvoro package 
#########################################
##### Get the number of waters per monomer
#######################################
import pyvoro

class pyvoroTesselations:

	def __init__(self,ts,Np, Nw, Na, boxsize, allCoordinates, weightingValues,indices):
		#### the constructor### 
		### we set all the values to the global variables ###
		self.ts =ts
		self.Np = Np  ### number of atoms in polymer
		self.Nw = Nw  ### number of waters 
		self.Na = Na  ### number of atoms in system( polymer and waters)
		self.boxsize=boxsize
		self.allCoordinates = allCoordinates
		self.weightingValues = weightingValues
		self.indices = indices
		
		### add the keys that we can use to access various data for the voronoi tesselations
		self.vorKeys = ['volume','faces','adjacency','original','vertices']
		### we care only about keys at indices 0(volume),1(faces),3(original)
		### volume gives volume of Voronoi polyhedra
		### faces gives the voronoi polyhedra faces and
		### used in conjuction with 'adjacent cell' gives the index of the adjacent cell
 		### original gives info about the Voronoi polyhedra cell and its center atom
		
		### Shell[index] = which shell the atom is occupying
		### they can be in the polymer/protein, first shell, second shell, or bulk
		Shells = []
		for x in range(0,Na):
			Shells.append("B") #assuming every atom is bulk water for the moment
	
		self.Shells = Shells

		### create a waternumber vector to keep track of the average number of waters per polymer atom

		waterCount ={}
		for x in range(0,Np):
			waterCount[x] =0

		self.waterCount = waterCount

		### field that states whether the tesselations were created is initialized to false
		self.created = False
		
	def voronoiTesselations(self):
		################
		### create the voronoi tesselation
		################
		try:
			self.vor = pyvoro.compute_voronoi(self.allCoordinates,self.boxsize, 2.0, radii = self.weightingValues, periodic = [True, True, True])
			## if the tesselation is created this becomes true ###
			self.created = True
			### unpack the volumes into a list for ease of further calculations
			self.__getVolumes()
			### determine the shell each atom belongs to 
			self.__determineShells()
				
		except Exception as e: 
			#### if there was an error in input, etc
			## the tesselation is not created
			self.created=False
			self.vor = 0
		
	def __getVolumes(self): 
		#### determine the volumes of all atoms
		self.Volume = []
		for atom in self.vor: 
			#add the volume of this atom to the list
			self.Volume.append(atom[self.vorKeys[0]])

	def __determineShells(self):
	
		#### Note that if i<=P than the atom is part
		#### of the selection we are interested in 
		for atom in range(0,self.Np):
			self.Shells[atom] = "P"					
					
		##### get waters that are in #####
		##### the first shell #####
		for atom in range(0,self.Na):
			if( atom>= self.Np):
				inFirstShell =False
				for x in self.vor[atom][self.vorKeys[1]]:
					### if the adjacent cell is part of the polymer
					## this atom is in the first shell
					if x['adjacent_cell'] < self.Np:
						inFirstShell =True
				if inFirstShell:
					self.Shells[atom] = "S1"
		##### get waters that are in #####
		##### the second shell #####				
		for atom in range(0, self.Na):
			### if the atom has not been assigned to a shell yet, determine it's shell
			if atom>=self.Np and self.Shells[atom]!="S1":
				inSecondShell = False
				for x in self.vor[atom][self.vorKeys[1]]:
					### if the acjacent cell is in the first shell
					## this atom is in the second shell
					if self.Shells[x['adjacent_cell']] == "S1":
						inSecondShell = True
						
				if inSecondShell:
					self.Shells[atom] = "S2"
		#### if the atom's shell has not been modified than it is in bulk
		
	def polymerVolume(self):
		### calculate the polymer's volume 
		polymerVolume=0
		for i in range( 0, len(self.Shells)):
			if(self.Shells[i] =='P'):
				polymerVolume += self.Volume[i]
				
		self.polymerVolume = polymerVolume
		self.__printToFile("mainMoleculeVolume.dat",polymerVolume)
		return polymerVolume
	
	def firstHydrationShellVolumeCount(self):
		### calculate the first shell's volume and number of waters in first shell
		firstShellVolume=0
		firstShellCount=0
		for i in range( 0, len(self.Shells)):
			if(self.Shells[i] == 'S1'):
				firstShellVolume += self.Volume[i]
				firstShellCount +=1
		
		self.firstShellVolume = firstShellVolume
		self.firstShellCount = firstShellCount
		self.__printToFile("FirstHydrationShell.dat", [firstShellVolume, firstShellCount])
		return [firstShellVolume, firstShellCount]
	
	def secondHydrationShellVolumeCount(self):
		### calculate the second shell's volume and number of waters in second shell
		secondShellVolume = 0
		secondShellCount =0
		for i in range(0, len(self.Shells)):
			if(self.Shells[i] == 'S2'):
				secondShellVolume += self.Volume[i]
				secondShellCount +=1
				
		self.secondShellVolume = secondShellVolume
		self.secondShellCount = secondShellCount
		self.__printToFile("SecondHydrationShell.dat",[secondShellVolume, secondShellCount])
		return [secondShellVolume, secondShellCount]

	def __dist(self,r1,r2):
		x= (float(r1[0])-float(r2[0]))**2
		y= (float(r1[1])-float(r2[1]))**2
		z= (float(r1[2])-float(r2[2]))**2
		return (x + y + z)
	
	def watersPerMonomerCount(self):
		### determine the closest polymer atom to each water molecule within the first shell
		waterCount={}
		for i in range(0,self.Np):
			waterCount[i]=0
		
		nearNeigh = {}
		for i in range(0, len(self.Shells)):
			if(self.Shells[i] == 'S1'):
				junklist = []
				for nn in self.vor[i][self.vorKeys[1]]:
					nn_index = nn['adjacent_cell']
					if nn_index not in junklist:
						if self.Shells[nn_index]=='P': 
							junklist.append(nn['adjacent_cell'])
				nearNeigh[i]=junklist
		for i in nearNeigh:
			original = self.vor[i][self.vorKeys[3]]
			min=100
			min_index = -1
			for j in nearNeigh[i]:
				dist = self.__dist(self.vor[j][self.vorKeys[3]],original)
				if dist < min:
					min=dist
					min_index = j	
			nearNeigh[i] = min_index
				
		for i in nearNeigh:
			if nearNeigh[i]!= -1:
				waterCount[nearNeigh[i]] = waterCount[nearNeigh[i]] +1
			

		self.waterCount = waterCount
		self.__printToFileType2("WaterCountPerMonomer.dat",waterCount)
		return [waterCount]
	
	def partialMolarVolume(self):
		### calculate the average volume of a bulk water molecule ###
		bulkVolume =0
		bulkCount =0
		for i in range(0, len(self.Shells)):
			if( self.Shells[i] == 'B'):
				bulkVolume += self.Volume[i]
				bulkCount +=1
		avgbulkVolume = bulkVolume/(bulkCount*1.0) #### get the average volume of a single water in bulk
		
		## use the avgBulkVolume and multiply by the number of waters in first shell 
		## this gives the bulk volume for a number of waters (N)
		NbulkVolume = avgbulkVolume * self.firstShellCount  
		
		### calculat the partial molar volume using equation from Voloshin paper
		### this will be cited
		PMV =  (self.polymerVolume + self.firstShellVolume - NbulkVolume )
		self.PMV = PMV
		self.__printToFile("PMV.dat",PMV)
		return PMV
	
	def __printToFile(self,filename, data):	
		if type(data) is float:
			with open(filename, 'a') as file:
				file.write("%d," %self.ts)
				file.write("%f \n" %data)
		elif type(data) is list:
			with open(filename, 'a') as file:
				file.write("%d," %self.ts)
				file.write("%f," %data[0])
				file.write("%f \n" %data[1])
		#else:
		#	print "There is an error in your data format"
	
	def __printToFileType2(self,filename, data):	
			with open(filename, 'a') as file:
				file.write("%d," %self.ts)
				for el in data:
					file.write("%f," %data[el])
				file.write("\n")
		
	def getFirstShellIndividualVolumes(self):
		histfile = open("FirstShellVolumesAtAvg.dat","a")
		vmdFile = open("FirstShellIndices_VMD.dat","w")
                for x in range(self.Np,len(self.Shells)): # only going through the waters
	                if self.Shells[x]=='S1':
         	                histfile.write("%d," %self.indices[x])
                                histfile.write("%f \n" % self.Volume[x])
				vmdFile.write("%d " %self.indices[x])
		vmdFile.close()
		histfile.close()

        def getSecondShellIndividualVolumes(self):
                histfile = open("SecondShellVolumesAtAvg.dat","a")
                vmdFile = open("SecondShellIndices_VMD.dat","w")
                for x in range(self.Np,len(self.Shells)): # only going through the waters
                        if self.Shells[x]=='S2':
                                histfile.write("%d," %self.indices[x])
                                histfile.write("%f \n" % self.Volume[x])
                                vmdFile.write("%d " %self.indices[x])
                vmdFile.close()
                histfile.close()

	def getBulkWaterIndividualVolumes(self):
		with open("BulkVolumesAtAvg.dat","a") as file:
			for x in range(self.Np,len(self.Shells)): # only going through the waters
				if self.Shells[x]=='B':
					file.write("%d," %self.indices[x])
					file.write("%f \n" % self.Volume[x])
				
	def getPolymerIndividualVolumes(self):
		with open("mainMoleculeVolumesAtAvg.dat","a") as file:
			for x in range(0,self.Np):
				file.write("%d," %self.indices[x])
				file.write("%f \n" %self.Volume[x])
