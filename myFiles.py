# -*- coding: utf-8 -*-
#########################################
##### @author: Madeline Galbraith #####
#####   Last Modified: Oct/16/2015   #####
#########################################
## This class gets the files from user ##
####  then returns the main program  ####
#########################################
class myFiles:
	
	def getfilename(self):
		
		#print "What is the name of the *.psf or *.pdb file for your system?"
		#filename = raw_input("> ")
		
		filename = "methyl270.a.0.psf"### use for debugging
		return filename
		
	def getfilelist(self):
		#print "What is the rootname of the file? (i.e. for dynamics1.d.200.dcd give dynamics1)"
		#rootname = raw_input("> ")
		rootname = "methyl270.e."	

		filelist = []
		for x in range(30,51):
			filelist.append(rootname + str(x) + "700.coor.dcd")
			print(x)
			print(filelist)	
		return filelist

	def initializeOutputFiles(self):
		###### put in headers for the files we will print to ########
		list = range(1,31) ## puts the numbers 1 to 30 in the list
		with open('Rgroup_Tetrahedrality.dat','w') as file:
			#file.write('%s \n ' %"timestep \t\t\t\t\t\t\t unit")
			#file.write("\t")
			#file.write("\t")
			file.write(" \t" .join(str(p) for p in list))
			file.write("\n")
		with open('Carbonyl_Tetrahedrality.dat','w') as file:
			#file.write('%s \n ' %"timestep \t\t\t\t\t\t\t unit")
			#file.write("\t")
			#file.write("\t")
			file.write(" \t" .join(str(p) for p in list))
			file.write("\n")
		with open('Amide_Tetrahedrality.dat','w') as file:
			#file.write('%s \n ' %"timestep \t\t\t\t\t\t\t unit")
			#file.write("\t")
			#file.write("\t")
			file.write(" \t" .join(str(p) for p in list))
			file.write("\n")
		with open('RgroupNumberOfNearWaters.dat','w') as file:
			#file.write('%s \n ' %"timestep \t\t\t\t\t\t\t unit")
			#file.write("\t")
			#file.write("\t")
			file.write(" \t" .join(str(p) for p in list))
			file.write("\n")
		with open('CarbonylNumberOfNearWaters.dat','w') as file:
			#file.write('%s \n ' %"timestep \t\t\t\t\t\t\t unit")
			#file.write("\t")
			#file.write("\t")
			file.write(" \t" .join(str(p) for p in list))
			file.write("\n")
		with open('AmideNumberOfNearWaters.dat','w') as file:
			#file.write('%s \n ' %"timestep \t\t\t\t\t\t\t unit")
			#file.write("\t")
			#file.write("\t")
			file.write(" \t" .join(str(p) for p in list))
			file.write("\n")
		with open('WaterCountPerMonomer.dat') as file:
			file.write("," %("timestep"))
			file.write("," .join(str(p) for p in list))
			file.write("\n")
		with open('mainMoleculeVolume.dat','w') as file:
			file.write('%s \n ' %"timestep,volume")
		with open('FirstHydrationShell.dat','w') as file:
			file.write('%s \n' %"timestep,volume,count")
		with open('SecondHydrationShell.dat','w') as file:
			file.write('%s \n' %"timestep,volume,count")
		with open('PMV.dat','w') as file:
			file.write('%s \n' %"timestep,PMV")
		with open("rg.dat","w") as file:
			file.truncate()		
	def initializeIndividualVolumeFiles(self,ts):
		with open("Frame.dat","w") as file:
			file.write("The frame = %d" %ts)
		with open("FirstShellVolumesAtAvg.dat","w") as file:
			file.write("%s \n" %"Index,Volume")		
		with open("SecondShellVolumesAtAvg.dat","w") as file:
			file.write("%s \n" %"Index,Volume")
		with open("BulkVolumesAtAvg.dat","w") as file:
			file.write("%s \n" %"Index,Volume")
		with open("mainMoleculeVolumesAtAvg.dat","w") as file:
			file.write("%s \n" %"Index,Volume")
	def writeRadius(self,num):
		with open("rg.dat","a") as file:
			file.write("%f \n" %num)
