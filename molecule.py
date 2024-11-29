#!/usr/bin/env python
#
# 
# TODO: transform all numerical lists into arrays
# TODO: separate general stuff to a module
#
# Sat Dec  3 17:41:01 EST 2005 - transforming into general 
# 
# Tue Aug  9 18:46:51 EDT 2005 - Molecule class introduced 
#                              - Angular properties, etc. are calulated inside the Molecule class
#                                care should be taken to move the mol to the center of mass orogin befor calulating angular stuff   

import sys
import string
import re
from numpy import * 
import copy  
#from Matrix import * 
from numpy.linalg import *

#from LinearAlgebra import *

#======
#
tokcal = 627.5096 
toang = 0.5291772489999999
       
#
range3 = range(3)
#
# magic coeff - to kcal/mol
C1 = 4.184E-4
#
#======
# various math stuff
# 
def kronecker(a,b):
    if a == b:
       return 1
    else:
       return 0

def cross(v, w):
    "Calcualates cross product "
    x = v[1]*w[2] - v[2]*w[1]
    y = v[2]*w[0] - v[0]*w[2]
    z = v[0]*w[1] - v[1]*w[0]

    return [x, y, z]

#========
# this stuff is stolen from Venus
masses = [1.008,4.003,6.941,9.012,10.810,12.011,
14.007,15.999,18.998,20.170,22.990,24.305,
26.980,28.090,30.974,32.060,35.453,39.948,
39.102,39.962,44.960,47.900,50.941,51.996,
54.938,55.847,59.933,58.710,63.546,65.380,
69.720,72.590,74.922,78.960,79.904,83.890,
85.467,87.620,88.906,91.220,92.906,95.940,
98.906,101.070,102.906,106.040,107.868,112.410,
114.820,118.690,121.750,127.600,129.900,126.900,
131.300,132.905,137.330,138.905,140.120,140.908,
144.240,145.000,150.400,151.960,157.250,158.925,
162.500,164.930,167.260,168.934,173.040,174.970,
178.490,180.947,183.850,186.200,190.200,192.220,
195.090,196.966,200.590,204.370,207.200,208.981,
209.000,210.000,222.000,223.000,226.025,227.000,
232.038,231.036,238.029,237.048,244.000,243.000,
247.000,247.000,251.000,254.000,257.000,258.000,
259.000,260.000,78.9183,80.9163 ]

symbolAndMass = ['H ','He','Li','Be','B ','C ','N ','O ','F ',
'Ne','Na','Mg','Al','Si','P ','S ','Cl','Al',
'K ','Ca','Sc','Ti','V ','Cr','Mn','Fe','Co',
'Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr',
'Rb','Sr','Y ','Zr','Nb','Mo','Tc','Ru','Rh',
'Pd','Ag','Cd','In','Sn','Sb','Te','Te','I ',
'Xe','Cs','Ba','La','Ce','Pr','Nd','Pm','Sm',
'Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu',
'Hf','Ta','W ','Re','Os','Ir','Pt','Au','Hg',
'Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac',
'Th','Pa','U ','Np','Pu','Am','Cm','Bk','Cf',
'Es','Fm','Md','No','Lr','Br','Br' ]

sym2mass = {}
for i in range(len(symbolAndMass)):
   sym2mass[symbolAndMass[i]] = masses[i]
   sym2mass[symbolAndMass[i].lower()] = masses[i]

#========
# this - from PyQaunte (?)

symbol = [
'X ','H ','He',
'Li','Be','B ','C ','N ','O ','F ','Ne',
'Na','Mg','Al','Si','P ','S ','Cl','Ar',
'K ', 'Ca', 'Sc', 'Ti', 'V ', 'Cr', 'Mn', 'Fe',
'Co', 'Ni', 'Cu', 'Zn',
'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr',
'Rb', 'Sr', 'Y ', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru',
'Rh', 'Pd', 'Ag', 'Cd',
'In', 'Sn', 'Sb', 'Te', 'I ', 'Xe']


sym2no = {}
for i in range(len(symbol)):
   sym2no[symbol[i]] = i
   sym2no[symbol[i].lower()] = i

#------------------------------
# Atom

class Atom:

      def __init__(self):
          "Just init"
          self.n=0
          self.Zz=0
	  self.Symb=""
	  self.mass=0.0
	  # x,y,z coords
          self.q = []
	  # gradients 
          self.dedq = []
	  # velocites 
          self.dqdt = []
	  for i in range(3):
	      self.q.append(0.0)
	      self.dedq.append(0.0)
	      self.dqdt.append(0.0)

      def Rsq(self):
          "Square of the dist to the origin"
	  #return self.x**2+self.y**2+self.z**2
	  return self.q[0]**2+self.q[1]**2+self.q[2]**2

      def coord(self,i):
          if i > 2:
	     print "error in Atom.coord()"
	     sys.exit()
	  #foo = (self.x,self.y,self.z)
	  return self.q[i]

      def InitGeom(self,str):
          "Read the geom"
          L = string.split(str)
          self.n=int(L[0])
	  self.Symb=L[1]
          self.Zz=float(L[2])
          self.q[0] = float(L[3])*toang
          self.q[1] = float(L[4])*toang
          self.q[2] = float(L[5])*toang

      def setGradStr(self,str):
          "Read the grad; depricated; remove it later"
          L = string.split(str)
#          self.n=int(L[0])
          self.dedq[0] = float(L[1])*tokcal/toang
          self.dedq[1] = float(L[2])*tokcal/toang
          self.dedq[2] = float(L[3])*tokcal/toang

      def setGrad(self,gx,gy,gz,convToKcal=False,convToAng=False):
          """Set the grads; the units are not followed yet; 
	  However, setting convToAng and convToKcal to True will do the nesessary transformations """
	  c_ = 1.0
	  if convToKcal: c_ = c_ * tokcal
	  if convToAng: c_ = c_ / toang
          self.dedq[0] = gx * c_
          self.dedq[1] = gy * c_
          self.dedq[2] = gz * c_
      #
      # Verlet-related stuff
      # 
      def setHVerlet(self,dt_):
          """Calculates and remembers the H = dt**2*C1/2m coefficient to be used in Verlet
	     Yes, this is a mess and it should not be here; I need to get rid of `Atom` alltogether. 
	     May be?"""
	  self.HVerlet = dt_ * dt_ * C1 /2 /self.mass   
	  return self.HVerlet

      def VerletEOM(self,dt_):
          """ Does vel Verlet EOM  """
	  for i in range3:
	      self.q[i] += ( self.dqdt[i] * dt_ - self.HVerlet * self.dedq[i] )
	      self.dqdt[i] -= self.HVerlet * self.dedq[i] / dt_ 


      def VerletVelocitiesUpdate(self,dt_,quench_=1.0):
          """ ajusts the velocities half-way  """ 
	  for i in range3:
	      self.dqdt[i] -= self.HVerlet * self.dedq[i] / dt_ * quench_
	      

      def getEnergyKin(self):
          """ returns kinetic energy in kcal/mol  """ 
	  T_ = 0.0
	  for i in range3:
	      T_ += self.mass * self.dqdt[i]**2 / C1 / 2 
          return T_

      def InitExplicit(self,x,y,z,vx=0.0,vy=0.0,vz=0.0,no=0,Symb="",Zz=0,mass=0.0):
          "Init explicitly, give only symbol or charge, mass will be found, velocities optinal"
          self.q[0] = x
          self.q[1] = y
          self.q[2] = z
	  # velocites   
          self.dqdt[0] = vx
          self.dqdt[1] = vy
          self.dqdt[2] = vz
	  # 
          self.n = no
	  #
          self.Zz = Zz
	  self.mass = mass

	  if len(Symb) < 2:
	     Symb = "%s " % Symb 
	  self.Symb = Symb
	  
          # whatever is given, find out the rest
	  if (self.Zz == 0) and (self.Symb != " "):
	      self.Zz = sym2no[Symb]
	  if (self.Zz != 0) and (self.Symb == " "):
	      self.Symb = symbol[Zz]
	  if (self.mass == 0.0) and (self.Symb != " "):
              self.mass = sym2mass[self.Symb]

      def Move(self,dq,ddqdt=[]):
          "Shifts coords according to dq[3] array, if the velocities are given, modify them too"
	  for i in range3:
	      self.q[i] += dq[i]
	  #
	  if (ddqdt!=[]):
	      for i in range3:
	          self.dqdt[i] += ddqdt[i]
          #self.q[1] += dy
          #self.q[2] += dz
	      

      def __str__(self):
           return '%s\t%4i\t%7.3f\t%12.7f\t%12.7f\t%12.7f\t%12.7f\t%12.7f\t%12.7f' % (self.Symb, self.Zz, self.mass, self.q[0], self.q[1], self.q[2], self.dqdt[0], self.dqdt[1], self.dqdt[2])

      def XYZ(self):
           return '%s\t%12.7f\t%12.7f\t%12.7f' % (self.Symb, self.q[0], self.q[1], self.q[2])

      def  Prn(self):
           return '%s\t%4i\t%7.3f\t%12.7f\t%12.7f\t%12.7f' % (self.Symb, self.Zz, self.mass, self.q[0], self.q[1], self.q[2])

      def  Prn_XYZ(self):
           return '%s %12.7f %12.7f %12.7f' % (self.Symb, self.q[0], self.q[1], self.q[2])

      def  Prn_Vel(self):
           return '%s %12.7f %12.7f %12.7f' % (self.Symb, self.dqdt[0], self.dqdt[1], self.dqdt[2])

      def  PrnGrad(self):
           return '%15.9f\t%15.9f\t%15.9f' % (self.dedq[0], self.dedq[1], self.dedq[2])

      def  PrnGradGaussianExternal(self):
           return '%20.12E%20.12E%20.12E' % (self.dedq[0], self.dedq[1], self.dedq[2])

      def  PrnForce(self):
           return '%15.9f\t%15.9f\t%15.9f' % (-self.dedq[0], -self.dedq[1], -self.dedq[2])

      def  PrnVecs(self):
           scle=0.1
           return '%12.7f %12.7f %12.7f' % ((self.q[0]-self.dedq[0])*scle, (self.q[1]-self.dedq[1])*scle, (self.q[2]-self.dedq[2])*scle)

#------------------------------
#  Molecule class	
# Sun Dec  4 12:10:20 EST 2005 - TODO: I come up to think that the way it is done here - 
#                                Atom as a separte object really sux performance-wize
#                                A better strategy would be to store mega-arrays with 
#                                coords (x1,y1,z1,x2,y2,z2....), velocities, gradients, etc...
#                                - the way it is done in sensible progz
#                                For this I will need to write a set of the calls to set/get vaious atom properties 
#

class Molecule:

      def __init__(self):
          "Just init"
          self.atoms = []
          self.mass = 0.0;
          # the center of mass
          # coords and velocities
          self.qcm = [0.0,0.0,0.0]
          self.vcm = [0.0,0.0,0.0]
	  #
	  self.comment = "" # <joke> mozhno "pasrat' w kAmentakh", "kAmenty rulyat", anyway </joke>
	  # potential energy
	  self.energy = 0.0
	  # normalization factor for potential energy
	  self.energyNorm = 0.0
	  #
	  self.energyKin = 0.0
	  self.energyTot = 0.0
	  self.energyUnit = ""


      def __str__(self):
          "Returns the string (well, a number of strings, in fact) representaion of the molecule"
          str_ = "" 
          for curr_at in self.atoms:
               str_ += curr_at.__str__() 
	       str_ += "\n"
	  return str_     

      def setComment(self,comm_=""):
	   self.comment=comm_

      def getComment(self):
	   return self.comment

      def initHessian(self):
        "create a NumPy matrix based on the current number of atoms"
    	nAtoms = self.noOfAtoms()
    	# numpy 
    	self.Hessian = zeros((nAtoms*3,nAtoms*3)) 
        # mass-weigheted Hessian
        self.HessianMW = zeros((nAtoms*3,nAtoms*3))
    
        for ii in range(nAtoms):
          for jj in range(nAtoms):
              for kk in range(3):
                  for ll in range(3):
                      self.HessianMW[ii*3+kk][jj*3+ll] = self.HessianMW[ii*3+ll][jj*3+kk] = 1.0 / math.sqrt(self.atoms[ii].mass*self.atoms[jj].mass)  / 1822.866 #amu to au - but following MOLPRO we do not do it now (?) 

      def averageHessian(self):
        "average Hij and Hji" 
    	nAtoms = self.noOfAtoms()
        for kk in range(nAtoms*3):
            for ll in range(kk+1):
                self.Hessian[kk][ll] = self.Hessian[ll][kk] = 0.5 * (self.Hessian[kk][ll] + self.Hessian[ll][kk]) 

      def massWeightHessian(self):
        # mass-weighting coeffs are already in HessianMW
    	nAtoms = self.noOfAtoms()
        for kk in range(nAtoms*3):
            for ll in range(kk+1):
                self.HessianMW[kk][ll] = self.HessianMW[ll][kk] = self.Hessian[kk][ll] * self.HessianMW[kk][ll]

      def setEnergy(self,enn_=0.0,convToKcal=False,normalize=False):
           "sets energy, if nesessary does Hartree -> kcal/mol transformation"
	   c_ = 1.0
	   n_ = 0.0
	   if convToKcal:
	      self.energyUnit = "KCAL/MOL"
	      c_ = tokcal
	   if normalize:
	      n_ = self.energyNorm
	   self.energy = enn_ * c_ + n_

      def getEnergy(self):
           "returns energy"
	   return self.energy,self.energyUnit

      def setEnergyNorm(self,enn_=0.0,convToKcal=False):
           "sets energy normalization, if nesessary does Hartree -> kcal/mol transformation"
	   c_ = -1.0
	   if convToKcal:
	      self.energyUnit = "KCAL/MOL"
	      c_ = -1 * tokcal
	   self.energyNorm = enn_ * c_

      def normalizeEnergy(self):
          """ takes normalization factor and applies it to the current potential energy """
	  self.energy += self.energyNorm

      def getEnergyKin(self):
           "calulates, sets and returns kinetic energy"
           T = 0.0
           for atom_ in self.atoms:
	       T += atom_.getEnergyKin()
	   self.energyKin = T
	   return T

      def getEnergyTot(self):
           "calulates, sets and returns energy"
	   self.getEnergyKin()
	   self.energyTot = self.energy + self.energyKin
	   return self.energyTot

      def addAtom(self,newAtom):
           "Adds a new atom to atoms list"
	   # trying to be smart here if Python does not have type checking
	   if isinstance(newAtom,Atom):
	      self.atoms.append(newAtom)
	      self.mass += newAtom.mass
	   else:
	      print "Attempt to append non-atom to a molecule"
	      sys.exit()

      def setAtomGrad(self,atomN,gx,gy,gz,convToKcal=False,convToAng=False):        
          "sets gradients atom for the given atom"
	  if atomN > self.noOfAtoms()-1:
	     print "trying to set gradient for non-exicting atom No.", atomN
	     sys.exit()
	  self.atoms[atomN].setGrad(gx,gy,gz,convToKcal,convToAng)   
	  
      def readXYZ(self,fileName):
      
          xyz = open(fileName)

          try:
             xyz = open(fileName)
          except:
             print "molecule.readXYZ: Can't open ", fileName 
             sys.exit()

          line = xyz.readline()
          
          try:
             noAtoms = int(line)
          except:
             print "molecule.readXYZ: Can't read number of atoms from ", fileName 
             sys.exit()
          
          self.setComment(xyz.readline().strip())
          
          for ii in range(noAtoms):
              L = xyz.readline().split()
              curr_at=Atom()
              tmpXX = float(L[1])
              tmpYY = float(L[2])
              tmpZZ = float(L[3])
              # read optional velocities, if present
              if len(L) == 7:
                tmpVX = float(L[4])
                tmpVY = float(L[5])
                tmpVZ = float(L[6])
                curr_at.InitExplicit(tmpXX,tmpYY,tmpZZ,tmpVX,tmpVY,tmpVZ,Symb=L[0])
              else:
                curr_at.InitExplicit(tmpXX,tmpYY,tmpZZ,Symb=L[0])
                 
              self.addAtom(curr_at)
      
          return noAtoms

      def readGaussianExternal(self,fileName):
          """Reads the molecule from an input file created by Gaussian 'external' keyword.
          returns the number of atoms and the type of the derivs requested: 0 - energies, 
          1 - gradient, 2 - Hessian"""
          ##############QT
          try:
            gauExt = open(fileName)
          except: 
            raise #"molecule.readMoleculeGauExt: cannot open", fileName
          
          L = gauExt.readline().split()
          
          try:
             noAtoms = int(L[0])
          except:
             raise #"readMoleculeGauExt: Can't read number of atoms from ", fileName 
          
          derivativesRequested = int(L[1])    
          
          for ii in range(noAtoms):
              L = gauExt.readline().split()
              curr_at=Atom()
              tmpXX = float(L[1]) * toang
              tmpYY = float(L[2]) * toang
              tmpZZ = float(L[3]) * toang
              curr_at.InitExplicit(tmpXX,tmpYY,tmpZZ,Symb=symbol[int(L[0])])
              self.addAtom(curr_at)

          #gauExt.close()
          
          return noAtoms, derivativesRequested
      
      def getAtomGradStr(self):        
          "returns a set of strings with gradients "
          str_ = "" 
          for curr_at in self.atoms:
               str_ += curr_at.PrnGrad() 
	       str_ += "\n"
	  return str_     

      def getHessianMWStr(self, half=True):        
          "returns a set of strings with the mass-weighted Hessian "
          str_ = "" 
          if self.HessianMW is None:
	  	     return str_     
          nAtoms = self.noOfAtoms()
          for kk in range(nAtoms*3):
             if half:
                llRange = kk+1 
             else:
                llRange = nAtoms*3
             for ll in range(llRange):
                 str_ +=  "%10.7f" % self.HessianMW[kk][ll]
             str_ += "\n"
          return str_     
          
      def getHessianStr(self, half=True):        
          "returns a set of strings with the Hessian "
          str_ = "" 
          if self.Hessian is None:
	  	     return str_     
          nAtoms = self.noOfAtoms()
          for kk in range(nAtoms*3):
             if half:
                llRange = kk+1 
             else:
                llRange = nAtoms*3
             for ll in range(llRange):
                 str_ +=  "%10.7f" % self.Hessian[kk][ll]
             str_ += "\n"
          return str_     

      def getAtomForceStr(self):        
          "returns a set of strings with -1*gradients "
          str_ = "" 
          for curr_at in self.atoms:
               str_ += curr_at.PrnForce() 
	       str_ += "\n"
	  return str_     

      def noOfAtoms(self):
          "Returns subj."
	  return len(self.atoms)

      def moveAtom(self,ii,dq):
          str_ = self.atoms[ii].XYZ() + " --> "
          self.atoms[ii].Move(dq)
          str_ += self.atoms[ii].XYZ() + "\n"
	  return str_

      ##################################################################################
      def GaussianExternalOutput(self):
          "Returns the Gaussian External format output file"
          str_ = "%20.12E%20.12E%20.12E%20.12E\n" % (self.getEnergy()[0], 0.0, 0.0, 0.0 )
          for curr_at in self.atoms:
               str_ += curr_at.PrnGradGaussianExternal() 
	       str_ += "\n"
          #for curr_at in self.atoms:
          #     str_ += curr_at.XYZ() 
	      #     str_ += "\n"

          # polarizability      Polar(I), I=1,6     3D20.12
          str_ += "%20.12E%20.12E%20.12E\n" % ( 0.0, 0.0, 0.0 )
          str_ += "%20.12E%20.12E%20.12E\n" % ( 0.0, 0.0, 0.0 )

          # polarizability      Polar(I), I=1,6     3D20.12
          for i in range(self.noOfAtoms()*3):
               str_ += "%20.12E%20.12E%20.12E\n" % ( 0.0, 0.0, 0.0 )

          try:
            if self.Hessian is None:
               print "No Hess" 
               return str_.replace('E', 'D')
          except:
              # print "No Hess" 
               return str_.replace('E', 'D')

          # force constants     FFX(I), I=1,(3*NAtoms*(3*NAtoms+1))/2       3D20.12
          # the Hessian is given in lower triangular form: alpha_ij, i=1 to N, j=1 to i
          count = 0 
          for kk in range(self.noOfAtoms()*3):
             llRange = kk+1 
             for ll in range(llRange):
                 str_ += "%20.12E" % self.Hessian[kk][ll]
                 count += 1
                 if count == 3:
                    count = 0
                    str_ += "\n"
          return str_     

	  return str_.replace('E', 'D')     
      ##################################################################################

      def XYZ(self):
          "Returns a string with XYZ format representaion"
          str_ = "%s\n" % self.noOfAtoms() 
	  str_ += self.comment
          str_ += "\n"
          for curr_at in self.atoms:
               str_ += curr_at.XYZ() 
	       str_ += "\n"
	  return str_     

      def XYZArranged(self):
          """Returns a string with XYZ format representaion
	  with atoms of the same kind grouped together so MOLPRO would no screw-up anything 
	  also returns dictionary where the key is the number of an atom in rearranged XYZ
	  and the value is the real number """

          str_ = "%s\n" % self.noOfAtoms() 
	  str_ += self.comment
          str_ += "\n"

	  newVSold= {} # to get the real number for the atom using its displaced number 
	  
	  # go through atoms looking for the duplicating atom types
	  ii = 0
	  elements = {}  # dictionary Symbol:[list of atom numbers]

          for curr_at in self.atoms:
               
	       curr_symb = curr_at.Symb
	       if elements.has_key(curr_symb):
	          elements[curr_symb].append(ii)
	       else:
	          elements[curr_symb] = [ ii ]
	       
	       ii += 1

          jj = 0 # counting atoms in new arragenment
          for curr_el in elements.keys():
	       for kk in elements[curr_el]:
#		   print self.atoms[kk].XYZ()
                   str_ += self.atoms[kk].XYZ() 
                   str_ += "\n"
		   newVSold[jj] = kk 
		   jj += 1
#         print newVSold

          return str_, newVSold    

      def XYZonly(self):
          "Returns the string with XYZ format representaion"
          str_ = ""
          for curr_at in self.atoms:
               str_ += curr_at.XYZ() 
	       str_ += "\n"
	  return str_     

      def calcCM(self):
           """Calculates center of mass coords and velocity 
	      calulate once, later use self.qcm,self.vcm"""                  
           self.qcm = [0.0,0.0,0.0]
           self.vcm = [0.0,0.0,0.0]
           for curr_at in self.atoms:
               for i in range3:
           	self.qcm[i] += curr_at.q[i] * curr_at.mass  
           	self.vcm[i] += curr_at.dqdt[i] * curr_at.mass
           
           for i in range3:
               self.qcm[i] = self.qcm[i] / self.mass  
               self.vcm[i] = self.vcm[i] / self.mass  
	   return self.qcm,self.vcm    
           
      def moveToCM(self):
           "Moves the molecule to its center of mass and ajusts velocities"
	   self.calcCM()
           for curr_at in self.atoms:  
	       curr_at.Move([-foo for foo in self.qcm],[-foo for foo in self.vcm])
	   self.calcCM()

      def calcAngMom(self):
          """Calculates the angular momentum  L = r x p = r x v*m
	     Note that the coordinates and velocities are reset to the center of mass 
	     Uses the cross() function to calculte the cross product"""
	  # need Numeric module (or sci_core or whatever is to come as Numeric subsitute
          self.angMom = array([0.0,0.0,0.0])
	  # we need to move the mol to the CM as the origin
	  self.moveToCM()
          for curr_at in self.atoms:
              self.angMom = self.angMom + array(cross(curr_at.q,[curr_at.mass*vv for vv in curr_at.dqdt]))
          return self.angMom

      def calcMomOfInert(self):
           """Moments of Inertia momOfInert[3,3]
              momOfInert is also array object
	      calculate once, later use self.momOfInert"""
	   # uses kronecker function   
	   if len(self.atoms) < 2:
	      # make this more graceful later, by now abort brutally
	      print "ACHTUNG! calcMomOfInert() for monoatomic or empty molecule"
	      sys.exit()
	   #   
           self.momOfInert = resize(0.0,(3,3)) 
           for i in range3:
               for j in range3:
           	 Ipq = 0.0;
           	 for curr_at in self.atoms:
           	     Ipq = Ipq + curr_at.mass * (curr_at.Rsq() * kronecker(i,j) - curr_at.coord(i) * curr_at.coord(j))
       	         self.momOfInert[i,j] = Ipq    
           return self.momOfInert

      def calcPrincipal(self):
           """Calculates principal moments of inrtia and principal axes """
           self.moveToCM()
	   self.calcMomOfInert()
           self.Iprinc, self.princAxes = eigenvectors(self.momOfInert)
	   print  self.Iprinc
	   return self.Iprinc, self.princAxes


      def calcRotEn(self):
           """Calculates angular velocities and ROTational ENergy in kcal/mol
	      returns the latter.
	      NOTE: this moves the molecule to the center of mass and adjusts the velocities!
	      So, it is adviced to create a copy of a molecule w/copy.deepcopy (not just copy!)
	      and then call calcRotEn(), and delete the 
	      copy afterwards, unless you want the molecule to be placed at the origin"""
	   self.moveToCM()
	   self.calcMomOfInert()
	   self.calcAngMom()
           UU = inverse(self.momOfInert)
	   # omega - angular velocities
           self.omega = matrixmultiply(self.angMom,UU)
	   # Erot =  omega * L / 2 
           self.eRot = dot(self.omega,self.angMom)/2/C1
           return self.eRot


#-----------------------------

if __name__ == "__main__": 

   if len(sys.argv) < 2:
      print "Input file needed"
      sys.exit()
   
   thisMolecule = Molecule()

   xyz = open(sys.argv[1])
   
   line = xyz.readline()
   
   try:
      noAtoms = int(line)
   except:
      print "Can't read number of atoms"
      sys.exit()
   
   thisMolecule.setComment(xyz.readline().strip())
   
   for ii in range(noAtoms):
       L = xyz.readline().split()
       curr_at=Atom()
       tmpXX = float(L[1])
       tmpYY = float(L[2])
       tmpZZ = float(L[3])
       # read optional velocities, if present
       if len(L) == 7:
         tmpVX = float(L[4])
         tmpVY = float(L[5])
         tmpVZ = float(L[6])
         curr_at.InitExplicit(tmpXX,tmpYY,tmpZZ,tmpVX,tmpVY,tmpVZ,Symb=L[0])
       else:
         curr_at.InitExplicit(tmpXX,tmpYY,tmpZZ,Symb=L[0])
          
       thisMolecule.addAtom(curr_at)

   print thisMolecule
   print thisMolecule.XYZ()
   print "Grouping atoms by element:"
   XYZtxt, dict = thisMolecule.XYZArranged()
   print XYZtxt
   print dict
      
      
