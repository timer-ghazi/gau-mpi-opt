#!/usr/bin/env python
# 
#
#  ------- Embarrassingly Parallel MRCC-Gaussian Geometry Optimization ------- 
#                      MRCC Energies + Gaussian Optimizer 
#                        Q.T., Marquette U., 2018--2022
#
#
# * Sept 2 2022
#   - extended to Molpro 
#  
# * Aug 7-8 2022
#   - filePrefix for generated input names
#   - Warn flag for proposeSymmConstraints()
#   - '=' sign in symmConstraints 
#   - batch processing and argument parsing 
# * Nov 25 2023
#   - Symmetric displacements to keep the symmetry of the molecule 
# 
# 

import sys
import re   
import os
from copy import deepcopy
from colored import fg, bg, attr

import molecule
from molecule import Molecule,Atom

import mrcc 
import molpro

ZERO = 0.0
X = 0
Y = 1
Z = 2
SAMESIGN = 1
NEGATIVE = -1

axes = (X, Y, Z)
axisNames = ['X','Y','Z']


##############################################################
#
# Numerical gradient evalulator class
# 
##############################################################

class GradientEvaluator:

    def __init__(self, mol_, delta=0.01*molecule.toang, symmConstraints=[], filePrefix="", codeName="MRCC", estRunScript=""):
        ''' delta is in Ang! '''

        knownCodes = ("MOLPRO", "MRCC")
        runScripts = {}
        runScripts["MOLPRO"] = "./run-molpro.sh"
        runScripts["MRCC"] = "./run-mrcc.sh"

        self.mol_ = deepcopy(mol_)
        self.delta = delta
        self.filePrefix = filePrefix
        self.codeName = codeName.upper()

        if not self.codeName in knownCodes:
            print "Unkown external code name: %s \n Known codes: %s" % (self.codeName, str(knownCodes) ) 

        if estRunScript == "":
            self.estRunScript = runScripts[self.codeName]

        self.symmConstraints = symmConstraints # list of user-supplied lines like:
                                               # H    X09 = X08    1.0151130000   Y09 = Y08   -1.7743660000   Z09 = -Z08    0.8017580

        self.variableNames = []        #  serialized variable names 'X01', 'Y01', etc. 
        self.atom2variables = {}       #  0: 'X00', 'Y00', 'Z00'

        self.symmConstraintDic = {}    #  'X02' = 'X02', 'X02' = 'ZERO', or 'X02' = '-X01' kind of thing 
        self.variableDic = {}          #  var name to ( atom#, coord#)
        self.uniqueVariables = []      #  essentially pruned self.variableNames for gradient eval 

        self.dependentVars = {} #  what vars a unique var projects into with the sign: 
                                       #  vUar : [ (uVar,sign), (sVar1,sign), ...] 
                                       #  e.g. 'X01': [ ('X01', SAMESIGN), ('X02', NEGATIVE)]
                                       #  needs to be populated after self.symmConstraintDic is populated 

        self.varDisplacementScale = {} #  varN: scaleFactor dictionary, 
                                       #  scaleFactor is usually 1.0; -1.0 for negatively dependent variables, 
                                       #  also may be used to scale symmetric displacements 

        self.var2uniqVarNoSign = {}    #  varN: uniqueVarN dictionary, 
                                       #  similar to self.symmConstraintDic, but no zero variables and no signs                                
                                       #  the sign can be found in varDisplacementScale dictionary

        self.displacementNames = [ 'MOLECULE' ]   # X00-plus, X00-mins, etc.
        self.displaced_ = {}                      # molecule instances with displaced geoms and also main geometry (i.e. 'MOLECULE') accessed via displacement names 
        self.displaced_['MOLECULE'] = self.mol_   # 

        self.unique2Displacement = {}  # "X00": "X00-plus", "X00-mins"

        self.inputFiles = []           # input files for whatever elelctronic structure code is used to evaluate energies 
        self.outputFiles = {}          # stores output file name for each displacement name (self.displacementNames), usually just displacementName + ".out"

        if self.symmConstraints == []:
            self.generateDefaultVariableNames()
        
    def generateDefaultVariableNames(self):
        nAtoms = self.mol_.noOfAtoms()
        for i in range(nAtoms):
            atV = []
            for anAxis in axes:
                coordN = axisNames[anAxis] + ("%02i" % i) 
                atV.append(coordN)
                self.variableNames.append(coordN)
                self.variableDic[coordN] = (i, anAxis)
            self.atom2variables[i] = atV
        return self.variableNames


    def proposeSymmConstraints(self, Color=True, TemplateOnly=False, Accept=False, Warn=True):
        '''Tries to figure coordinates to drop from gradient evaluation.
        Very primitive, would work pretty much for Cs symmetry only: 
        * if an atom coordinate is zero, it's assumed to be in the Cs symm plane, so the gradient would be zero
        * if an atom has a coodinate equal (or a negative) to the same axis coorinate of another atom, the coordinates are assumed to be symmetrical

        For each coordinate, outputs the serialized coordinate ID and, if a symmetry constraint detected, the linked coordinate (with or w/o '-') or ZERO, 
        or just the coordinate ID again, followed by the current coordinate value. 

        Generated plain-txt or ANSI-colored output.
        
        Use it to generate a template symmetry constraint string for human editing and later processing by parseSymmConstraints()
        If the generated constraints are OK (CHECK VERY CAREFYLLY!), they can be accepted for gradient evaluation with Accept=True

        TemplateOnly skips the symmetry detection and just returns a template for human editing. 
        '''

        nAtoms = self.mol_.noOfAtoms()

        coords_X = {}
        coords_Y = {}
        coords_Z = {}
        knownCoords = [ coords_X, coords_Y, coords_Z ]

        symmConstraintDic = {} # 'X02' = 'X02', 'X02' = 'ZERO', or 'X02' = '-X01' kind of thing 

        ##################################
        # just a template: 
        if TemplateOnly:
            outStr = ""
            for ii in range(nAtoms):

                atom = self.mol_.atoms[ii] 
                str_ = '!' + atom.Symb

                for ax in axes:
                    varN = self.variableNames[ 3*ii + ax]
                    str_ += "%6s=%6s % 15.10f" % (varN, varN, atom.coord(ax))

                str_ += "\n"
                outStr += str_
            return outStr
        #
        ##################################

        for aVar in self.variableNames:
            n, ax = self.variableDic[aVar]
            atom = self.mol_.atoms[n]
            # print aVar, n, ax, atom.coord(ax)

            if atom.coord(ax) == ZERO:
                symmConstraintDic[aVar] = 'ZERO'
                continue 

            if atom.coord(ax) in knownCoords[ax].keys():
                symmConstraintDic[aVar] = knownCoords[ax][atom.coord(ax)]
                continue 

            if (-1 * atom.coord(ax)) in knownCoords[ax].keys():
                symmConstraintDic[aVar] = '-' + knownCoords[ax][-1 * atom.coord(ax)]
                continue 
                
            knownCoords[ax][atom.coord(ax)] = aVar
            symmConstraintDic[aVar] = aVar

        outStr = ""
        outStrColor = ""

        for ii in range(nAtoms):
            atom = self.mol_.atoms[ii] 

            str_ = '!' + atom.Symb
            strC_ = attr('dim') + '!' + attr('reset') + atom.Symb

            for ax in axes:
                varN = self.variableNames[ 3*ii + ax]
                uniq = symmConstraintDic[self.variableNames[ 3*ii + ax]]

                str_  += "%6s =" % varN
                strC_ += "%6s =" % varN

                if uniq == 'ZERO':
                    strC_ += fg('dark_red_2')
                elif uniq[0] == '-':
                    strC_ += fg('blue')
                elif varN != uniq:
                    strC_ += fg('light_blue')

                str_ += "%5s % 15.10f" % (uniq, atom.coord(ax))
                strC_ += "%5s % 15.10f" % (uniq, atom.coord(ax))

                if varN != uniq:
                    strC_ += attr('reset')

            str_ += "\n"
            strC_ += "\n"

            outStr += str_
            outStrColor += strC_

        if Accept:
            if Warn:
                print "Warning: symmetry constraints generated by proposeSymmConstraints() will be used for gradient evaluation!"
                print "These are likely wrong and need to be corrected."
            self.symmConstraintDic = symmConstraintDic

        if Color:
            return outStrColor
        else:
            return outStr

    def parseSymmConstraints(self):
        '''Read symmetry constraints from self.symmConstraints lines of the same format as generated by proposeSymmConstraints().
           The following workflow is assumed: proposeSymmConstraints() is used to generate a template for editing.
           proposeSymmConstraints(TemplateOnly=True) can be used if generated constraints are all wrong.
        '''

        nAtoms = self.mol_.noOfAtoms()
        cLen = len(self.symmConstraints)

        if  self.symmConstraints == []:
            return "GradientEvaluator.parseSymmConstraints(): No self.symmConstraints defined"
        elif  cLen < nAtoms:
            # there could be empty lines 
            return "GradientEvaluator.parseSymmConstraints(): User-supplied symmetry constraint file does not match number of atoms: %i vs %i" % (cLen, nAtoms)

        validRefs = []
        validRefs.append('ZERO')

        for ii in range(nAtoms):

            L = self.symmConstraints[ii].replace("=", '').split()
            
            # variable names
            vX = L[1]
            vY = L[4]
            vZ = L[7]

            # constraints 
            cX = L[2]
            cY = L[5]
            cZ = L[8]
            
            self.variableNames.append(vX)
            self.variableNames.append(vY)
            self.variableNames.append(vZ)

            self.variableDic[vX] = ( ii, X ) 
            self.variableDic[vY] = ( ii, Y ) 
            self.variableDic[vZ] = ( ii, Z ) 

            self.atom2variables[ii] = ( vX, vY, vZ ) 

            self.symmConstraintDic[vX] = cX
            self.symmConstraintDic[vY] = cY
            self.symmConstraintDic[vZ] = cZ

            for v in ( vX, vY, vZ ):
                if v == self.symmConstraintDic[v]:
                    validRefs.append(v)
                    validRefs.append('-'+v)
                else:
                    if self.symmConstraintDic[v] not in validRefs:
                        return "GradientEvaluator.parseSymmConstraints(): %s refers to %s not previously defined (or ZERO)." % ( v, self.symmConstraintDic[v]) 
        return ""
                    
    def printSymmConstraints(self, Color=True):

        if self.symmConstraintDic == {}:
            return "" 

        nAtoms = self.mol_.noOfAtoms()
        outStr = ""
        outStrColor = ""

        for ii in range(nAtoms):
            atom = self.mol_.atoms[ii] 

            str_ = '!' + atom.Symb
            strC_ = attr('dim') + '!' + attr('reset') + atom.Symb

            for ax in axes:
                varN = self.variableNames[ 3*ii + ax]
                uniq = self.symmConstraintDic[self.variableNames[ 3*ii + ax]]

                str_  += "%6s =" % varN
                strC_ += "%6s =" % varN

                if uniq == 'ZERO':
                    strC_ += fg('dark_red_2')
                elif uniq[0] == '-':
                    strC_ += fg('blue')
                elif varN != uniq:
                    strC_ += fg('light_blue')

                str_ += "%6s % 15.10f" % (uniq, atom.coord(ax))
                strC_ += "%6s % 15.10f" % (uniq, atom.coord(ax))

                if varN != uniq:
                    strC_ += attr('reset')

            str_ += "\n"
            strC_ += "\n"

            outStr += str_
            outStrColor += strC_

        if Color:
            return outStrColor
        else:
            return outStr

    def generateDependentVarsList(self, scaleDegenrateVars=False):
        '''self.dependentVars is needed when generating displacements to symmetrically shift not only the unique varitself, 
        but also whatever variables that depend on it '''
        # 
        if self.symmConstraintDic == {}:
            print "GradientEvaluator.generateDependentVarsList: Error--empty symmConstraintDic"
            sys.exit()

        # for debug purposes, since the function doesnt need to return anything might as well do a report
        reportStr = ""
        
        # A unique variable at least affects itself, so we populate the dictionary 
        # with this default: 
        for uniqueVar in self.uniqueVariables:
            self.dependentVars[uniqueVar] = [ ( uniqueVar, SAMESIGN) ]
            reportStr += str(uniqueVar) + ": " + str(self.dependentVars[uniqueVar]) + "\n"
        

        # appending dependentVars entries with the names of the actual dependent vars
        for varN in self.variableNames:

            if self.symmConstraintDic[varN] == "ZERO":
                continue
            
            if self.symmConstraintDic[varN] == varN:
                continue

            sign = SAMESIGN  # SAMESIGN = 1, NEGATIVE = -1
            uniqueVar = self.symmConstraintDic[varN]
            dependentVarN = varN

            if uniqueVar.startswith('-'):
                sign = NEGATIVE
                uniqueVar = uniqueVar[1:]  

            self.dependentVars[uniqueVar].append( ( dependentVarN, sign) )
            reportStr += str(uniqueVar) + ": " + str(self.dependentVars[uniqueVar]) + "\n"

        # determining how many dependent variables a uniqueVar has ("degeneracy")
        # and populating varDisplacementScale and var2uniqVarNoSign dictionaries
        for uniqueVar in self.uniqueVariables:
            uniqueVarDegeneracy =  len(self.dependentVars[uniqueVar])
            
            if scaleDegenrateVars:
                scale = 1.0/uniqueVarDegeneracy
            else:
                scale = 1.0

            for dependentVarTuple in self.dependentVars[uniqueVar]:
                varN, sign = dependentVarTuple
                self.varDisplacementScale[varN] = scale * sign
                self.var2uniqVarNoSign[varN] = uniqueVar

        return reportStr
       

    def generateDisplacements(self):
        '''Creates self.displaced_ = {} with Molecule() instances accesible via X00-plus type of names'''

        if self.symmConstraintDic == {}:
            self.uniqueVariables = self.variableNames 
        else:
           for varN in self.variableNames:
                if varN == self.symmConstraintDic[varN]:
                    self.uniqueVariables.append(varN)

        # moved to init:
        # self.displaced_ = {}
        # self.displaced_['MOLECULE'] = self.mol_
        # self.displacementNames = [ 'MOLECULE' ]

        # self.generateDependentVarsList(scaleDegenrateVars=False)
        self.generateDependentVarsList(scaleDegenrateVars=True)

        Zr = 0.0 
        deltaPlus = self.delta 
        deltaMinus = self.delta * -1
            
        for uniqueVar in self.uniqueVariables:

            plus  = uniqueVar + "-plus"
            minus = uniqueVar + "-mins"

            self.unique2Displacement[uniqueVar] = ( plus, minus )

            self.displacementNames.append(plus)
            self.displacementNames.append(minus)

            self.displaced_[plus] = deepcopy(self.mol_)
            self.displaced_[minus] = deepcopy(self.mol_)

            for dependentVarTuple in self.dependentVars[uniqueVar]:

                # we dont use sign here, it is already in varDisplacementScale[VarN] 
                varN, sign = dependentVarTuple
                ii, ax = self.variableDic[varN] 
                
                Dp = deltaPlus  * self.varDisplacementScale[varN]
                Dm = deltaMinus * self.varDisplacementScale[varN]

                deltaVectorsUp   = [ [ Dp, Zr, Zr ], [ Zr, Dp, Zr], [ Zr, Zr, Dp ] ]
                deltaVectorsDown = [ [ Dm, Zr, Zr ], [ Zr, Dm, Zr], [ Zr, Zr, Dm ] ]

                self.displaced_[plus].moveAtom(ii,deltaVectorsUp[ax])
                self.displaced_[minus].moveAtom(ii,deltaVectorsDown[ax])

                #--  print "\t", varN, sign, ii, ax
                #--  print deltaVectorsUp[ax]
                #--  print deltaVectorsDown[ax]

        return len(self.displaced_)

    def showDisplacementEnergies(self):
        outS = ""

        en_, unit = self.mol_.getEnergy()
        outS += "Molecule energy:" + str(en_) + unit  + "\n" # unit in Molecule is empty by default, but just in case 

        for dispN in self.displacementNames:
            en_, unit = self.displaced_[dispN].getEnergy()
            outS += dispN + " energy:" + str(en_) + unit  + "\n" # unit in Molecule is empty by default, but just in case 
        return outS

    def calculateGradients(self):
        
        uniqGrads = {} # X00, -X00, and ZERO keys
        uniqGrads['ZERO'] = ZERO

        for varN in self.uniqueVariables:
            plus, minus =  self.unique2Displacement[varN]

            enPlus = self.displaced_[plus].getEnergy()[0]
            enMinus = self.displaced_[minus].getEnergy()[0]

            grad = (enPlus - enMinus) / (2*self.delta / molecule.toang) # delta is in Ang, so needs to be converted into Bohr

            uniqGrads[varN] = grad

            varMinus = "-" + varN
            uniqGrads[varMinus] = -1 * grad

        allGrads = {}

        if self.symmConstraintDic == {}:
            allGrads = uniqGrads
        else:
            for varN in self.variableNames:
                v = self.symmConstraintDic[varN]
                allGrads[varN] = uniqGrads[v]

        nAtoms = self.mol_.noOfAtoms()
        for ii in range(nAtoms):
            vs = self.atom2variables[ii]
            gX = allGrads[ vs[X] ]
            gY = allGrads[ vs[Y] ]
            gZ = allGrads[ vs[Z] ]
            self.mol_.setAtomGrad(ii,gX,gY,gZ)
        
    def defaultDisplacementNumber(self):
        nAtoms = self.mol_.noOfAtoms()
        return 6*nAtoms + 1

    def actualDisplacementNumber(self):
        return len(self.displacementNames)

    #
    #   General Interface stuff
    #
    def generateESTinputs(self,template):
        '''Generates input files for each displacement, populates self.inputFiles[] and self.outputFiles{}
           Execution is left to be performed outside of GradientEvaluator
        '''
        if self.codeName == "MRCC":
            self.generateMRCCinputs(template)
        elif self.codeName == "MOLPRO":
            self.generateMolproInputs(template)

        return self.inputFiles
    
    def readESToutputs(self,method=''):
        
        str_ = ""
        if self.codeName == "MRCC":
            str_ = self.readMRCCoutputs(method)
        elif self.codeName == "MOLPRO":
            str_ = self.readMolproOutputs()

        return str_

    def runEST(self):
        if self.codeName == "MRCC":
            return mrcc.executeCode
        elif self.codeName == "MOLPRO":
            return molpro.executeCode
        
    #
    #   MRCC Interface stuff
    #
    def generateMRCCinputs(self,template):
        '''Generates input files for each displacement, populates self.inputFiles[] and self.outputFiles{}
           Execution is left to be performed outside of GradientEvaluator
        '''
        print "-------------------------------------------------------"
        print "GradientEvaluator.generateMRCCinputs():"
        print "Generating MRCC MINP files using the following template:\n"
        print template,
        print "-------------------------------------------------------"
        for dispN in self.displacementNames:
            baseFlName = "%s-%s" % (self.filePrefix, dispN)
            ifn = mrcc.generateInputEnergy(self.displaced_[dispN],baseFlName=baseFlName,headerStr=template)
            self.inputFiles.append(ifn)
            self.outputFiles[dispN] = (baseFlName + ".out")

        return self.inputFiles
    
    def readMRCCoutputs(self,method):
        
        if self.outputFiles == "{}":
            print "GradientEvaluator.readMRCCoutputs(): empty output file list"
            return 

        zeroEns = [] # collect files w/errors

        for dispN in self.displacementNames:
            m_ = self.displaced_[dispN]
            o_ = self.outputFiles[dispN]
            e_ = mrcc.updateEnergyFromOutput(m_,o_,method)
            if e_ == 0.0: 
                zeroEns.append(dispN)
            # print o_ + " " + method + " energy is ", e_ 

        if zeroEns == []:
            return ""
        else:
            return "Could not read energies from: " + str(zeroEns)

    #
    #   Molpro Interface stuff
    #
    def generateMolproInputs(self,header):
        '''Generates input files for each displacement, populates self.inputFiles[] and self.outputFiles{}
           Execution is left to be performed outside of GradientEvaluator
        '''
        print "-------------------------------------------------------"
        print "GradientEvaluator.generateMolproInputs():"
        print "Generating Molpro input files using the following template:\n"
        print header,
        print "-------------------------------------------------------"
        for dispN in self.displacementNames:
            baseFlName = "%s-%s" % (self.filePrefix, dispN)
            ifn = molpro.generateInputEnergy(self.displaced_[dispN],baseFlName=baseFlName,headerStr=header)
            self.inputFiles.append(ifn)
            self.outputFiles[dispN] = (baseFlName + ".xyz")

        return self.inputFiles
    
    def readMolproOutputs(self,method=''):
        
        if self.outputFiles == "{}":
            print "GradientEvaluator.readMolproOutputs(): empty output file list"
            return 

        zeroEns = [] # collect files w/errors

        for dispN in self.displacementNames:
            m_ = self.displaced_[dispN]
            o_ = self.outputFiles[dispN].lower() # Molpro makes things lowercase 
            e_ = molpro.updateEnergyFromOutput(m_,o_)
            if e_ == 0.0: 
                zeroEns.append(dispN)
            # print o_ + " " + method + " energy is ", e_ 

        if zeroEns == []:
            return ""
        else:
            return "Could not read energies from: " + str(zeroEns)
    
##############################################################
#
# Module testing 
# 
##############################################################

MRCCinpTest = """basis=6-31G
calc=HF
geom=xyz
"""
MOLPROinpTest = """***,SinglePt
memory,3,G
gthresh,energy=1.d-8,orbital=1.d-8,grid=1.d-8

NOORIENT;
gparam,dsyevd=0

geomtyp=xyz
! the name is set by the script; do not modify:
GEOMETRY={%s
}

BASIS
default=6-31G
END

hf;
show, Energy
!  the names are set by the script; do not modify:
put,xyz,%s.xyz;
"""


if __name__ == "__main__": 

    #codeName = "MRCC" 
    codeName = "MOLPRO" 

    if codeName == "MRCC":
        inputTemplate = MRCCinpTest
    elif codeName == "MOLPRO":
        inputTemplate = MOLPROinpTest


    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument('xyzfile', metavar='xyzfile.xyz', type=str, help='XYZ file in standard orientation')
    parser.add_argument('-s', metavar='xyzfile.symm',type=str, help='Symmetry constraints file')
    parser.add_argument('--generate_inputs', action='store_true', help='Generate %s input files for displaced geometries' % codeName)
    parser.add_argument('--run', action='store_true', help='Generate and run the EST code input files for displaced geometries and calculate gradients')
    parser.add_argument('--read', action='store_true', help='Read output files for displaced geometries and calculate gradients')
    parser.add_argument('--batch', action='store_true', help='Batch processing')

    args = parser.parse_args()

    Batch = args.batch
    xyzFNinput =args.xyzfile
    symFNinput = args.s
    baseN = xyzFNinput.replace('.xyz','')
    Propose = Parse = Generate = Run = Read = False 

    if symFNinput == None:
        Propose = True
    else:
        Parse = True

    if Parse: 
        if args.generate_inputs:
            Generate = True
        elif args.run:
            Run = True
        elif args.read:
            Read = True

    if Propose:
        thisMolecule = Molecule()
        thisMolecule.readXYZ(xyzFNinput)
        
        thisGradientEv = GradientEvaluator(thisMolecule,codeName=codeName)
        
        symFN = baseN + ".symm"
        if not Batch:
            print " * Generated symmetry constraints, also dumped into %s:" % symFN
        print thisGradientEv.proposeSymmConstraints(Color=True,Accept=False)
        
        symF = open(symFN,'w')
        proposedSym = thisGradientEv.proposeSymmConstraints(Color=False,Accept=True,Warn=False)
        symF.write(proposedSym)
        symF.close()

        if Batch:
            pass 
        else:
            thisGradientEv.generateDisplacements()
            s1 = " * Energy evaluations with proposed constraints applied:" 
            s2 = attr('dim') + " * Energy evaluations without constraints applied:" 
            print "%-58s%4i" % ( s1, thisGradientEv.actualDisplacementNumber() )
            print "%-53s%4i%s" % ( s2, thisGradientEv.defaultDisplacementNumber(), attr('reset') )

    elif Parse:
        thisMolecule = Molecule()
        thisMolecule.readXYZ(xyzFNinput)

        with open(symFNinput, 'r') as f:
            userCs = f.readlines()
        
        thisGradientEv = GradientEvaluator(thisMolecule,symmConstraints=userCs,filePrefix=baseN,codeName=codeName)
        errMsg = thisGradientEv.parseSymmConstraints()
        if errMsg == "":
            thisGradientEv.generateDisplacements()
            if Batch:
                print thisGradientEv.actualDisplacementNumber()
            else:
                print " * Parsed symmetry constraints read from %s:" % symFNinput 
                print thisGradientEv.printSymmConstraints(),
                s1 = " * Energy evaluations with constraints applied:" 
                s2 = attr('dim') + " * Energy evaluations without constraints applied:" 
                print "%-53s%4i" % ( s1, thisGradientEv.actualDisplacementNumber() )
                print "%-53s%4i%s" % ( s2, thisGradientEv.defaultDisplacementNumber(), attr('reset') )
        else:
            print errMsg 

        if Generate:
            fileL = thisGradientEv.generateESTinputs(inputTemplate)
            print " * %s input files generated: %i" % ( codeName, len(fileL) )
    
        if Run:
            fileL = thisGradientEv.generateESTinputs(inputTemplate)
            print " * %s input files generated: %i" % ( codeName, len(fileL) )
            for fl in fileL:
                thisGradientEv.runEST()(fl)
            print thisGradientEv.readESToutputs("HF")

            thisGradientEv.calculateGradients()
            print " * Calculated gradients:"
            print thisGradientEv.mol_.getAtomGradStr()
    
        if Read:
            fileL = thisGradientEv.generateESTinputs(inputTemplate)
            print " * %s input files generated: %i" % ( codeName, len(fileL) )
            print thisGradientEv.readESToutputs("HF")

            thisGradientEv.calculateGradients()
            print " * Calculated gradients:"
            print thisGradientEv.mol_.getAtomGradStr()
    


