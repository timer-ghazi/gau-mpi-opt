#!/usr/bin/env python
# 
#
#  ------- Embarrassingly Parallel MRCC-Gaussian Geometry Optimization ------- 
#                   MRCC/MOLPRO Energies + Gaussian Optimizer 
#                        Q.T., Marquette U., 2018--2022
#

import string
import sys
import re   
import os
import molecule
from molecule import Molecule,Atom



def generateInputEnergy(mol_,baseFlName,headerStr):    
    "Generate MOLPRO input based on template string and the molecule object" 

    minp = headerStr % ( mol_.XYZ(), baseFlName )

    estCodeInputName = baseFlName+".mpi"
    estCodeInput = open(estCodeInputName,'w')
    estCodeInput.write(minp)
    estCodeInput.close()
    
    return estCodeInputName

def executeCode(flName,estCodeScript="./run-molpro.sh",Verbose=True):    
    "Executes Bash script" 
    if Verbose:
        print "Running Molpro script for ", flName
    errcode = os.system(estCodeScript + " " + flName)
    if Verbose:
        if errcode != 0:
            print "Error running ", flName, " ", errcode
        else:
            print "Done running Molpro script for ", flName, " ", errcode
    return errcode
    #  return parseOutputEnergy(outFlName=baseFlName+".xyz")


def parseOutputEnergy(outFlName):    
    "Extract energy from the XYZ file dumped by Molpro; keeping dictionary for compatibility" 

    energies = { }

    #with open(outFlName, 'r') as f:
    #    out = f.read()
    try:
        f = open(outFlName, 'r')
    except:
        print "ERROR: parseOutputEnergy() could not open", outFlName, ", returning empty energy dictionary"
        return energies
    out = f.readlines()
    f.close()

    EstrL = out[1].split()
    Method = EstrL[0]
    energies[Method] = float(EstrL[2])

    return energies

def updateEnergyFromOutput(mol_,outFlName,method=""):    
    '''If 0.0 comes back, something went wrong; keeping method key for compatibility'''
    energies = parseOutputEnergy(outFlName)
    en_ = 0.0

    method = energies.keys()[0]
    en_ = energies[method]

    mol_.setEnergy(enn_ = en_)

    return en_


    
##############################################################
#
# Module testing 
# 
##############################################################


minpTest = """***,SinglePt
memory,3,G
gthresh,energy=1.d-8,orbital=1.d-8,grid=1.d-8

NOORIENT;
gparam,dsyevd=0

geomtyp=xyz
! the name is set by the script; do not modify:
GEOMETRY={%s
}

BASIS
default=vdz
END

hf;
show, Energy
!  the names are set by the script; do not modify:
put,xyz,%s.xyz;
"""

if __name__ == "__main__": 

    if len(sys.argv) < 2:
       print "Input file needed"
       sys.exit()
    
    thisMolecule = Molecule()
    thisMolecule.readXYZ(sys.argv[1])
    
    baseN = sys.argv[1].replace('.xyz','') + "-mlp" # so Molpro doesn't overwrite the original XYZ
    
    print "\n * generateInputEnergy(): Generating %s.mpi file with a template:" % baseN
    print minpTest 
    flName = generateInputEnergy(thisMolecule,baseN,minpTest)
    
    scr = "./run-molpro.sh" 
    print "\n * executeCode(): Running %s script for %s.mpi:" % (scr, baseN)
    #executeCode(flName,estCodeScript=scr)
    
    oF = baseN + ".xyz" 
    print "\n * parseOutputEnergy(): Extracting 'method':energy dictionary from %s:" % (oF)
    enD = parseOutputEnergy(oF)
    for mth in enD.keys():
        print "%10s:\t% 12.9f" % (mth, enD[mth])
    
    print "\n * updateEnergyFromOutput(): Updating a Molecule() instance energy with %s energy from %s:" % (mth, oF)
    updateEnergyFromOutput(thisMolecule,oF)
    print thisMolecule.getEnergy()

    
    # for fl in fileL:
    #      estCodeExecute(fl)
    # print thisGradientEv.readMRCCoutputs("HF")
    # print thisGradientEv.showDisplacementEnergies()
    # thisGradientEv.calculateGradients()
    


