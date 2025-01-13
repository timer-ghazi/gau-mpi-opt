#!/usr/bin/env python
# 
#
#  ------- Embarrassingly Parallel MRCC-Gaussian Geometry Optimization ------- 
#                      MRCC Energies + Gaussian Optimizer 
#                        Q.T., Marquette U., 2018--2022
#

import string
import sys
import re   
import os
import molecule
from molecule import Molecule,Atom



def generateInputEnergy(mol_,baseFlName,headerStr):    
    "Generate MRCC input based on template string and the molecule object" 

    minp = headerStr % mol_.XYZ()

    estCodeInputName = baseFlName+".MINP"
    estCodeInput = open(estCodeInputName,'w')
    estCodeInput.write(minp)
    estCodeInput.close()
    
    return estCodeInputName

def executeCode(flName,estCodeScript="./run-mrcc.sh",Verbose=True):    
    "Executes Bash script" 
    if Verbose:
        print "Running MRCC script for ", flName
    errcode = os.system(estCodeScript + " " + flName)
    if Verbose:
        if errcode != 0:
            print "Error running MRCC", flName, " ", errcode
        else:
            print "Done running MRCC script for ", flName, " ", errcode
    return errcode
    #return estCodeParseOutputEnergy(outFlName=baseFlName+".out")


def parseOutputEnergy(outFlName):    
    "Extract energies; returns a dictionary" 

    energies = { }

    #with open(outFlName, 'r') as f:
    #    out = f.read()
    try:
        f = open(outFlName, 'r')
    except:
        print "ERROR: estCodeParseOutputEnergy() could not open", outFlName, ", returning empty energy dictionary"
        return energies
    out = f.read()
    f.close()

    HFstr = re.findall("(.*?FINAL HARTREE-FOCK ENERGY.*?)\n",out)
    HF = HFstr[-1].split()[-2]
    energies["HF"] = float(HF)
    postHFlist = re.findall("( Total .*?energy.*?)\n",out)

    for e in postHFlist:
        l = e.split()
        method = l[1]
        en = l[4]
        energies[ method ] = float(en)

    return energies

def updateEnergyFromOutput(mol_,outFlName,method="HF"):    
    '''If 0.0 comes back, something went wrong'''
    energies = parseOutputEnergy(outFlName)
    en_ = 0.0
    if method in energies.keys():
        en_ = energies[method]
    else:
        en_ = 0.0
        print "ERROR: updateEnergyFromOutput(): no '" + method + "' energy found in", energies, "from", fileName, "Setting energy to 0.00"

    mol_.setEnergy(enn_ = en_)

    return en_


    
##############################################################
#
# Module testing 
# 
##############################################################

minpHeader = """basis=cc-pVDZ
calc=HF
geom=xyz
"""

minpTest = """basis=6-31G
calc=ccsd(t)
mem=500MB
geom=xyz
"""

if __name__ == "__main__": 

    if len(sys.argv) < 2:
       print "Input file needed"
       sys.exit()
    
    thisMolecule = Molecule()
    thisMolecule.readXYZ(sys.argv[1])
    
    baseN = sys.argv[1].replace('.xyz','')
    
    print "\n * generateInputEnergy(): Generating %s.MINP file with a header:" % baseN
    print minpTest 
    flName = generateInputEnergy(thisMolecule,baseN,minpTest)
    
    scr = "./run-mrcc.sh" 
    print "\n * executeCode(): Running %s script for %s.MINP:" % (scr, baseN)
    # executeCode(flName,estCodeScript=scr)
    
    oF = baseN + ".out" 
    print "\n * parseOutputEnergy(): Extracting 'method':energy dictionary from %s:" % (oF)
    enD = parseOutputEnergy(oF)
    for mth in enD.keys():
        print "%10s:\t% 12.9f" % (mth, enD[mth])
    
    mth = "CCSD"
    print "\n * updateEnergyFromOutput(): Updating a Molecule() instance energy with %s energy from %s:" % (mth, oF)
    updateEnergyFromOutput(thisMolecule,oF,method=mth)
    print thisMolecule.getEnergy()

    
    # for fl in fileL:
    #      estCodeExecute(fl)
    # print thisGradientEv.readMRCCoutputs("HF")
    # print thisGradientEv.showDisplacementEnergies()
    # thisGradientEv.calculateGradients()
    


