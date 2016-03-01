from __future__ import print_function
import numpy as np
from .sir import *
import os.path
import shutil

__all__ = ["listLinesSIR", "initializeSIR", "synthesizeSIR"]

def listLinesSIR():
    """List the lines available in SIR for synthesis
        
    """
    if (not os.path.exists('LINEAS')):
        local = str(__file__).split('/')
        sdir = '/'.join(local[0:-2])+'/data'
        shutil.copy(sdir+'/LINEAS', os.getcwd())

    f = open('LINEAS', 'r')
    lines = f.readlines()
    f.close()

    print("Available lines:")
    for l in lines[:-1]:
        print(l)

def initializeSIR(lines):
    """Initialize the SIR synthesis code for a set of spectral lines
    
    Args:
        lines (list): list of lists containing the information for the lines to be synthesized
                    Each region is defined is defined by a list containing the following elements:
                     - A string defining which lines are synthesized in the region. E.g. '1,2,3'
                     - Initial wavelength displacement in mA
                     - Step in mA
                     - Final wavelength displacement in mA

                     E.g. lines = [['1', -500.0, 10.0, 1500.0], ['2', -750.0, 10.0, 1300.0]]
    
    Returns:
        int: the number of wavelength points to be synthesized

    """
    f = open('malla.grid', 'w')
    f.write("IMPORTANT: a) All items must be separated by commas.                 \n")
    f.write("           b) The first six characters of the last line                \n")
    f.write("          in the header (if any) must contain the symbol ---       \n")
    f.write("\n")                                                                       
    f.write("Line and blends indices   :   Initial lambda     Step     Final lambda \n")
    f.write("(in this order)                    (mA)          (mA)         (mA)     \n")
    f.write("-----------------------------------------------------------------------\n")

    for i in range(len(lines)):
        f.write("{0}    :  {1}, {2}, {3}".format(lines[i][0], lines[i][1], lines[i][2], lines[i][3]))
    f.close()

    if (not os.path.exists('LINEAS')):
        local = str(__file__).split('/')
        sdir = '/'.join(local[0:-2])+'/data'
        shutil.copy(sdir+'/LINEAS', os.getcwd())

    if (not os.path.exists('THEVENIN')):
        local = str(__file__).split('/')
        sdir = '/'.join(local[0:-2])+'/data'
        shutil.copy(sdir+'/THEVENIN', os.getcwd())
        
    return init()    

def synthesizeSIR(model, macroturbulence=0.0, fillingFactor=1.0, stray=0.0):
    """Carry out the synthesis and returns the Stokes parameters and the response 
    functions to all physical variables at all depths
    
    Args:
        model (float array): an array of size [nDepth x 8], where the columns contain the depth stratification of:
            - log tau
            - Temperature [K]
            - Electron pressure [dyn cm^-2]
            - Microturbulent velocity [km/s]
            - Magnetic field strength [G]
            - Line-of-sight velocity [km/s]
            - Magnetic field inclination [deg]
            - Magnetic field azimuth [deg]
        macroturbulence (float, optional): macroturbulence velocity [km/s]. Default: 0
        fillingFactor (float, optional): filling factor. Default: 1
        stray (float, optional): stray light in %. Default: 0
    
    Returns:
        stokes: (float array) Stokes parameters, with the first index containing the wavelength displacement and the remaining
                                containing I, Q, U and V. Size (5,nLambda)
        rf: (float array) Response functions to T, Pe, vmic, B, v, theta, phi. Size (4,nLambda,nDepth)
    """
    stokes, rf = sir.synth(model, macroturbulence, fillingFactor, stray)

    return stokes, rf

if (__name__ == "__main__"):
    listLinesSIR()
    l = [['1',-500.,10.,1500.]]
    nLambda = initializeSIR(l)
    out = np.loadtxt('../test/modelo.mod', dtype=np.float32, skiprows=1)[:,0:8]
    stokes, rf = synthesizeSIR(out)