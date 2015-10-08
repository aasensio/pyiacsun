# Author: cdiazbas@iac.es
# Date: 07.10.2015
# Code: Translation of IDL ftsread:
#       /usr/pkg/rsi/idl_local/data/ftsread.pro


def ftsread(ini,endi,ftsdir=None):

    """
    Extract spectral data from the Kitt Peak FTS-Spectral-Atlas
    as provided by H. Neckel, Hamburg.

    OUTPUT: Array with stepwidth 2mA . It is based on (interpolated)  
    Fourier-Transform-Spectra from the Kitt Peak Observatory 
    taken by J. Brault et al.

    CALL: [atlas,xlam] = ftsread(ini = waveIni ,endi = waveEndi)
    with a wavelength range (3290 - 12508 A).
    """

    import os
    from struct import unpack
    import numpy as np

    # Important message:
    print('Wavelength range (3290 - 12508 A)')
    
    # Atlas directory
    if ftsdir == None:
        # sdir = '/usr/pkg/rsi/idl_local/data'
        import sys
        # Path from __init__.py
        sdir = sys.path[0]
    else:
        sdir = ftsdir

    # Work directory
    tdir = os.getcwd()

    # Copy from the data a temporary range
    skip  = (ini-3290)
    count = abs(ini-endi)
    os.system('dd if='+sdir+'/fts_cent.dat of='+tdir+'/tmp bs=1000 skip='+str(skip)+' count='+str(count))
    lmbda = np.arange(ini,endi,0.002) 
    
    # Read the tmp file
    FileTmp   = open('tmp', 'rb')
    FileRead  = FileTmp.read(count*500*2)
    varDecode = unpack('!'+str(int(count*500))+'h',FileRead)
    varFinal  = np.array(varDecode)

    # Delete the tmp file
    os.remove(tdir+'/tmp')

    return [varFinal,lmbda]

if __name__ == '__main__':
    
    import matplotlib.pyplot as plt
    [atlas,xlam] = ftsread(ini = 10825,endi = 10843)

    plt.plot(xlam,atlas/1.E+4)
    plt.title('Kitt Peak FTS-Spectral-Atlas')
    plt.xlabel('Wavelength [A]')
    plt.ylim(0.,1.)
    plt.show()
    