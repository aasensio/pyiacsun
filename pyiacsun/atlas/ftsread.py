# cdiazbas@iac.es


def ftsread(ini, endi, ftsdir=None):
    """
    Extract spectral data from the interpolated disk-center 
    intensity atlas recorded at the Kitt-Peak National
    Observatory: Neckel and Labs (1984)

    Wavelength range: 3290 - 12508 A
    Wavelength step: 0.002 A

    CALL: atlas,xlam = ftsread(ini = waveIni ,endi = waveEndi)

    Args:
        ini (int): Initial wavelength
        endi (int): Final wavelength
        ftsdir (string, optional): FTS directory

    Returns:
        specFinal, lmbda: spectrum and wavelength
    """

    import os
    from struct import unpack
    import numpy as np

    # Atlas directory
    if ftsdir is None:
        ftsdir = str(__file__).split('/')
        sdir = '/'.join(ftsdir[0:-2])+'/data'
    else:
        sdir = ftsdir

    # Work directory
    tdir = os.getcwd()

    # Copy from the data a temporary range
    skip = (ini-3290)
    count = abs(ini-endi)
    os.system('dd if='+sdir+'/fts_cent.dat of='+tdir +
              '/tmp bs=1000 skip='+str(skip)+' count='+str(count))
    lmbda = np.arange(ini, endi, 0.002)

    # Read the tmp file
    FileTmp = open('tmp', 'rb')
    FileRead = FileTmp.read(count*500*2)
    varDecode = unpack('!'+str(int(count*500))+'h', FileRead)
    varFinal = np.array(varDecode)

    # Delete the tmp file
    os.remove(tdir+'/tmp')

    return varFinal/1e4, lmbda
