# cdiazbas@iac.es


def Delbouille73(ini, endi, atlasdir=None):
    """
    Extract spectral data from the disk-center intensity atlas:
    Delbouille L., Neven L., Roland G. (1973)

    Wavelength range: 3000 - 10.000 A 
    Wavelength step (visible): 0.002 A

    Downloaded from:
    http://bass2000.obspm.fr/solar_spect.php

    Args:
        ini (int): Initial wavelength
        endi (int): Final wavelength
        atlasdir (string, optional): Atlas directory
    """

    import numpy as np

    # Atlas directory
    if atlasdir is None:
        atlasdir = str(__file__).split('/')
        sdir = '/'.join(atlasdir[0:-2])+'/data'
    else:
        sdir = atlasdir

    file0 = np.load(sdir + '/Delbouille73.npy')
    lmbda0 = np.arange(3000., 10000., 0.002)
    iniI = np.argmin(abs(ini - lmbda0))
    endiI = np.argmin(abs(endi - lmbda0))

    lmbda = lmbda0[iniI:endiI]
    varFinal = file0[iniI:endiI]

    return [varFinal / 1e4, lmbda]

if __name__ == '__main__':

    import matplotlib.pyplot as plt
    [atlas, xlam] = Delbouille73(ini=6300, endi=6303)

    plt.plot(xlam, atlas / max(atlas))
    plt.title('Delbouille73 Atlas')
    plt.xlabel('Wavelength [A]')
    plt.ylim(0., 1.)
    plt.show()
