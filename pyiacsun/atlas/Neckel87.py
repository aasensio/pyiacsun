# cdiazbas@iac.es


def Neckel87(ini, endi, atlasdir=None):
    """
    Extract spectral data from the disk-center intensity atlas:
    Brault & Neckel (1987)

    Wavelength range: 3290 - 12510 A
    Wavelength step: from 0.004 A to 0.02

    Downloaded from:
    ftp://ftp.hs.uni-hamburg.de/pub/outgoing/FTS-Atlas

    Args:
        ini (int): Initial wavelength
        endi (int): Final wavelength
        ftsdir (string, optional): Atlas directory
    """

    import numpy as np

    # Atlas directory
    if atlasdir is None:
        atlasdir = str(__file__).split('/')
        sdir = '/'.join(atlasdir[0:-2])+'/data'
    else:
        sdir = atlasdir

    file0 = np.load(sdir + '/Neckel87.npy')
    iniI = np.argmin(abs(ini - file0[:, 0]))
    endiI = np.argmin(abs(endi - file0[:, 0]))

    lmbda = file0[iniI:endiI, 0]
    varFinal = file0[iniI:endiI, 1]

    return [varFinal, lmbda]

if __name__ == '__main__':

    import matplotlib.pyplot as plt
    [atlas, xlam] = Neckel87(ini=6300, endi=6303)

    plt.plot(xlam, atlas / max(atlas))
    plt.title('Neckel87 Atlas')
    plt.xlabel('Wavelength [A]')
    plt.ylim(0., 1.)
    plt.show()
