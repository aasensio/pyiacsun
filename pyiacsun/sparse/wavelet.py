# Aims to provide functions for fast periodic daubechies transforms (forward and inverse) in 2D
# https://github.com/amitgroup/amitgroup/tree/master/amitgroup/util/wavelet

import numpy as np
import scipy
import scipy.sparse

SPARSITY_THRESHOLD = 256

def _populate(W, filtr, yoffset):
    N = len(filtr)
    for i in range(W.shape[1]//2):
        for j in range(N):
            W[yoffset+i, (-(N-2)//2+2*i+j)%W.shape[1]] += filtr[j]

def _create_W(shape, level, filter_low, filter_high):
    d = 1 << (level-1)
    sh = (shape[0]//d, shape[1]//d)
    if np.min(sh) >= SPARSITY_THRESHOLD:
        W = scipy.sparse.lil_matrix(sh)
    else:
        W = np.asmatrix(np.zeros(sh))

    _populate(W, filter_low, 0)
    _populate(W, filter_high, shape[0]//(2*d))
    if scipy.sparse.issparse(W):
        return W.tocsr() # Faster if you're not changing it anymore
    else:
        return W

def _create_single(shape, level, filtr):
    d = 1 << (level-1)
    sh = (shape[0]//(2*d), shape[1]//d)
    if np.min(sh) >= SPARSITY_THRESHOLD:
        GH = scipy.sparse.lil_matrix(sh)
    else:
        GH = np.asmatrix(np.zeros(sh))
    _populate(GH, filtr, 0)
    if scipy.sparse.issparse(GH):
        return GH.tocsr() # Faster if you're not changing it anymore
    else:
        return GH 

def _qdot(X, A):
    return X * A * X.T

_db_kernels = {
    'haar': np.array([0.70710678118654757, 0.70710678118654757]),
    'db1': np.array([0.70710678118654757, 0.70710678118654757]),
    'db2': np.array([0.48296291314469025, 0.83651630373746899, 0.22414386804185735, -0.12940952255092145]), # D4
    'db3': np.array([0.33267055295095688, 0.80689150931333875, 0.45987750211933132, -0.13501102001039084, -0.085441273882241486, 0.035226291882100656]),
    'db4': np.array([0.23037781330885523, 0.71484657055254153, 0.63088076792959036, -0.027983769416983849, -0.18703481171888114, 0.030841381835986965, 0.032883011666982945, -0.010597401784997278]),
    'db5': np.array([0.16010239797412501, 0.60382926979747287, 0.72430852843857441, 0.13842814590110342, -0.24229488706619015, -0.03224486958502952, 0.077571493840065148, -0.0062414902130117052, -0.012580751999015526, 0.0033357252850015492]),
    'db6': np.array([0.11154074335008017, 0.49462389039838539, 0.75113390802157753, 0.3152503517092432, -0.22626469396516913, -0.12976686756709563, 0.097501605587079362, 0.027522865530016288, -0.031582039318031156, 0.0005538422009938016, 0.0047772575110106514, -0.0010773010849955799]),
    'db7': np.array([0.077852054085062364, 0.39653931948230575, 0.72913209084655506, 0.4697822874053586, -0.14390600392910627, -0.22403618499416572, 0.071309219267050042, 0.080612609151065898, -0.038029936935034633, -0.01657454163101562, 0.012550998556013784, 0.00042957797300470274, -0.0018016407039998328, 0.00035371380000103988]),
    'db8': np.array([0.054415842243081609, 0.31287159091446592, 0.67563073629801285, 0.58535468365486909, -0.015829105256023893, -0.28401554296242809, 0.00047248457399797254, 0.12874742662018601, -0.017369301002022108, -0.044088253931064719, 0.013981027917015516, 0.0087460940470156547, -0.0048703529930106603, -0.00039174037299597711, 0.00067544940599855677, -0.00011747678400228192]),
    'db9': np.array([0.038077947363167282, 0.24383467463766728, 0.6048231236767786, 0.65728807803663891, 0.13319738582208895, -0.29327378327258685, -0.096840783220879037, 0.14854074933476008, 0.030725681478322865, -0.067632829059523988, 0.00025094711499193845, 0.022361662123515244, -0.004723204757894831, -0.0042815036819047227, 0.0018476468829611268, 0.00023038576399541288, -0.00025196318899817888, 3.9347319995026124e-005]),
    'db10': np.array([0.026670057900950818, 0.18817680007762133, 0.52720118893091983, 0.68845903945259213, 0.28117234366042648, -0.24984642432648865, -0.19594627437659665, 0.12736934033574265, 0.093057364603806592, -0.071394147165860775, -0.029457536821945671, 0.033212674058933238, 0.0036065535669883944, -0.010733175482979604, 0.0013953517469940798, 0.0019924052949908499, -0.00068585669500468248, -0.0001164668549943862, 9.3588670001089845e-005, -1.3264203002354869e-005]),
    'db11': np.array([0.018694297761470441, 0.14406702115061959, 0.44989976435603013, 0.68568677491617847, 0.41196436894789695, -0.16227524502747828, -0.27423084681792875, 0.066043588196690886, 0.14981201246638268, -0.04647995511667613, -0.066438785695020222, 0.031335090219045313, 0.020840904360180039, -0.015364820906201324, -0.0033408588730145018, 0.0049284176560587777, -0.00030859285881515924, -0.00089302325066623663, 0.00024915252355281426, 5.4439074699366381e-005, -3.4634984186983789e-005, 4.4942742772363519e-006]),
    'db12': np.array([0.013112257957229239, 0.10956627282118277, 0.37735513521420411, 0.65719872257929113, 0.51588647842780067, -0.044763885653777619, -0.31617845375277914, -0.023779257256064865, 0.18247860592758275, 0.0053595696743599965, -0.09643212009649671, 0.010849130255828966, 0.041546277495087637, -0.01221864906974642, -0.012840825198299882, 0.0067114990087955486, 0.0022486072409952287, -0.0021795036186277044, 6.5451282125215034e-006, 0.00038865306282092672, -8.8504109208203182e-005, -2.4241545757030318e-005, 1.2776952219379579e-005, -1.5290717580684923e-006]),
    'db13': np.array([0.0092021335389622788, 0.082861243872901946, 0.31199632216043488, 0.61105585115878114, 0.58888957043121193, 0.086985726179645007, -0.31497290771138414, -0.12457673075080665, 0.17947607942935084, 0.072948933656788742, -0.10580761818792761, -0.026488406475345658, 0.056139477100276156, 0.0023799722540522269, -0.023831420710327809, 0.0039239414487955773, 0.0072555894016171187, -0.002761911234656831, -0.0013156739118922766, 0.00093232613086724904, 4.9251525126285676e-005, -0.00016512898855650571, 3.0678537579324358e-005, 1.0441930571407941e-005, -4.7004164793608082e-006, 5.2200350984547998e-007]),
    'db14': np.array([0.0064611534600864905, 0.062364758849384874, 0.25485026779256437, 0.55430561794077093, 0.63118784910471981, 0.21867068775886594, -0.27168855227867705, -0.21803352999321651, 0.13839521386479153, 0.13998901658445695, -0.086748411568110598, -0.071548955503983505, 0.05523712625925082, 0.026981408307947971, -0.030185351540353976, -0.0056150495303375755, 0.012789493266340071, -0.00074621898926387534, -0.003849638868019787, 0.001061691085606874, 0.00070802115423540481, -0.00038683194731287514, -4.1777245770370672e-005, 6.875504252695734e-005, -1.0337209184568496e-005, -4.3897049017804176e-006, 1.7249946753674012e-006, -1.7871399683109222e-007]),
    'db15': np.array([0.0045385373615773762, 0.046743394892750617, 0.20602386398692688, 0.49263177170797529, 0.64581314035721027, 0.33900253545462167, -0.19320413960907623, -0.28888259656686216, 0.065282952848765688, 0.19014671400708816, -0.039666176555733602, -0.11112093603713753, 0.033877143923563204, 0.054780550584559995, -0.025767007328366939, -0.020810050169636805, 0.015083918027862582, 0.0051010003604228726, -0.0064877345603061454, -0.00024175649075894543, 0.0019433239803823459, -0.00037348235413726472, -0.00035956524436229364, 0.00015589648992055726, 2.579269915531323e-005, -2.8133296266037558e-005, 3.3629871817363823e-006, 1.8112704079399406e-006, -6.3168823258794506e-007, 6.1333599133037138e-008]),
    'db16': np.array([0.0031892209253436892, 0.034907714323629047, 0.1650642834886438, 0.43031272284545874, 0.63735633208298326, 0.44029025688580486, -0.089751089402363524, -0.32706331052747578, -0.02791820813292813, 0.21119069394696974, 0.027340263752899923, -0.13238830556335474, -0.0062397227521562536, 0.075924236044457791, -0.0075889743686425939, -0.036888397691556774, 0.010297659641009963, 0.013993768859843242, -0.0069900145633907508, -0.0036442796214883506, 0.00312802338120381, 0.00040789698084934395, -0.00094102174935854332, 0.00011424152003843815, 0.00017478724522506327, -6.103596621404321e-005, -1.394566898819319e-005, 1.133660866126152e-005, -1.0435713423102517e-006, -7.3636567854418147e-007, 2.3087840868545578e-007, -2.1093396300980412e-008]),
    'db17': np.array([0.0022418070010387899, 0.025985393703623173, 0.13121490330791097, 0.37035072415288578, 0.61099661568502728, 0.5183157640572823, 0.027314970403312946, -0.32832074836418546, -0.12659975221599248, 0.19731058956508457, 0.10113548917744287, -0.12681569177849797, -0.057091419631858077, 0.081105986654080822, 0.022312336178011833, -0.046922438389378908, -0.0032709555358783646, 0.022733676583919053, -0.0030429899813869555, -0.0086029215203478147, 0.0029679966915180638, 0.0023012052421511474, -0.001436845304805, -0.00032813251941022427, 0.00043946542776894542, -2.5610109566546042e-005, -8.2048032024582121e-005, 2.3186813798761639e-005, 6.9906009850812941e-006, -4.5059424772259631e-006, 3.0165496099963414e-007, 2.9577009333187617e-007, -8.4239484460081536e-008, 7.2674929685663697e-009]),
    'db18': np.array([0.0015763102184365595, 0.019288531724094969, 0.10358846582214751, 0.31467894133619284, 0.57182680776508177, 0.57180165488712198, 0.14722311196952223, -0.29365404073579809, -0.21648093400458224, 0.14953397556500755, 0.16708131276294505, -0.092331884150304119, -0.10675224665906288, 0.064887216212358198, 0.057051247739058272, -0.04452614190225633, -0.023733210395336858, 0.026670705926689853, 0.0062621679544386608, -0.013051480946517112, 0.00011863003387493042, 0.0049433436054565939, -0.0011187326669886426, -0.0013405962983313922, 0.00062846568296447147, 0.0002135815619103188, -0.00019864855231101547, -1.5359171230213409e-007, 3.7412378807308472e-005, -8.5206025374234635e-006, -3.3326344788769603e-006, 1.7687129836228861e-006, -7.691632689865049e-008, -1.1760987670250871e-007, 3.0688358630370302e-008, -2.5079344549419292e-009]),
    'db19': np.array([0.0011086697631864314, 0.01428109845082521, 0.08127811326580564, 0.26438843174202237, 0.52443637746688621, 0.60170454913009164, 0.26089495265212009, -0.22809139421653665, -0.28583863175723145, 0.074652269708066474, 0.21234974330662043, -0.033518541903202262, -0.14278569504021468, 0.027584350624887129, 0.086906755555450702, -0.026501236250778635, -0.045674226277784918, 0.021623767409452484, 0.019375549889114482, -0.013988388678695632, -0.0058669222811121953, 0.0070407473670804953, 0.00076895435922424884, -0.0026875518007344408, 0.00034180865344939543, 0.0007358025205041731, -0.00026067613568119951, -0.00012460079173506306, 8.7112704672504432e-005, 5.1059504870906939e-006, -1.6640176297224622e-005, 3.0109643163099385e-006, 1.5319314766978769e-006, -6.8627556577981102e-007, 1.4470882988040879e-008, 4.6369377758023682e-008, -1.1164020670405678e-008, 8.6668488390344833e-010]),
    'db20': np.array([0.00077995361366591117, 0.010549394624937735, 0.063423780459005291, 0.21994211355113222, 0.47269618531033147, 0.61049323893785579, 0.36150229873889705, -0.13921208801128787, -0.32678680043353758, -0.016727088308801888, 0.22829105082013823, 0.039850246458519104, -0.15545875070604531, -0.024716827337521424, 0.10229171917513397, 0.0056322468576854544, -0.061722899624668884, 0.0058746818113949465, 0.032294299530119162, -0.0087893249245557647, -0.013810526137727442, 0.0067216273018096935, 0.0044205423867663502, -0.003581494259744107, -0.00083156217287724745, 0.0013925596193045254, -5.3497598443404532e-005, -0.0003851047486990061, 0.00010153288973669777, 6.7742808283730477e-005, -3.7105861833906152e-005, -4.3761438621821972e-006, 7.2412482876637907e-006, -1.0119940100181473e-006, -6.847079596993149e-007, 2.633924226266962e-007, 2.0143220235374613e-010, -1.8148432482976221e-008, 4.05612705554717e-009, -2.9988364896157532e-010]),
}

def _get_filters(wavelet):
    global _db_kernels
    try:
        lowpass = _db_kernels[wavelet]
    except KeyError:
        raise ValueError("Wavelet type not supported: ('{0}')".format(wavelet))
    highpass = lowpass[::-1].copy()
    highpass[1::2] *= -1
    return lowpass, highpass

def _arrange_filter_matrices(shape, wavelet):
    filter_low, filter_high = _get_filters(wavelet)

    assert len(shape) == 2
    assert shape[0] == shape[1], "Shape must be square (at least for now)"

    levels_list = range(1, int(np.log2(shape[0]))+1)

    max_level = int(np.log2(max(shape)))

    # Setup matrices
    Ws = [_create_W(shape, level, filter_low, filter_high) for level in levels_list]
    Wgs = []

    # Combine all the matrices for the steps where we throw away the coefficients.
    Wg = np.asmatrix(np.eye(shape[0], shape[1]))
    for l in range(0, max_level):
        new_M = Ws[l] * Wg
        if np.min(new_M.shape) >= SPARSITY_THRESHOLD:
            new_M = scipy.sparse.csr_matrix(new_M)
        Wgs.append(new_M)
        Wg = _create_single(shape, l+1, filter_low) * Wg
    Wgs.append(Wg)

    Wgs = Wgs[::-1]

    WsT = [W.T for W in Ws]
    WgsT = [Wg.T for Wg in Wgs]
    return Wgs, WgsT, Ws, WsT,  max_level

def daubechies_factory(shape, wavelet='db2'):
    """
    Creates a forward and an inverse discrete wavelet transform function.

    The function is specialized for a specific size and wavelet.

    .. seealso::

        :ref:`wavelet2d`

    Parameters
    ----------
    shape : tuple 
        A tuple describing the size of the input, for instance ``(32, 32)``. Values must be powers of two.
    wavelet : str
        Type of wavelet described as a string. Supported values are ``'db1'``, ``'db2'``, ... ``'db20'``. What is called, for instance, `db3` is what is generally called the D6 wavelet (since it uses a kernel of size 6). The string ``'haar'`` is a synonym for ``'db1'``.

    Returns
    -------
    wavedec2 : func(A[, levels])
        Returns a function that takes an argument, `A`, input data that must be of the size specified above. It also takes an optional argument `levels`, where you can specify how many coefficient levels you plan to use. It will return an array with the coefficients of shape ``(2**levels, 2**levels)``.
    waverec2 : func(coefs)
        Returns a function that takes a single argument, `coefs`, the coefficients to use to reconstruct the spatial information. 

    Examples
    --------
    >>> import amitgroup as ag
    >>> import amitgroup.util.wavelet
    >>> import matplotlib.pylab as plt
    >>> face = ag.io.load_example('faces')[0]

    To compress a face and then inspect the results, let's first create the transform functions:

    >>> wavedec2, waverec2 = ag.util.wavelet.daubechies_factory(face.shape, 'db8')

    And then deconstruct a face to coefficients and the reconstruct it again. Since we only used 4 coefficient levels, information will be lost. 
    
    >>> new_face = waverec2(wavedec2(face, levels=4))
    >>> ag.plot.images([face, new_face])
    >>> plt.show()
    """
    if isinstance(shape, int): # One dimensional!
        Wgs, WgsT, Ws, WsT, max_level = _arrange_filter_matrices((shape, shape), wavelet) 
        def wavedec(A, levels=np.inf):
            A = A.reshape((len(A), 1))
            levels = min(max_level, levels)
            coefs = Wgs[levels] * A 
            for l in range(levels-1, 0, -1):
                N = 1 << l 
                coefs[:N] = Ws[max_level-l] * coefs[:N]
            return np.asarray(coefs).flatten()
        def waverec(coefs):
            levels = int(np.log2(coefs.shape[0]))
            A = coefs.reshape((len(coefs), 1)).copy()
            for l in range(1, levels):
                N = 1 << l 
                A[:N] = WsT[max_level-l] * A[:N]
            A = WgsT[levels] * A
            return np.asarray(A).flatten()
        return wavedec, waverec 

    elif isinstance(shape, tuple) and len(shape) == 2: # 2 dimensional!
        Wgs, WgsT, Ws, WsT, max_level = _arrange_filter_matrices(shape, wavelet) 
        def wavedec2(A, levels=np.inf):
            levels = min(max_level, levels)
            coefs = Wgs[levels] * A * WgsT[levels]
            for l in range(levels-1, 0, -1):
                N = 1 << l 
                L = max_level-l
                coefs[:N,:N] = Ws[L] * coefs[:N,:N] * WsT[L]
            return np.asarray(coefs)
        def waverec2(coefs, levels=np.inf):
            #print coefs.shape
            levels = int(np.log2(coefs.shape[0]))
            #levels = min(max_level, levels)
            A = coefs.copy()
            for l in range(1, levels):
                N = 1 << l 
                L = max_level-l
                A[:N,:N] = WsT[L] * A[:N,:N] * Ws[L]
            return np.asarray(WgsT[levels] * A * Wgs[levels])
        return wavedec2, waverec2 
    else:
        raise ValueError("Shape must be either integer or tuple of size two")


# CACHED 1-D
################################################################################

_db_wavedec_cache = {}
_db_waverec_cache = {}

def wavedec(A, wavelet='db2', levels=np.inf, length=None):
    """
    Performs a 1D wavelet decomposition (forward transform).

    .. note::

        This function runs :func:`daubechies_factory` for you and caches the value. This means that first time you call it,
        performance will be slower than expected. You will also incur a dictionary lookup, which might cost you 100 ns. 

    .. seealso::

        :ref:`wavelet1d`

    Parameters
    ----------
    A : ndarray
        1D input data. Length must be powers of two.
    wavelet : str
        Wavelet type. See :func:`daubechies_factory`.
    levels : int
        Specify how many levels of coefficients you plan to use. The default is ``np.inf``, which will default to the maximum number possible, which will make the coefficient array the same length as `A`. Notice that `levels` is zero-based, in the sense that entering 0 is valid and will the transform operate only on the energy-level coefficient.
    """
    global _db_wavedec_cache, _db_waverec_cache
    tup = (length or len(A), wavelet)
    try:
        dec = _db_wavedec_cache[tup]
    except KeyError:
        dec, rec = daubechies_factory(*tup)
        _db_wavedec_cache[tup] = dec
        _db_waverec_cache[tup] = rec 

    return dec(A, levels) 

def waverec(coefs, wavelet='db2', length=None):
    """
    Performs a 1D wavelet reconstruction (inverse transform).

    In :func:`wavedec`, you specify `levels`, which is not done in this function since it can be inferred from the shape of `coefs`.
    
    .. note::

        This function runs :func:`daubechies_factory` for you and caches the value. This means that first time you call it,
        performance will be slower than expected. You will also incur a dictionary lookup, which might cost you 100 ns. 

    .. seealso::

        :ref:`wavelet1d`

    Parameters
    ----------
    A : ndarray
        1D input data. Length must be powers of two and square.
    wavelet : str
        Wavelet type. See :func:`daubechies_factory`.
    """
    global _db_wavedec_cache, _db_waverec_cache
    tup = (length or len(coefs), wavelet)
    try:
        rec = _db_waverec_cache[tup]
    except KeyError:
        dec, rec = daubechies_factory(*tup)
        _db_wavedec_cache[tup] = dec
        _db_waverec_cache[tup] = rec

    return rec(coefs) 

# CACHED 2-D
################################################################################

_db_wavedec2_cache = {}
_db_waverec2_cache = {}

def wavedec2(A, wavelet='db2', levels=np.inf, shape=None):
    """
    Performs a 2D wavelet decomposition (forward transform).

    .. note::

        This function runs :func:`daubechies_factory` for you and caches the value. This means that first time you call it,
        performance will be slower than expected. You will also incur a dictionary lookup, which might cost you 100 ns. 

    .. seealso::

        :ref:`wavelet2d`

    Parameters
    ----------
    A : ndarray
        2D input data. Shape must be powers of two and square.
    wavelet : str
        Wavelet type. See :func:`daubechies_factory`.
    levels : int
        Specify how many levels of coefficients you plan to use. The default is ``np.inf``, which will default to the maximum number possible, which will make the coefficient array the same size as `A`. Notice that `levels` is zero-based, in the sense that entering 0 is valid and will the transform operate only on the energy-level coefficient.
    """
    global _db_wavedec2_cache, _db_waverec2_cache
    tup = (shape or A.shape, wavelet)
    try:
        dec = _db_wavedec2_cache[tup]
    except KeyError:
        dec, rec = daubechies_factory(*tup) 
        _db_wavedec2_cache[tup] = dec
        _db_waverec2_cache[tup] = rec

    return dec(A, levels) 

def waverec2(coefs, wavelet='db2', shape=None):
    """
    Performs a 2D wavelet reconstruction (inverse transform).

    In :func:`wavedec2`, you specify `levels`, which is not done in this function since it can be inferred from the shape of `coefs`.
    
    .. note::

        This function runs :func:`daubechies_factory` for you and caches the value. This means that first time you call it,
        performance will be slower than expected. You will also incur a dictionary lookup, which might cost you 100 ns. 

    .. seealso::

        :ref:`wavelet2d`

    Parameters
    ----------
    A : ndarray
        2D input data. Shape must be powers of two and square.
    wavelet : str
        Wavelet type. See :func:`daubechies_factory`.
    """
    global _db_wavedec2_cache, _db_waverec2_cache
    tup = (shape or coefs.shape, wavelet)
    try:
        rec = _db_waverec2_cache[tup]
    except KeyError:
        dec, rec = daubechies_factory(*tup) 
        _db_wavedec2_cache[tup] = dec
        _db_waverec2_cache[tup] = rec 

    return rec(coefs) 

# HELPER FUNCTIONS
################################################################################

def smart_flatten(coefficients):
    """
    This flattens 2D coefficients in a smart way, so that all coefficients levels are grouped into contiguous blocks, starting from the low-frequency coefficients going to the high-frequency ones.

    Notice that 1D coefficients are already flat and sorted.

    Parameters
    ----------
    coefficients : ndarray
        Wavelet coefficients returned by :func:`wavedec2`.
    """
    assert coefficients.shape == (8, 8), "TODO: Has not been generalized, only works with shape (8, 8), not {0}".format(coefficients.shape)
    olds = np.zeros(64)
    olds[0] = coefficients[0,0]
    olds[1] = coefficients[1,0]
    olds[2] = coefficients[0,1]
    olds[3] = coefficients[1,1]
    olds[4:8] = coefficients[2:4,0:2].flatten()
    olds[8:12] = coefficients[0:2,2:4].flatten()
    olds[12:16] = coefficients[2:4,2:4].flatten()
    olds[16:32] = coefficients[4:8,0:4].flatten()
    olds[32:48] = coefficients[0:4,4:8].flatten()
    olds[48:64] = coefficients[4:8,4:8].flatten()
    return olds 

def smart_deflatten(flatcoefs):
    """
    Inverse function of :func:`smart_flatten`.

    Parameters
    ----------
    flatcoefs : ndarray
        Flat array of coefficients returned by :func:`smart_flatten`.
    """
    N = int(np.sqrt(len(flatcoefs)))
    A = np.arange(N*N, dtype=int).reshape(N, N)
    indices = new2old(A).astype(int)
    new_indices = np.empty(indices.shape, dtype=int)
    for i, index in enumerate(indices):
        new_indices[index] = i
    news = flatcoefs[new_indices].reshape(8, 8).copy()
    return news 

def structured_to_contiguous(structured_coefs):
    """
    Converts a structured list-of-tuples-of-arrays-of-coefficients to a contiguous block.

    The input format follows `PyWavelets <http://www.pybytes.com/pywavelets/>`_.

    Works for both 1D and 2D coefficients.

    Parameters
    ----------
    structured_coefs : list
        List of coefficients. 
    """
    in2d = structured_coefs[0][0].ndim == 2 
    
    if in2d:
        N = 1 << (len(structured_coefs)-1)
        u = np.zeros((N, N))
        u[0,0] = float(structured_coefs[0])
        for level, c in enumerate(structured_coefs): 
            if level != 0:
                S = len(c[0]) 
                u[S:2*S,:S] = c[0]
                u[:S,S:2*S] = c[1]
                u[S:2*S,S:2*S] = c[2]
        return u
    else:
        N = 1 << (len(structured_coefs)-1)
        u = np.zeros(N)
        u[0] = float(structured_coefs[0])
        for level, c in enumerate(structured_coefs): 
            if level != 0:
                S = len(c) 
                u[S:2*S] = c
        return u

def contiguous_to_structured(contiguous_coefs, levels=np.inf):
    """
    Convert from continguous array to a structured format (identical to the one used in PyWavelets).

    Works for both 1D and 2D coefficients.

    Parameters
    ---------- 
    contiguous_coefs : ndarray
        Coefficients as returned by our wavelet functions.
    levels : int, optional
        If you don't want all levels, you can set this value to specify how many you want. Notice that this
        refers to levels of `wavelet` coefficients, which means that the scaling coefficient is not included and
        will always be returned, even if `levels` is set to zero. 
    """
    u = contiguous_coefs
    in2d = contiguous_coefs.ndim == 2 
    N = int(np.log2(len(contiguous_coefs)))
    
    coefs = []
    if in2d:
        coefs.append( contiguous_coefs[:1,:1] )
        for level in range(min(levels, N)):
            S = 1 << level
            coefs.append( (contiguous_coefs[S:2*S,:S], contiguous_coefs[:S,S:2*S], contiguous_coefs[S:2*S,S:2*S]) )
    else:
        coefs.append( contiguous_coefs[:1] )
        for level in range(min(levels, N)):
            S = 1 << level
            coefs.append( contiguous_coefs[S:2*S] ) 
        
    return coefs 