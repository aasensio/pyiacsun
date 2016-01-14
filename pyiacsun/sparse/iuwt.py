import numpy as np
import ctypes

# Usage
# a = np.random.randn(64,64)
# nLevels = 3
# scaleAdjust = 0
# res = iuwt.ser_iuwt_decomposition(a, nLevels, scaleAdjust, True)
# out = iuwt.ser_iuwt_recomposition(res[0], scaleAdjust, res[1])
# out-a = 0

def iuwt_decomposition(in1, scale_count, scale_adjust, store_smoothed):
    """
    This function calls the a trous algorithm code to decompose the input into its wavelet coefficients. This is
    the isotropic undecimated wavelet transform implemented for a single CPU core.

    INPUTS:
    in1                 (no default):   Array on which the decomposition is to be performed.
    scale_count         (no default):   Maximum scale to be considered.
    scale_adjust        (default=0):    Adjustment to scale value if first scales are of no interest.
    store_smoothed      (default=False):Boolean specifier for whether the smoothed image is stored or not.

    OUTPUTS:
    detail_coeffs                       Array containing the detail coefficients.
    C0                  (optional):     Array containing the smoothest version of the input.
    """

    wavelet_filter = (1./16)*np.array([1,4,6,4,1])      # Filter-bank for use in the a trous algorithm.

    # Initialises an empty array to store the coefficients.

    detail_coeffs = np.empty([scale_count-scale_adjust, in1.shape[0], in1.shape[1]])

    C0 = in1    # Sets the initial value to be the input array.

    # The following loop, which iterates up to scale_adjust, applies the a trous algorithm to the scales which are
    # considered insignificant. This is important as each set of wavelet coefficients depends on the last smoothed
    # version of the input.

    if scale_adjust>0:
        for i in range(0, scale_adjust):
            C0 = ser_a_trous(C0, wavelet_filter, i)

    # The meat of the algorithm - two sequential applications fo the a trous followed by determination and storing of
    # the detail coefficients. C0 is reassigned the value of C on each loop - C0 is always the smoothest version of the
    # input image.

    for i in range(scale_adjust,scale_count):
        C = ser_a_trous(C0, wavelet_filter, i)                                  # Approximation coefficients.
        C1 = ser_a_trous(C, wavelet_filter, i)                                  # Approximation coefficients.
        detail_coeffs[i-scale_adjust,:,:] = C0 - C1                             # Detail coefficients.
        C0 = C

    if store_smoothed:
        return detail_coeffs, C0
    else:
        return detail_coeffs

def iuwt_recomposition(in1, scale_adjust, smoothed_array):
    """
    This function calls the a trous algorithm code to recompose the input into a single array. This is the
    implementation of the isotropic undecimated wavelet transform recomposition for a single CPU core.

    INPUTS:
    in1             (no default):   Array containing wavelet coefficients.
    scale_adjust    (no default):   Indicates the number of truncated array pages.
    smoothed_array  (default=None): For a complete inverse transform, this must be the smoothest approximation.

    OUTPUTS:
    recomposition                   Array containing the reconstructed image.
    """

    wavelet_filter = (1./16)*np.array([1,4,6,4,1])      # Filter-bank for use in the a trous algorithm.

    # Determines scale with adjustment and creates a zero array to store the output, unless smoothed_array is given.

    max_scale = in1.shape[0] + scale_adjust

    if smoothed_array is None:
        recomposition = np.zeros([in1.shape[1], in1.shape[2]])
    else:
        recomposition = smoothed_array

    # The following loops call the a trous algorithm code to recompose the input. The first loop assumes that there are
    # non-zero wavelet coefficients at scales above scale_adjust, while the second loop completes the recomposition
    # on the scales less than scale_adjust.

    for i in range(max_scale-1, scale_adjust-1, -1):
        recomposition = ser_a_trous(recomposition, wavelet_filter, i) + in1[i-scale_adjust,:,:]

    if scale_adjust>0:
        for i in range(scale_adjust-1, -1, -1):
            recomposition = ser_a_trous(recomposition, wavelet_filter, i)

    return recomposition

def ser_a_trous(C0, filter, scale):
    """
    The following is a serial implementation of the a trous algorithm. Accepts the following parameters:

    INPUTS:
    filter      (no default):   The filter-bank which is applied to the components of the transform.
    C0          (no default):   The current array on which filtering is to be performed.
    scale       (no default):   The scale for which the decomposition is being carried out.

    OUTPUTS:
    C1                          The result of applying the a trous algorithm to the input.
    """
    tmp = filter[2]*C0

    tmp[(2**(scale+1)):,:] += filter[0]*C0[:-(2**(scale+1)),:]
    tmp[:(2**(scale+1)),:] += filter[0]*C0[(2**(scale+1))-1::-1,:]

    tmp[(2**scale):,:] += filter[1]*C0[:-(2**scale),:]
    tmp[:(2**scale),:] += filter[1]*C0[(2**scale)-1::-1,:]

    tmp[:-(2**scale),:] += filter[3]*C0[(2**scale):,:]
    tmp[-(2**scale):,:] += filter[3]*C0[:-(2**scale)-1:-1,:]

    tmp[:-(2**(scale+1)),:] += filter[4]*C0[(2**(scale+1)):,:]
    tmp[-(2**(scale+1)):,:] += filter[4]*C0[:-(2**(scale+1))-1:-1,:]

    C1 = filter[2]*tmp

    C1[:,(2**(scale+1)):] += filter[0]*tmp[:,:-(2**(scale+1))]
    C1[:,:(2**(scale+1))] += filter[0]*tmp[:,(2**(scale+1))-1::-1]

    C1[:,(2**scale):] += filter[1]*tmp[:,:-(2**scale)]
    C1[:,:(2**scale)] += filter[1]*tmp[:,(2**scale)-1::-1]

    C1[:,:-(2**scale)] += filter[3]*tmp[:,(2**scale):]
    C1[:,-(2**scale):] += filter[3]*tmp[:,:-(2**scale)-1:-1]

    C1[:,:-(2**(scale+1))] += filter[4]*tmp[:,(2**(scale+1)):]
    C1[:,-(2**(scale+1)):] += filter[4]*tmp[:,:-(2**(scale+1))-1:-1]

    return C1
 
 
def iuwt_threshold(in1, sigma_level=4):
    """
    This function performs the thresholding of the values in array in1 based on the estimated standard deviation
    given by the MAD (median absolute deviation) estimator about zero.
    INPUTS:
    in1             (no default):   The array which is to be thresholded.
    sigma_level     (no default):   The number of estimated deviations at which thresholding is to occur.
    OUTPUTS:
    out1                            An thresholded version of in1.
    """

    out1 = np.empty_like(in1)

    # The conditional here ensures that the function works even when only one scale is considered. Both cases are the
    # same: the MAD estimator is calculated and then the resulting value is used to threshold the input. NOTE: This
    # discards all negative coefficients.

    if len(in1.shape)==2:
        threshold_level = np.median(np.abs(in1))/0.6745                     # MAD estimator for normal distribution.
        out1 = (in1>(sigma_level*threshold_level))*in1
    else:
        for i in range(in1.shape[0]):
            threshold_level = np.median(np.abs(in1[i,:,:]))/0.6745          # MAD estimator for normal distribution.
            out1[i,:,:] = (in1[i,:,:]>(sigma_level*threshold_level))*in1[i,:,:]

    return out1


# a = np.random.randn(64,1)
# nLevels = 3
# scaleAdjust = 0
# res = ser_iuwt_decomposition(a, nLevels, scaleAdjust, True)
# out = ser_iuwt_recomposition(res[0], scaleAdjust, res[1])
# print(np.allclose(a,out))

# resThresholded = threshold(res[0])
# outThresholded = ser_iuwt_recomposition(resThresholded, scaleAdjust, res[1])