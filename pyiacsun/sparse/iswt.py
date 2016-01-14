import pywt
import numpy as np

def iswt(coefficients, wavelet): 
    """ 
    Input parameters: 

        coefficients 
          approx and detail coefficients, arranged in level value 
          exactly as output from swt: 
          e.g. [(cA1, cD1), (cA2, cD2), ..., (cAn, cDn)] 

        wavelet 
          Either the name of a wavelet or a Wavelet object 

    """ 
    output = coefficients[0][0].copy() # Avoid modification of input data 

  #num_levels, equivalent to the decomposition level, n 
    num_levels = len(coefficients) 
    for j in range(num_levels,0,-1): 
        step_size = int(pow(2, j-1)) 
        last_index = step_size 
        _, cD = coefficients[num_levels - j] 
        for first in range(last_index): # 0 to last_index - 1 
      # Getting the indices that we will transform 
          indices = np.arange(first, len(cD), step_size) 
          #print first, indices 

          # select the even indices 
          even_indices = indices[0::2] 
          # select the odd indices 
          odd_indices = indices[1::2] 

          # perform the inverse dwt on the selected indices, 
          # making sure to use periodic boundary conditions 
          x1 = pywt.idwt(output[even_indices], cD[even_indices], wavelet, 'per') 
          x2 = pywt.idwt(output[odd_indices], cD[odd_indices], wavelet, 'per') 
          # perform a circular shift right 

          # original: 
          #x2 = roll(x2, 1) 
          # average and insert into the correct indices 
          #output[indices] = (x1 + x2)/2. 

          #modified to allow exact reconstruction of original data, if swt2 is used 
          # with start_level = 0, and wavelet is haar or db1 
          output[even_indices] = x1[0::2] 
          output[odd_indices] = x2[0::2] 

    return output 