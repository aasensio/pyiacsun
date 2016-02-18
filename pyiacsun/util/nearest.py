__all__ = ['nearest']

import numpy as np


def nearest(array,value):
    """
    Find the index of the array that is close to value
    
    Args:
        array (array): array to be tested
        value (float): value to be tested
    
    Returns:
        int: index    
    """
    return (np.abs(array-value)).argmin()