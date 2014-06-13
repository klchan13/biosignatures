"""
Useful tools for biosignature analyses.
"""

import numpy as np

def signature_map(sig_masks, xLen, yLen):
    """
    Create a color-coded signature map.
    """
    sig_map = np.zeros((yLen, xLen))
    for sig_idx in np.arange(len(sig_masks)):
        sig_mask = np.reshape(sig_masks[sig_idx], (yLen, xLen))
        sig_map[np.where(sig_mask)] = sig_idx + 1
                
    return sig_map

def pearsons_correlation_coeff(vector1, vector2):
    """
    Pearson's correlation coefficient between two vectors.
    """
    return np.dot(vector1, vector2)/(np.sqrt(np.sum(vector1**2))*np.sqrt(np.sum(vector2**2)))