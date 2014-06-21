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
    
def all_sig_mask(sig_mask, xLen=491, yLen=673):
    """
    Create a mask of the biological data from all the signatures.
    """
    aggre_sig_mask = np.zeros(np.squeeze(sig_mask[0]).shape[0])
    for s_idx in np.arange(len(sig_mask)):
        aggre_sig_mask[np.where(np.squeeze(sig_mask[s_idx]))] = 1
    
    aggre_sig_mask = np.reshape(aggre_sig_mask, (yLen, xLen))
    
    return aggre_sig_mask
def sig_list_to_mask(sig_list, xLen=491, yLen=673):
    """
    Converts a signature list to a mask.
    """
    sig_masks = []
    for sig in np.arange(len(sig_list)):
        sigs_arr = np.zeros((1, xLen*yLen))
        if np.squeeze(sig_list[sig]).shape:
            sigs_arr[0, (np.squeeze(sig_list[sig]),)] = 1
            sig_masks.append(sigs_arr)
            
    return sig_masks

def pearsons_correlation_coeff(vector1, vector2):
    """
    Pearson's correlation coefficient between two vectors.
    """
    return np.dot(vector1, vector2)/(np.sqrt(np.sum(vector1**2))*np.sqrt(np.sum(vector2**2)))