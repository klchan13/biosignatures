"""
Useful tools for biosignature analyses.
"""

import numpy as np
from math import factorial as f

def signature_map(sigs, xLen=491, yLen=673):
    """
    Create a color-coded signature map.
    """
    sig_map = np.zeros((yLen, xLen))
    
    # For keeping track of the number of actual signatures
    # (the ones with more than 1 pixel in them)
    big_sigs = np.zeros(len(sigs))
    non_sigs = False
    for i in np.arange(len(sigs)):
        if np.shape(np.squeeze(sigs[i])):
            big_sigs[i] = 1
        else:
            non_sigs = True
            
    count = 0
    inds_ref = np.where(np.ones((yLen, xLen)))
    for sig_idx in np.arange(len(sigs)):
        # The only numbers within the mask are 0s and 1s.  If not,
        # then it's a sig list and not a mask.
        uniq_sigs = np.unique(sigs[0])
        if (uniq_sigs[0] != 0) & (uniq_sigs[1] != 1):
            sig_inds = (inds_ref[0][sigs[sig_idx]], inds_ref[1][sigs[sig_idx]],)
        else:
            sig_mask = np.reshape(sigs[sig_idx], (yLen, xLen))
            sig_inds = np.where(sig_mask)
        
        # If there is only one pixel, give them the idx
        # of maximum number of real signatures.
        if np.shape(np.squeeze(sigs[sig_idx])):
            sig_map[sig_inds] = sig_idx+1#count + 1
            #count = count + 1
        else:
            sig_map[sig_inds] = len(np.where(big_sigs)[0])
    
    return sig_map, non_sigs
    
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
    
def nchoosek(n,k):
    """
    Finds all the number of unique combinations from choosing groups of k from a pool of n.
    
    Parameters
    ----------
    n: int
        Number of items in the pool you are choosing from
    k: int
        Size of the groups you are choosing from the pool
        
    n!/(k!*(n-k)!)
    """
    return f(n)/f(k)/f(n-k)

def pearsons_correlation_coeff(vector1, vector2):
    """
    Pearson's correlation coefficient between two vectors.
    """
    return np.dot(vector1, vector2)/(np.sqrt(np.sum(vector1**2))*np.sqrt(np.sum(vector2**2)))