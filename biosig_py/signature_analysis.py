"""
Functions for performing various analyses on the output signatures found through
the MATLAB program.
"""

import numpy as np

def signature(sig_masks, lin_data, s_idx):
    """
    Find the desired signature's elemental values.
    """
    return lin_data[:, np.where(np.squeeze(sig_masks[s_idx]))]
    
def signature_medians(sig_masks, lin_data):
    """
    Find the signature medians.
    """
    sig_meds = np.zeros((len(sig_masks), np.shape(lin_data_alt)[0]))
    for sig_idx in np.arange(len(sig_forwards)):
        sig_meds[sig_idx] = np.squeeze(np.median(np.squeeze(signature(sig_masks,
                                                         lin_data, s_idx))), -1)
                
    return sig_meds
    
def sig_iqr(lin_data, sigs):
    """
    Evaluating the interquartile ranges for each element of each signature.
    """
    q25s = sig_quartile(lin_data, sigs, 25)
    q75s = sig_quartile(lin_data, sigs, 75)
    
    return q75s - q25s
    
def sig_quartile(lin_data, sigs, percentile):
    """
    Evaluating the quartile for each element of each signature.
    """
    quartile = np.zeros((len(sigs), lin_data.shape[0]))
    for s_idx in np.arange(len(sigs)):
        quartile[s_idx] = stats.scoreatpercentile(np.squeeze(sa.signature(sigs, lin_data, s_idx)),
                                                  percentile, axis=-1)
        
    return quartile
    
def sig_reliability(lin_data_alt, sig_forwards, sig_backwards, sz_diff=0.4):
    """
    Measuring the reliability of the signature finding algorithm using Pearson's Correlation Coefficient.
    """
    sig_meds_forw = signature_medians(sig_forwards, lin_data_alt)
    sig_meds_forw = signature_medians(sig_forwards, lin_data_alt)

    sig_forw_sz = np.zeros(len(sig_forwards))
    sig_back_sz = np.zeros(len(sig_backwards))
    
    # Find all the sizes for each signature
    for sig_idx in np.arange(len(sig_forwards)):
        sig_forw_sz[sig_idx] = np.sum(sig_forwards[sig_idx])
        
    for sig_idx in np.arange(len(sig_backwards)):
        sig_back_sz[sig_idx] = np.sum(sig_backwards[sig_idx])
    
    # If unequal number of signatures in forwards and backwards,
    # find CC for the one with less signatures
    if (len(sig_forwards) < len(sig_backwards)) | (len(sig_forwards) == len(sig_backwards)):
        less_sigs = [sig_meds_forw, sig_forw_sz]
        more_sigs = [sig_meds_back, sig_back_sz]
    else:
        less_sigs = [sig_meds_back, sig_back_sz]
        more_sigs = [sig_meds_forw, sig_forw_sz]
        
    cc_arr = np.zeros(len(less_sigs[0]))
    sim_sig = np.zeros(len(less_sigs[0]))
    # Now find the reliability
    for less_idx in np.arange(len(less_sigs[0])):
        # Find signatures within more_sigs with comparable size to the current sig in less_sigs
        less_sz = less_sigs[1][less_idx]
        idx = abs(less_sz*np.ones(len(more_sigs[1])) - more_sigs[1])/float(less_sz) < sz_diff
        this_more_sig = more_sigs[0][idx]
        this_cc_arr = np.zeros(len(this_more_sig))
        
        # Find the signature in more_sig closest to the one in less_sigs
        for more_idx in np.arange(len(this_more_sig)):
            this_cc_arr[more_idx] = stats.pearsonr(less_sigs[0][less_idx], this_more_sig[more_idx])[0]
            
        cc_arr[less_idx] = max(this_cc_arr)
        sim_sig[less_idx] = np.where(idx)[0][np.where(this_cc_arr == cc_arr[less_idx])]
        
    return cc_arr, sim_sig