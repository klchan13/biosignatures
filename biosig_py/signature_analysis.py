"""
Functions for performing various analyses on the output signatures found through
the MATLAB program.
"""

import numpy as np
import scipy.stats.stats as stats
import biosignatures.utils as bsu

def signature(sig_masks, lin_data, s_idx):
    """
    Find the desired signature's elemental values.
    """
    return lin_data[:, np.where(np.squeeze(sig_masks[s_idx]))]
    
def signature_medians(sig_masks, lin_data):
    """
    Find the signature medians.
    """
    sig_meds = np.zeros((len(sig_masks), np.shape(lin_data)[0]))
    for sig_idx in np.arange(len(sig_masks)):
        sig_meds[sig_idx] = np.median(np.squeeze(signature(sig_masks,
                                                 lin_data, sig_idx)),-1)
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
        quartile[s_idx] = stats.scoreatpercentile(np.squeeze(signature(sigs, lin_data, s_idx)),
                                                  percentile, axis=-1)
        
    return quartile
    
def reclass(lin_data_alt, allSig, minReassignPix=115):
    """
    For reclassifying over and over and over...until the end of time!
    (Or just until you have reduced the IQR as much as you can.)
    
    This does the main computation and is re-adapted from the MATLAB code.
    """
    sigList = []
    for idx in np.arange(len(allSig)):
        sigList.append(np.where(np.squeeze(allSig[idx]))[0])
        
    new_sig_list = []    
    medians = signature_medians(allSig, lin_data_alt)
    for s_idx in np.arange(len(sigList)):
        new_sig_list.append(np.zeros(len(sigList[s_idx])))
    
    # Go through each pixels signature in each rough aggregate and find the
    # correlation between that and the median of each signature.
    for ref_idx in np.arange(len(sigList)):
        for pix_sig_idx in np.arange(len(sigList[ref_idx])):
            pixel = sigList[ref_idx][pix_sig_idx]
            pix_sig = lin_data_alt[:, pixel]
            
            # Keep track of the correlation between the current pixel signature
            # and the median of the other rough aggregate.
            p_ref_arr = np.zeros(len(sigList))
            for in_idx in np.arange(len(sigList)):
                # You don't want to be accidentally assigning pixels from other
                # signatures to a small signature that may be biased (due 
                # simply to its small amount).  Since these signatures may not
                # be real, bias them towards assigning them to larger
                # signatures and don't let pixels from larger signatures be
                # assigned to the small signatures.  However, keep these
                # signatures as they may be real.
                if (minReassignPix is None) | ((len(sigList[in_idx]) > minReassignPix) | (in_idx == ref_idx)):
                    p_ref_arr[in_idx] = bsu.pearsons_correlation_coeff(medians[in_idx], pix_sig)
                else:
                    p_ref_arr[in_idx] = 0
            new_sig_list[ref_idx][pix_sig_idx] = np.where(p_ref_arr == max(p_ref_arr))[0][0]
    
    # Now use some fancy indexing to find the misfits and prep the pixels
    # for reassignment.
    reassign = list(np.arange(len(new_sig_list)))
    for ai in np.arange(len(new_sig_list)):
        where_misfits = np.where(new_sig_list[ai]+1 != (ai+1)*np.ones(len(new_sig_list[ai])))
        misfits = new_sig_list[ai][where_misfits]
        for m_idx in np.arange(len(where_misfits[0])):
            if reassign[int(misfits[m_idx])].shape:
                reassign[int(misfits[m_idx])] = np.concatenate((reassign[int(misfits[m_idx])],
                                                               np.array([sigList[ai][where_misfits[0][m_idx]]])))
            else:
                reassign[int(misfits[m_idx])] = np.array([sigList[ai][where_misfits[0][m_idx]]])           

    # Now add the misfits to their appropriate signature
    fin_sig_list = []
    for ai in np.arange(len(sigList)):
        real_vox = np.where(new_sig_list[ai]+1 == (ai+1)*np.ones(len(new_sig_list[ai])))
        if reassign[ai].shape:
            fin_sig_list.append(np.squeeze(np.concatenate((sigList[ai][real_vox], reassign[ai]))))
        else:
            fin_sig_list.append(sigList[ai][real_vox])
        
    return fin_sig_list
    
def refine_reclass(sig_no_reclass, lin_data, err_thresh=0.065, minReassignPix=None):
    """
    Keep reclassifying until the error between iterations no longer changes.
    """
    sigs = np.copy(sig_no_reclass)
    err_arr = np.array([5]) # Place holder before adding on the real errors.

    itr = 0
    # Continue reclassification if the sum of the differences in IQR is larger
    # than zero, which indicates that it hasn't converged at an answer yet.
    while err_arr[len(err_arr)-1] > err_thresh:
        itr = itr + 1
        
        # Run through reclassification algorithm
        fin_sig_list = reclass(lin_data, sigs, minReassignPix=minReassignPix)
        mask_list = []
        for idx in np.arange(len(fin_sig_list)):
            this_mask = np.zeros((1, np.shape(med_forw_sigs[idx])[1]))
            this_mask[0, (fin_sig_list[idx],)] = 1
            mask_list.append(this_mask)
        
        # Find the difference between the IQR of the signatures from the
        # previous reclassification and the IQR of the signatures from the
        # current reclassification.
        prev_iqr = sig_iqr(lin_data, sigs)
        cur_iqr = sig_iqr(lin_data, mask_list)
        this_err = np.sum(abs(prev_iqr - cur_iqr))
        
        if itr == 0:
            err_arr = np.array([this_err])
        else:
            err_arr = np.concatenate((err_arr, np.array([this_err])))
        
        sigs = mask_list
        
    return sigs, err_arr, itr
    
def sig_reliability(lin_data_alt, sig_forwards, sig_backwards, sz_diff=None):
    """
    Measuring the reliability of the signature finding algorithm using Pearson's Correlation Coefficient.
    """
    sig_meds_forw = signature_medians(sig_forwards, lin_data_alt)
    sig_meds_back = signature_medians(sig_backwards, lin_data_alt)

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
        less_sigs_idx = 0
    else:
        less_sigs = [sig_meds_back, sig_back_sz]
        more_sigs = [sig_meds_forw, sig_forw_sz]
        less_sigs_idx = 1
        
    cc_arr = np.zeros(len(less_sigs[0]))
    sim_sig = np.zeros(len(less_sigs[0]))
    # Now find the reliability
    for less_idx in np.arange(len(less_sigs[0])):
        # Find signatures within more_sigs with comparable size to the current sig in less_sigs
        less_sz = less_sigs[1][less_idx]
        if sz_diff != None:
            idx = abs(less_sz*np.ones(len(more_sigs[1])) - more_sigs[1])/float(less_sz) < sz_diff
            this_more_sig = more_sigs[0][idx]
        else:
            this_more_sig = more_sigs[0]
            
        this_cc_arr = np.zeros(len(this_more_sig))
        
        # Find the signature in more_sig closest to the one in less_sigs
        for more_idx in np.arange(len(this_more_sig)):
            this_cc_arr[more_idx] = stats.pearsonr(less_sigs[0][less_idx], this_more_sig[more_idx])[0]
            
        cc_arr[less_idx] = max(this_cc_arr)
        if sz_diff != None:
            sim_sig[less_idx] = np.where(idx)[0][np.where(idx)[0][np.where(this_cc_arr == cc_arr[less_idx])]]
        else:
            sim_sig[less_idx] = np.where(this_cc_arr == cc_arr[less_idx])[0]
        
    return cc_arr, sim_sig, less_sigs_idx