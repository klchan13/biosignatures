"""
Functions for performing various analyses on the output signatures found through
the MATLAB program.
"""

import numpy as np
import scipy.stats.stats as stats
import biosignatures.utils as bsu
import sys
import time

def signature(lin_data, s_idx, sig_masks=None, sig_list=None):
    """
    Find the desired signature's elemental values.
    """
    if sig_masks != None:
        return lin_data[:, np.where(np.squeeze(sig_masks[s_idx]))]
    else:
        return lin_data[:, sig_list[s_idx]]
    
def signature_medians(lin_data, sig_masks=None, sig_list=None):
    """
    Find the signature medians.
    """
    # What is the format the signatures are in?
    if sig_masks != None:
        in_sigs = sig_masks
    else:
        in_sigs = sig_list
        
    sig_meds = np.zeros((len(in_sigs), np.shape(lin_data)[0]))
    for sig_idx in np.arange(len(in_sigs)):
        if sig_masks != None:
            sig_len = len(np.where(np.squeeze(sig_masks[sig_idx]))[0])
            this_sig = np.squeeze(signature(lin_data, sig_idx, sig_masks=sig_masks))
        else:
            sig_len = len(sig_list[sig_idx])
            this_sig = np.squeeze(signature(lin_data, sig_idx, sig_list=sig_list))
        
        # Take the signatures as is instead of taking the median if there is only
        # one pixel.
        if sig_len > 1:
            sig_meds[sig_idx] = np.median(this_sig,-1)
        else:
            sig_meds[sig_idx] = this_sig
                
    return sig_meds
    
def sig_iqr(lin_data, sig_masks=None, sig_list=None):
    """
    Evaluating the interquartile ranges for each element of each signature.
    """
    if sig_masks != None:
        q25s = sig_quartile(lin_data, 25, sig_masks=sig_masks)
        q75s = sig_quartile(lin_data, 75, sig_masks=sig_masks)
    else:
        q25s = sig_quartile(lin_data, 25, sig_list=sig_list)
        q75s = sig_quartile(lin_data, 75, sig_list=sig_list)    
    
    return q75s - q25s
    
def sig_quartile(lin_data, percentile, sig_masks=None, sig_list=None):
    """
    Evaluating the quartile for each element of each signature.
    """
    # What is the format the signatures are in?
    if sig_masks != None:
        in_sigs = sig_masks
    else:
        in_sigs = sig_list
        
    quartile = np.zeros((len(in_sigs), lin_data.shape[0]))
    for s_idx in np.arange(len(in_sigs)):
        if sig_masks != None:
            this_sig = np.squeeze(signature(lin_data, s_idx, sig_masks=sig_masks))
        else:
            this_sig = np.squeeze(signature(lin_data, s_idx, sig_list=sig_list))
        
        quartile[s_idx] = stats.scoreatpercentile(this_sig, percentile, axis=-1)
        
    return quartile
def reshape_sig_list(sig_list):
    """
    Reshapes the sig_list so it has the right size that can be indexed.
    """
    for idx in np.arange(len(sig_list)):
        if np.squeeze(sig_list[idx]).shape:
            sig_list[idx] = np.squeeze(sig_list[idx])
        else:
            sig_list[idx] = np.reshape(sig_list[idx], (1,))
    
    return sig_list
    
def reclass(lin_data_alt, sigList=None, allSig=None, minReassignPix=115):
    """
    For reclassifying over and over and over...until the end of time!
    (Or just until you have minimized the drift problem as much as you can.)
    
    This does the main computation and is re-adapted from the MATLAB code.
    """
    if sigList == None:
        sigList = []
        for idx in np.arange(len(allSig)):
            sigList.append(np.where(np.squeeze(allSig[idx]))[0])
    else:
        # Reshaping the signatures with only one pixel in them.
        sigList = reshape_sig_list(sigList)
        
    new_sig_list = []
    if allSig != None:
        medians = signature_medians(lin_data_alt, sig_masks=allSig)
    else:
        medians = signature_medians(lin_data_alt, sig_list=sigList)
        
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
                if (minReassignPix is not None) & (len(sigList[in_idx]) == 1):
                    p_ref_arr[in_idx] = 0
                elif (minReassignPix is None) | (len(sigList[in_idx]) > minReassignPix) | (in_idx == ref_idx):
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
                                                               np.array([sigList[ai][
                                                               where_misfits[0][m_idx]]])))
            else:
                reassign[int(misfits[m_idx])] = np.array([sigList[ai][where_misfits[0][m_idx]]])           

    # Now add the misfits to their appropriate signature
    fin_sig_list = []
    for ai in np.arange(len(sigList)):
        real_vox = np.where(new_sig_list[ai]+1 == (ai+1)*np.ones(len(new_sig_list[ai])))
        if reassign[ai].shape:
            fin_sig_list.append(np.squeeze(np.concatenate((sigList[ai][real_vox],
                                                                  reassign[ai]))))
        else:
            fin_sig_list.append(sigList[ai][real_vox])
    
    end_sig_list = []
    for sig in np.arange(len(fin_sig_list)):
        if len(fin_sig_list[sig]) > 0:
            end_sig_list.append(fin_sig_list[sig])
            
    return end_sig_list
    
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
    sig_meds_forw = signature_medians(lin_data_alt, sig_masks=sig_forwards)
    sig_meds_back = signature_medians(lin_data_alt, sig_masks=sig_backwards)
    
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
            idx = abs(less_sz - more_sigs[1])/float(less_sz) < sz_diff
            this_more_sig = more_sigs[0][idx]
        else:
            this_more_sig = more_sigs[0]
            
        this_cc_arr = np.zeros(len(this_more_sig))
        
        # Find the signature in more_sig closest to the one in less_sigs
        for more_idx in np.arange(len(this_more_sig)):
            this_cc_arr[more_idx] = bsu.pearsons_correlation_coeff(less_sigs[0][less_idx], this_more_sig[more_idx])#[0] stats.pearsonr
        
        cc_arr[less_idx] = max(this_cc_arr)
        if sz_diff != None:
            sim_sig[less_idx] = np.where(idx)[0][np.where(this_cc_arr == cc_arr[less_idx])]
        else:
            sim_sig[less_idx] = np.where(this_cc_arr == max(this_cc_arr))[0]
        
    return cc_arr, sim_sig, less_sigs_idx
    
def twin_sigs(lin_data_alt, sig1, sig2, sz_diff=None):
    """
    Finding twin signatures across different seeding locations for determining
    signature error.
    """
    import scipy.stats.stats as stats
    
    count = 0
    # Initialize with a rough estimation of what the twin signatures are and
    # their corresponding correlation.  This can result in two different
    # signatures having the same twin signatures though which is not
    # physically possible.
    cc_arr, sim_sig, less_sigs_idx = sig_reliability(lin_data_alt,
                                                          sig1, sig2,
                                                          sz_diff=sz_diff)
    sig_list = [sig1, sig2]
    twin_sigs = []
    cc_twins = []
    
    # Prepare a mask the length of both original signatures for
    # tracking purposes.
    more_inds_mask = np.ones(len(np.arange(len(sig_list[
                               1-less_sigs_idx]))))
    less_inds_mask = np.ones(len(np.copy(np.arange(len(sig_list[
                                        less_sigs_idx])))))
    
    # Run this loop as long as there are no more repeats of signatures
    # of the similar signatures array.
    while stats.mode(sim_sig)[1] > 1:
        if count == 0:
            count = count + 1
        elif count > 0:
            # Find the reliability and twin sigs of the reduced arrays
            cc_arr, sim_sig, _ = sig_reliability(lin_data_alt,
                                                             new_sig_list1,
                                                             new_sig_list2,
                                                             sz_diff=sz_diff)

        new_sig_list1 = []
        new_sig_list2 = []
        
        # Find the signature that's repeated and reduce the cc_arr to only
        # contain the rows where the repeat occurs.
        red_repeat_sig = int(stats.mode(sim_sig)[0])
        repeated_sig = np.where(more_inds_mask)[0][red_repeat_sig]
        
        more_inds_mask[repeated_sig] = 0 # Update more sigs tracking array
        red_cc_arr = cc_arr[np.where(sim_sig == red_repeat_sig)]
        
        # Find the signature of the lesser sigs that has the highest CC
        # with the repeated sig
        twin_to_rep_bool = cc_arr == np.max(red_cc_arr)
        cc_twins.append(cc_arr[np.where(twin_to_rep_bool)])
        
        where_twin_mask = np.where(less_inds_mask)[0][twin_to_rep_bool]
        twin_sigs.append(np.array([where_twin_mask[0],
                                   repeated_sig])[..., None])
        less_inds_mask[where_twin_mask] = 0

        # Find the reliability and twin sigs of the reduced arrays
        for sig_idx in np.where(less_inds_mask)[0]:
            new_sig_list1.append(sig_list[less_sigs_idx][sig_idx])
        
        for sig_idx in np.where(more_inds_mask)[0]:
            new_sig_list2.append(sig_list[1-less_sigs_idx][sig_idx])
    
    # Add the rest of the signatures that don't have duplicates
    # within sim_sig to the twin sig array.
    if count == 0: # This means there were no twin sigs to begin with
        mor_sigs = sim_sig.astype(int)[None]
    else:
        cc_arr, sim_sig, _ = sig_reliability(lin_data_alt, new_sig_list1,
                                         new_sig_list2, sz_diff=sz_diff)
        mor_sigs = np.squeeze(np.where(more_inds_mask)[0]
                             [sim_sig.astype(int)])[None]

    twin_sigs.append(np.concatenate((np.where(less_inds_mask)[0][None],
                                     mor_sigs)))
                                     
    cc_twins.append(cc_arr)
    cc_twins_arr = np.concatenate(cc_twins)
    twin_sigs_arr = np.concatenate(twin_sigs, -1)
    
    return cc_twins_arr, twin_sigs_arr, less_sigs_idx
    
def sig_diffs(lin_data_alt, sig1, sig2, sig3=None, sz_diff=None):
    """
    Error metric of signatures groups found using different signature
    classification algorithms.
    """
    sig_map1, non_sigs1 = bsu.signature_map(sig1)
    sig_map2, non_sigs2 = bsu.signature_map(sig2)
    
    uniq_sigs1 = np.unique(sig1[0])
    uniq_sigs2 = np.unique(sig2[0])
    
    # Reduce the signature list to the real signatures and change from a list
    # to a mask
    if (uniq_sigs1[0] != 0) & (uniq_sigs1[1] != 1):
        sig1 = bsu.sig_list_to_mask(bsu.reduce_sig_list(sig1))
                
    if (uniq_sigs2[0] != 0) & (uniq_sigs2[1] != 1):
        sig2 = bsu.sig_list_to_mask(bsu.reduce_sig_list(sig2))
        
    cc_twins, twin_sigs_fin, less_sigs_idx = twin_sigs(lin_data_alt, sig1, sig2)
    sig_map_list = [sig_map1, sig_map2]
    
    if sig3 == None:
        all_sigs_mask = bsu.all_sig_mask(sig1)
    else:
        all_sigs_mask = bsu.all_sig_mask(sig3)
  
    coords = np.where(all_sigs_mask)
    err = 0
    diff_map = np.zeros(sig_map1.shape)
    
    # For all the pixels considered biological, check every pixel to see if the 
    # signature in one image corresponds to what is considered it's twin signature
    # in the second image.  If not, it's pixel is added to the running error total.
    for c_idx in np.arange(np.sum(all_sigs_mask)):
        this_coord = (coords[0][c_idx], coords[1][c_idx])
        this_sig1 = sig_map_list[less_sigs_idx][this_coord] - 1
        if (sig_map_list[1 - less_sigs_idx][this_coord] - 1 !=
            twin_sigs_fin[1][np.where(twin_sigs_fin[0] == this_sig1)]):
            err = err + 1
            # For display purposes: 2 is error
            diff_map[this_coord] = 2
        else:
            # ...and 1 is considered no error (the signature in the 
            # pixel in image 1 does not have its corresponding twin
            # signature in the same pixel in image 2)
            diff_map[this_coord] = 1
            
    # Extra error is where both images have unclassified pixels
    if (non_sigs1 == True) & (non_sigs2 == True):
        where_extra_err = np.where((sig_map1 == np.max(sig_map1)) &
                                   (sig_map2 == np.max(sig_map2)))
        extra_err = len(where_extra_err[0])
        diff_map[where_extra_err] = 2
    else:
        extra_err = 0 
        
    return err + extra_err, diff_map

def reclass_multi_data(data_sets, xLen=491, yLen=673, max_itr=25, minReassignPix=10, save=True):
    """
    For reclassifying multiple data sets at once and looking at the error between them
    at different iterations.
    
    This takes a while, so keep checking the progress bar.  Maybe grab a coffee or two.
    """
    t1 = time.time()
    all_masks_list = []
    for d_idx, data in enumerate(data_sets):
        # Initial reclassification after data generation.
        # Note: rand_reclass is in sig list form to save memory.
        
        # Change from minutes to hours if current run time is longer than an hour
        t2 = time.time()
        if ((t2-t1)/60.) > 60.:
            n = 60.
            units = "hours"
        else:
            n = 1.
            units = "mins"
            
        sys.stdout.write('\r' + "Reclassification %s for data set %s.  Elapsed time: %s %s"%(1, d_idx+1,(t2-t1)/(60.*n), units))
        sys.stdout.flush() 
        
        all_masks = []
        rand_reclass = reclass(lin_data_alt, sigList=data,
                               minReassignPix=minReassignPix)
        all_masks.append(np.squeeze(np.array(rand_reclass)))
        # Initialize an array to include all the masks from the reassignment.
        # After the first reclassification, all the 1 pixel signatures should
        # be gone and the number of signatures should be constant.

        for itr in np.arange(1, max_itr):
            t2 = time.time()
            sys.stdout.write('\r' + "Reclassification %s for data set %s.  Elapsed time: %s mins"%(itr+2, d_idx+1,(t2-t1)/60.))
            sys.stdout.flush() 
            
            # Reclassify and change from signature list to mask.
            rand_reclass = reclass(lin_data_alt, sigList=rand_reclass,
                                        minReassignPix=minReassignPix)
            all_masks.append(np.squeeze(np.array(rand_reclass)))
        
        # Save the masks for this data set.
        all_masks_list.append(all_masks)
        if save == True:
            np.save("sig_list_d%s.npy"%d_idx, all_masks)        

    aggre_sig_diffs = []
    med_sig_diffs = np.zeros(max_itr)
    for itr in np.arange(max_itr):
        these_sig_diffs = np.zeros(bsu.nchoosek(len(data_sets),2))
        for d_inds, di in enumerate(itertools.combinations(np.arange(len(data_sets)), 2)):
            # Find the cluster differences between mask lists from different data sets
            # at the same iteration of reassignment.
            t2 = time.time()
            
            # Change from minutes to hours if current run time is longer than an hour
            if ((t2-t1)/60.) > 60.:
                n = 60.
                units = "hours"
            else:
                n = 1.
                units = "mins"
                
            sys.stdout.write('\r' + "Finding differences of reclassification %s for data pairs %s of %s.  Elapsed time: %s %s"%(itr,
                                                                        d_inds+1, bsu.nchoosek(len(data_sets),2),(t2-t1)/60., units))
            sys.stdout.flush()
            
            sig_diff, diff_map = sig_diffs(lin_data_alt,
                                           bsu.sig_list_to_mask(all_masks_list[di[0]][itr]),
                                           bsu.sig_list_to_mask(all_masks_list[di[1]][itr]))
            these_sig_diffs[d_inds] = sig_diff
        aggre_sig_diffs.append(these_sig_diffs)
        med_sig_diffs[itr] = np.median(these_sig_diffs)
        
    # Save the signature errors differences and their medians.
    if save == True:
        np.save("aggre_sig_diffs.npy", aggre_sig_diffs)
        np.save("med_sig_diffs.npy", med_sig_diffs)
    
    return all_masks_list, aggre_sig_diffs, med_sig_diffs