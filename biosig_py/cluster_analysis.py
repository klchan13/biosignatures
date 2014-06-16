
"""
Functions for performing various classification and analyses on the output clusters found through
the MATLAB program.
"""

import numpy as np
from scipy.ndimage.measurements import label

import biosignatures.signature_analysis as sa
import biosignatures.utils as bsu

def cluster_signatures(sig_out_in, bstrp_med, med_sig, p_thresh=0.05,
                       ele_order=np.array(["Fe", "Cu", "Zn", "Ca", "K", "S", "P", "Cl", "Si", "Mn"])):
    """
    Analyzes different signatures within a cluster group and finds the elemental trends.
    """
    for s_idx, sig in enumerate(sig_out_in):
        if s_idx != len(sig_out_in)-1:
            # Grab the two signatures of interest
            s1 = sig_out_in[s_idx]
            s2 = sig_out_in[s_idx+1]
            
            # Bootstrap data for the comparison of these two signatures
            this_bstrp = np.squeeze(bstrp_med[np.where((bstrp_med[:,0]==s1) & (bstrp_med[:,1]==s2))])
            
            # Find the elements within the two signatures with significant differences between them
            sig_diff_inds = np.squeeze(np.where(this_bstrp[2:] < p_thresh))
            red_s1 = med_sig[s1 - 1][sig_diff_inds]
            red_s2 = med_sig[s2 - 1][sig_diff_inds]
            
            bool_s1_larger = red_s1 > red_s2
            
            if bool_s1_larger.size > 1:
                if not all(~bool_s1_larger): # If not all of s2 > s1, then there are some elements where s1 > s2
                    s1_larger_ele = ele_order[sig_diff_inds[np.where(bool_s1_larger)]]
                    print "Significant decrease in %s inward from %s to %s"%(s1_larger_ele, s1, s2)
                
                if not all(bool_s1_larger): # If not all of s1 > s2, then there are some elements where s2 > s1
                    s2_larger_ele = ele_order[sig_diff_inds[np.where(~bool_s1_larger)]]
                    print "Significant increase in %s inward from %s to %s"%(s2_larger_ele, s1, s2)
            elif bool_s1_larger.size == 1:
                # If only one element in bool_s1_larger, then only one location of significant difference
                larger_ele = ele_order[sig_diff_inds]
                if bool_s1_larger:
                    print "Significant decrease in %s inward from %s to %s"%(larger_ele, s1, s2)
                else:
                    print "Significant increase in %s inward from %s to %s"%(larger_ele, s1, s2)
    return []   
    
def cluster_sizes(allSig_mask):
    """
    Find the sizes of each cluster.
    """
    
    labels, numL = label(allSig_mask)
    clust_arr = np.arange(numL)
    clust_sizes = np.zeros(len(clust_arr))
    for clust_idx, clust_num in enumerate(clust_arr):
        clust_sizes[clust_idx] = len(np.where(np.ravel(labels) == clust_num)[0])
        
    return clust_sizes, labels, numL
    
def clust_sigs(clust_num, separate_clusters, sig_map):
    """
    Find the signatures in each of the clusters.
    """
    clust_sig_list = []
    for clust_idx in np.arange(1, clust_num): # Don't include the background so start at 1
        clust_sig_map = sig_map[np.where(separate_clusters == clust_idx)]
        clust_sig_list.append(np.unique(clust_sig_map))
        
    return clust_sig_list
       
def classify_clusters(sig_map, separate_clusters, clust_num):
    """
    Code to automatically and empirically classify clusters.
    """
    # First find all the signatures within a cluster
    clust_sig_list = clust_sigs(clust_num, separate_clusters, sig_map)
    
    # Compare the signatures from each cluster to each other to create classes
    # of clusters with similar signatures.
    count = 0
    clust_class_list = []
    # For tracking if a cluster is already classified, non-zero if not already classified
    clust_track = np.arange(len(clust_sig_list))
    for sigs1_idx, sigs1 in enumerate(clust_track):
        if (sigs1 != 0) | (sigs1_idx == 0):
            for sigs2 in clust_track[(sigs1 + 1):]:
                if sigs2 != 0:                       
                    if np.shape(clust_sig_list[sigs1]) == np.shape(clust_sig_list[sigs2]):
                        # For equal number of signatures in each cluster:
                        # If they have the exact same signatures, put them in the same class
                        # and update the tracking array
                        if np.all(clust_sig_list[sigs1] == clust_sig_list[sigs2]):
                            if count == 0:
                                clust_class_arr = np.array([sigs1, sigs2])
                                count = count + 1
                            else:
                                clust_class_arr = np.concatenate((clust_class_arr, np.array([sigs2])))
                            clust_track[sigs1] = 0
                            clust_track[sigs2] = 0
                    elif abs(len(clust_sig_list[sigs1]) - len(clust_sig_list[sigs2])) == 1:
                        # For unequal number of sigs in each cluster:
                        # Find which cluster has more signatures.
                        if len(clust_sig_list[sigs1]) - len(clust_sig_list[sigs2]):
                            less_sigs = sigs1
                            more_sigs = sigs2
                        else:
                            less_sigs = sigs2
                            more_sigs = sigs1
                        
                        # Check to see if the signature with less signatures share all its signatures
                        # With the cluster with more signatures:
                        sig_compare = np.zeros(len(clust_sig_list[less_sigs]))
                        for sig_idx, this_sig in enumerate(clust_sig_list[less_sigs]):
                            if this_sig in clust_sig_list[more_sigs]:
                                sig_compare[sig_idx] = 1
                        
                        # Add to the same class if they do, and update tracking array
                        if all(sig_compare):
                            if count == 0:
                                clust_class_arr = np.array([sigs1, sigs2])
                                count = count + 1
                            else:
                                clust_class_arr = np.concatenate((clust_class_arr, np.array([sigs2])))
                            
                            clust_track[sigs1] = 0
                            clust_track[sigs2] = 0
                            
            count = 0
            if clust_class_list == []:
                clust_class_list.append(clust_class_arr)
            elif clust_class_list[len(clust_class_list)-1] is not clust_class_arr:
                clust_class_list.append(clust_class_arr)
            
    # Find the homeless clusters:
    homeless_clust = np.where(clust_track)
    
    return clust_class_list, homeless_clust, clust_sig_list
    
def clust_diff(sig1, sig2, lin_data, xLen=491, yLen=673, sz_diff=None):
    """
    A metric for the difference between the signatures found from two initial seed
    locations.
    """
    all_sig1 = bsu.all_sig_mask(sig1)
    all_sig2 = bsu.all_sig_mask(sig2)

    clust_sizes1, labels1, numL1 = cluster_sizes(all_sig1)
    clust_sizes2, labels2, numL2 = cluster_sizes(all_sig2)
    
    sig_map1 = bsu.signature_map(sig1, xLen, yLen)
    sig_map2 = bsu.signature_map(sig2, xLen, yLen)
    
    # First find all the signatures within a cluster
    clust_sig_lists = clust_sigs(numL1, labels1, sig_map1), clust_sigs(numL2, labels2, sig_map2)]
    
    clust_diffs_list = []
    clust_diffs = np.zeros(len(clust_sig_lists[0]))
    for clust in np.arange(1, len(clust_sig_lists[0])):
        # Make new signature maps containing just the signatures from the desired cluster.
        this_sig_map1 = []
        this_sig_map2 = []
        for sig in clust_sig_lists[0][clust]:
            this_sig_map1.append(sig1[sig-1])
            
        for sig in clust_sig_lists[1][clust]:
            this_sig_map2.append(sig2[sig-1])
        
        # Find the corresponding signatures in each image
        cc_arr, sim_sig_arr, less_idx = sig_reliability(lin_data, this_sig_map1, this_sig_map2, sz_diff=sz_diff)
        
        if less_idx == 0:
            more_idx = 1
        else:
            more_idx = 0
        
        clust_sig_diff = np.zeros(len(sim_sig_arr))
        for l_idx, sim_idx in enumerate(sim_sig_arr):
            clust_maps = [sig_map1[np.where(labels1 == clust)], sig_map2[np.where(labels2 == clust)]]
            less_sig_amnt = len(np.where(clust_maps[less_idx] == clust_sig_lists[less_idx][clust][l_idx])[0])
            sim_sig_amnt = len(np.where(clust_maps[more_idx] == clust_sig_lists[more_idx][clust][sim_idx])[0])
            
            clust_sig_diff[l_idx] = abs(less_sig_amnt - sim_sig_amnt)
        
        clust_diffs[clust] = np.sum(clust_sig_diff)
        clust_diffs_list.append(clust_sig_diff)
        
    return clust_diffs, clust_diffs_list