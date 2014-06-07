import numpy as np
import csv
from scipy.ndimage.measurements import label

def cluster_signatures(sig_out_in, bstrp_med, med_sig, p_thresh=0.05,
                       ele_order=np.array(["Fe", "Cu", "Zn", "Ca", "K", "S", "P", "Cl", "Si", "Mn"])):
    """
    Analyzes different signatures within a cluster group and finds the elemental trends.
    """
    for s_idx, sig in enumerate(sig_out_in):
        if s_idx != len(sig_out_in)-1:
            # Grab the two signatures of interst
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