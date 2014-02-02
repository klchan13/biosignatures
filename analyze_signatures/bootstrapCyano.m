function [matSig_pairs, p_array] = bootstrapCyano(linearData, allSig, k)
% Kimberly Chan
% Last edited 6/22/13
%
% This program boostraps the different signatures to confirm that
% differences between the groups are statistically significant
%
% Parameters
% ----------
% linearData: 1 dimensional array
%       Elemental data in linear form
% allSig: cell array
%       Binary mask of different signatures
% k: int
%       Number of bootstrap samples
% Returns
% -------
% matSig_pairs: 1 dimensional array
%       Pairs of signatures to be evaluated by boostrapping
% p_array: 2 dimensional array
%       An array with all the different p values between signatures at each
%       pair of elements.

% Find all possible pairs of signatures using n choose k
numSigs = 1:length(allSig);
matSig_pairs = nchoosek(numSigs, 2);

% Find size of linear data to get number of elements
szEleNum = size(linearData);
idx_eleNum = 1:szEleNum(1);

p_array = [];
for sp = 1:length(matSig_pairs)
    tempArray = [];
    for en = idx_eleNum
        % Find the populations corresponding to each of the signatures
        pop1 = linearData(en,find(allSig{matSig_pairs(sp, 1)}));
        pop2 = linearData(en,find(allSig{matSig_pairs(sp, 2)}));
        
        % Draw random samples from the larger population of the two
        % populations so that the two populations have the same amount
        % of samples.
        if length(pop1) > length(pop2)
            sa = randsample(pop1,length(pop2),1);
            sb = pop2;
        else
            sa = pop1;
            sb = randsample(pop2,length(pop1),1);
        end
        
        % Bootstrap k samples and find the median values
        m1 = bootstrp(k, @median, [pop1, pop2]);
        m2 = bootstrp(k, @median, [pop1, pop2]);
        fprintf('\rBootstrapping element %d of pair %d of %d pairs.\r',[en, sp, length(matSig_pairs)])
        
        % Find the actual mean and find p value of the difference in means
        d_empirical = median(sa) - median(sb);
        p = sum(abs(m1 - m2) >= abs(d_empirical))/k;
        tempArray = [tempArray, p];
    end
    p_array = [p_array; tempArray];
end
end
