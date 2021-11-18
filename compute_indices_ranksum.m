function cellind = compute_indices_ranksum(A, B, Indperm, leftThresh, rightThresh)
%
% INPUT:
% leftThresh, rightThresh - paired threshold values used to select the
%                           indices of significant columns with stat values
%                           in the lower/upper distribution tails.
%
% OUTPUT:
% cellind - a cell array containing the indices of significant columns for
%           the ranksum statistics according the input thresholds.
%           cellind is Nr x 2 x N, where Nr is the number of Monte Carlo
%           simulations or resampling and N is the number of 

if ~isvector(leftThresh) || (numel(leftThresh) ~= numel(rightThresh))
    error('The left and right threshold inputs should be vectors containing the same number of elements.');
end
if (max(leftThresh) > min(rightThresh))
    error('The left and right threshold inputs allow to select the indices for the features for which the statistics produced values in the lower and upper tail, respectively.');
end
[mA,N] = size(A);
if (size(B,2) ~= N)
    error('A and B matrices must have the same number of columns.');
end
mB = size(B,1);
M = mA + mB;
if (min(min(Indperm)) < 1)
    error('Permutation indices should start at 1, following Matlab notation.');
end
if (max(max(Indperm)) > M)
    error('Permutation indices should not exceed %d, following Matlab notation.', M);
end
Indperm = int32(Indperm)-1; % now integers start in 0 following C notation
X = [A; B];
[Xs,ii] = sort(X);
ii = int32(ii)-1;
cellind = compute_indices_ranksum_cppopt(Xs, ii, Indperm, int32(mA), [leftThresh;rightThresh]);
cellind = reshape(cellind, size(Indperm,2), 2, []);