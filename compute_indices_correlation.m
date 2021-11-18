function cellind = compute_indices_correlation(X, y, Indperm, type, leftThresh, rightThresh)
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

if isempty(type)
    type = 'classical';
elseif ~strcmp(type, 'classical') && ~strcmp(type, 'Spearman')
    error('The type of correlation stat must be either ''classical'' or ''Spearman''.');
end
if ~isvector(leftThresh) || (numel(leftThresh) ~= numel(rightThresh))
    error('The left and right threshold inputs should be vectors containing the same number of elements.');
end
if (max(leftThresh) > min(rightThresh))
    error('The left and right threshold inputs allow to select the indices for the features for which the statistics produced values in the lower and upper tail, respectively.');
end
M = size(X,1);
if (numel(y) ~= M) || (size(y,2) ~= 1)
    error('y must be a column vector with as many elements and number of rows in X.');
end
if (min(min(Indperm)) < 1)
    error('Permutation indices should start at 1, following Matlab notation.');
end
if (max(max(Indperm)) > M)
    error('Permutation indices should not exceed %d, following Matlab notation.', M);
end
Indperm = int32(Indperm)-1; % now integers start in 0 following C notation

switch type
    case 'classical'
        cellind = compute_indices_correlation_cppopt(X, y, Indperm, int32(1), [leftThresh;rightThresh]);
    case 'Spearman'
        [Xs,ix] = sort(X);
        ix = int32(ix)-1;
        [ys,iy] = sort(y);
        iy = int32(iy)-1;
        cellind = compute_indices_correlation_cppopt(Xs, ys, Indperm, int32(2), [leftThresh;rightThresh], ix, iy);
end
cellind = reshape(cellind, size(Indperm,2), 2, []);