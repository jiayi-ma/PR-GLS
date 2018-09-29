function s_matrix = ICG_DistanceMatrixNormalization(distance_matrix, kth)
% ICG_DISTANCEMATRIXNORMALIZATION normalizes a distance matrix (the higher 
%   the more different) by the method proposed in "Self-tuning spectral clustering" 
%   from Zelnik-Manor and Perona, in NIPS 2004
%
%   Parameter
%   ---------
%   distance_matrix ... The input matrix (the higher the more different)
%   kth ... The kth nearest neigbor to be considered
%
%   Returns
%   -------
%   s_matrix ... Normalized affinity matrix (values between 0 (differnt) and 1 (similar))
%
%   Example
%   -------
%   s_matrix = ICG_DistanceMatrixNormalization(distance_matrix, 8)

    % sort the distances and choose the kth nearest neigbor for normalization
    sA = sort(distance_matrix,2);
    sigmas_i = sA(:, kth + 1);
    W = sigmas_i*sigmas_i';

    s_matrix = exp(-distance_matrix.^2./W);