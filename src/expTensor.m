function tensor= expTensor(log_vec)
%EXPTENSOR calculates the inverse of the Log-Euclidean transformation (i.e.
% exp(log(S)) of a diffusion tensor following procedures described in: 
%
% Arsigny, V., Fillard, P., Pennec, X., Ayache, N., 2006. Log-Euclidean
% metrics for fast and simple calculus on diffusion tensors. Magn Reson Med
% 56 (2),p 411-421.http://doi.org/10.1002/mrm.20965
%
% USAGE: tensor = expTensor(log_vec)
% INPUT:
% log_vec: nx6 vectors of Log-Euclidean tensor in vector representation.
% See logTensor.m for details of the format.
% OUTPUT:
% tensor : 1D, 2D 3D tensor field. The last two dimensions will be 3x3,
%     containing the tensor for that voxel. Array dimensions will be 2D,
%     3D, 4D or 5D for a 0D (single tensor), 1D, 2D or 3D tensor field,
%     respectively.
%
% Bart Bolsterlee
% Neuroscience Research Australia
% May 2021

% Get input dimensions.
nd = ndims(log_vec);
sz = size(log_vec);

% Reshape vector field into one long array of vectors for easy processing.
n = prod(sz(1:nd-1));
log_vec = reshape(log_vec,[n,6]);

tensor = NaN(n,3,3);
for i = 1 : n
    % Put vector back in tensor format
    logS = [log_vec(i,1) log_vec(i,4)/sqrt(2) log_vec(i,5)/sqrt(2);...
            log_vec(i,4)/sqrt(2) log_vec(i,2) log_vec(i,6)/sqrt(2);...
            log_vec(i,5)/sqrt(2) log_vec(i,6)/sqrt(2) log_vec(i,3)];
        
    if any(isnan(logS(:)))
        tensor(i,:,:) = NaN;
        continue
    end
    % Extract eigenvalues and rotation matrix
    [R,Dstar] = eig(logS);
    
    % Transform back using the exponential operation
    D = diag(exp(diag(Dstar)));
    
    % Rotate back
    tensor(i,:,:) = R * D * R';
end

tensor = reshape(tensor,[sz(1:nd-1),3,3]);

end

