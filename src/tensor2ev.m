function [ev1,ev2,ev3] = tensor2ev(tensor)
%TENSOR2FA calculates the eigenvectors from a diffusion tensor
%field.
%
% USAGE: [ev1,ev2,ev3] = tensor2ev(tensor)
% INPUT: 
% tensor : 1D, 2D 3D tensor field. The last two dimensions should be 3x3,
%     containing the tensor for that voxel. Array dimensions will be 2D,
%     3D, 4D or 5D for a 0D (single tensor), 1D, 2D or 3D tensor field,
%     respectively.
% 
% OUTPUT:
% fa: (i,j,k)x3 field with the primary eigenvector. 
%
% Bart Bolsterlee
% Neuroscience Research Australia
% May 2021

% Get input dimensions.
nd = ndims(tensor);
sz = size(tensor);

% Reshape tensor field into one long array of tensors for easy processing.
n = prod(sz(1:nd-2));
tensor = reshape(tensor,[n,3,3]);

ev1 = NaN(n,3);
ev2 = NaN(n,3);
ev3 = NaN(n,3);
for i = 1 : n
%     if mod(i,100)==0
%         fprintf('%d of %d\n',i,n)
%     end
    S = squeeze(tensor(i,:,:));
    if any(isnan(S(:)))
        continue
    end
    
    % Decompose into rotation matrix v and eigenvalue matrix l.
    [v,l] = eig(S);
    
    if any(diag(l) <= 0)
        continue
    end    
    
    % Get primary eigenvector.
    [~,idx] = sort(diag(l),'descend');
    ev1(i,:) = v(:,idx(1))';
    ev2(i,:) = v(:,idx(2))';
    ev3(i,:) = v(:,idx(3))';
end

% Put tensor field back in original (input) dimensions.
if nd ==2
    ev1 = squeeze(ev1);
    ev2 = squeeze(ev2);
    ev3 = squeeze(ev3);
else
    ev1    = reshape(ev1,[sz(1:nd-2) 3]);
    ev2    = reshape(ev2,[sz(1:nd-2) 3]);
    ev3    = reshape(ev3,[sz(1:nd-2) 3]);
end


end

