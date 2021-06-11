function [l1,l2,l3,md] = tensor2lambda(tensor)
%TENSOR2LAMBDA calculates the eigenvalues from a diffusion tensor
%field.
%
% USAGE: [l1,l2,l3] = tensor2lambda(tensor)
% INPUT: 
% tensor : 1D, 2D 3D tensor field. The last two dimensions should be 3x3,
%     containing the tensor for that voxel. Array dimensions will be 2D,
%     3D, 4D or 5D for a 0D (single tensor), 1D, 2D or 3D tensor field,
%     respectively.
% 
% OUTPUT:
% lambda1: (i,j,k)x3 field with the primary eigenvalue. 
% lambda2: (i,j,k)x3 field with the secondary eigenvalue. 
% lambda3: (i,j,k)x3 field with the tertiary eigenvalue. 
% lambda3: (i,j,k)x3 field with mean diffusivity.
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

l1 = NaN(n,3);
l2 = NaN(n,3);
l3 = NaN(n,3);
md = NaN(n,3);
for i = 1 : n
%     if mod(i,100)==0
%         fprintf('%d of %d\n',i,n)
%     end
    S = squeeze(tensor(i,:,:));
    if any(isnan(S(:)))
        continue
    end
    
    % Decompose into rotation matrix v and eigenvalue matrix l.
    [~,l] = eig(S);
    
    if any(diag(l) <= 0)
        continue
    end    
    
    % Get primary eigenvector.
    l = sort(diag(l),'descend');
    l1(i) = l(1);
    l2(i) = l(2);
    l3(i) = l(3);
    md(i) = mean(l);
end

% Put tensor field back in original (input) dimensions.
if nd ==2
    l1 = squeeze(l1);
    l2 = squeeze(l2);
    l3 = squeeze(l3);
    md = squeeze(md);
else
    l1    = reshape(l1,[sz(1:nd-2) 3]);
    l2    = reshape(l2,[sz(1:nd-2) 3]);
    l3    = reshape(l3,[sz(1:nd-2) 3]);
    md    = reshape(md,[sz(1:nd-2) 3]);
end


end

