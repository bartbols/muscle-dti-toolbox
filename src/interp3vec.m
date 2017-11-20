function [ v ] = interp3vec(V, coords )
%%INTERP3VEC linearly interpolates the vectorfield V (4-D matrix, with
% x-y-z components stored in 4th dimension at the locations provided in
% 'coords' (3 x n matrix with voxel coordinates). The interpolation is done
% such that the direction of the vector does not matter, because it doesn't
% matter which direction the primary eigenvector points towards. The output
% is the normalised direction at the locations provided in 'coords'.
%
% Bart Bolsterlee, Neuroscience Research Australia
% November 2017

% Define some variables
nP        = size(coords,2);
v         = zeros(3,nP);
dim       = size(V);
vox_idx   = zeros(3,2);
W         = zeros(2,2,2);

% Add 1 to account for indexing from 0 in the input and indexing from 1 in
% MATLAB.
coords = coords + 1; 
for i = 1 : nP
    
    % Get the voxel index of the current points and the adjacent
    % voxels that are used for interpolation
    
    vox_idx(1:3,1) = round(coords(:,i));
%     vox_idx(1:3,2) = vox_idx(:,1) + sign(coords(:,i) - round(coords(:,i)));
    vox_idx(1:3,2) = vox_idx(:,1) + sign(coords(:,i) - round(coords(:,i)));
    
    
    if any(vox_idx(:) < 1) ||...
       any(vox_idx(:,1)' > dim([2 1 3])) ||...
       any(vox_idx(:,2)' > dim([2 1 3]))
        % outside field
        continue
    end

    
    % Make the vectors point in the same direction by checking the
    % sign of the inner product between the direction vector in the
    % voxel nearest to the current point and the direction vectors
    % of all neighbouring voxels. Flip if sign is negative.
    sgn = sign(repmat(V(vox_idx(2,1),vox_idx(1,1),vox_idx(3,1),1),2,2,2) .*...
                      V(vox_idx(2,:),vox_idx(1,:),vox_idx(3,:),1) +...
               repmat(V(vox_idx(2,1),vox_idx(1,1),vox_idx(3,1),2),2,2,2) .*...
                      V(vox_idx(2,:),vox_idx(1,:),vox_idx(3,:),2) +...
               repmat(V(vox_idx(2,1),vox_idx(1,1),vox_idx(3,1),3),2,2,2) .*...
                      V(vox_idx(2,:),vox_idx(1,:),vox_idx(3,:),3) );

    % If adjacent vectors are the same, the projection and the sign
    % are zero. In that case, keep the original direction by
    % setting 'sgn' to 1.
    sgn(sgn==0)=1;
    
    % Flip the sign of the vectors in the neighbouring voxels that points 
    % in the opposite  direction (= negative projection) to the vector in 
    % the nearest voxel
    X = V(vox_idx(1,:),vox_idx(2,:),vox_idx(3,:),1).* sgn;
    Y = V(vox_idx(1,:),vox_idx(2,:),vox_idx(3,:),2).* sgn;
    Z = V(vox_idx(1,:),vox_idx(2,:),vox_idx(3,:),3).* sgn;
    
    % Weight-matrix for interpolation. The weight for a voxel is
    % based on the distance of the centre of that voxel to the
    % point that is interpolated, so that points that are near the
    % point are more heavily weighed than points that are far away.
    % (This is just trilinear interpolation.)
    r = abs(coords(:,i) - vox_idx(:,1));
    
    W(1,1,1) = (1-r(2)) * (1-r(1)) * (1-r(3));
    W(1,2,1) =    r(2)  * (1-r(1)) * (1-r(3));
    W(2,1,1) = (1-r(2)) *    r(1)  * (1-r(3));
    W(2,2,1) =    r(2)  *    r(1)  * (1-r(3));
    W(1,1,2) = (1-r(2)) * (1-r(1)) *    r(3);
    W(1,2,2) =    r(2)  * (1-r(1)) *    r(3);
    W(2,1,2) = (1-r(2)) *    r(1)  *    r(3);
    W(2,2,2) =    r(2)  *    r(1)  *    r(3);
    
    % Interpolate
    v(1,i) = sum((X(:) .* W(:))) / sum(W(:));
    v(2,i) = sum((Y(:) .* W(:))) / sum(W(:));
    v(3,i) = sum((Z(:) .* W(:))) / sum(W(:));
    
    
%     A(i) = all(sgn(:)>0);
end

% make unit vector
v = v ./ (ones(3,1) * sqrt(v(1,:).^2 + v(2,:).^2 + v(3,:).^2));
v(:,all(isnan(v),1)) = 0;

% % This code was used for checking the interpolation. When all signs of
% neighbouring voxels are the same, the results should coincide with the
% results from interp3
% v2(1,:) = interp3(V(:,:,:,1),coords(2,:),coords(1,:),coords(3,:));
% v2(2,:) = interp3(V(:,:,:,2),coords(2,:),coords(1,:),coords(3,:));
% v2(3,:) = interp3(V(:,:,:,3),coords(2,:),coords(1,:),coords(3,:));
% v2 = v2 ./ (ones(3,1) * sqrt(v2(1,:).^2 + v2(2,:).^2 + v2(3,:).^2));
% [v;v2;A]
end


