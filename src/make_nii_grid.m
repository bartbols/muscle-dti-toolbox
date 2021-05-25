function [X,Y,Z] = make_nii_grid(nii)
%MAKE_NII_GRID Returns the X,Y and Z coordinates (in physical coordinaes)
% of the NIfTI image.
% USAGE:
% [X,Y,Z] = make_nii_grid(nii)
%
% INPUT
% nii: nifti structure of nifti filename
%
% OUTPUT:
% [X,Y,Z] : arrays containing the x, y and z-coordinates of voxels in the
%           input image.
%
% Bart Bolsterlee
% Neuroscience Research Australia
% May 2021
%

if ~isstruct(nii)
    % Filename is provided. Load the NIfTI file.
    nii = load_untouch_nii(nii);
end

% Make grid with voxel indices
[Xv,Yv,Zv] = meshgrid(...
    0:nii.hdr.dime.dim(2)-1,...
    0:nii.hdr.dime.dim(3)-1,...
    0:nii.hdr.dime.dim(4)-1);

% Transform to physical coordinates using the transformation stored in the
% NIfTI header.
T = makeT_from_quat(nii);

X = T(1,1)*Yv + T(1,2)*Xv + T(1,3)*Zv + T(1,4);
Y = T(2,1)*Yv + T(2,2)*Xv + T(2,3)*Zv + T(2,4);
Z = T(3,1)*Yv + T(3,2)*Xv + T(3,3)*Zv + T(3,4);

end

