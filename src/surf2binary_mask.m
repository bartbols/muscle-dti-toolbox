function surf2binary_mask(nii,surface,mask_filename)
%%SURF2BINARY_MASK Creates a binary mask with all points inside the surface
% of the dimensions provided in NIfTI image 'nii'. The resulting mask is 
% saved as a NIfTI filewith filename 'mask_filename'.
%
% Bart Bolsterlee, Neuroscience Research Australia
% November 2017
%
% Uses the function 'VOXELISE' by Adam A, downloaded from MATLAB central on
% 1/11/2017:
% http://au.mathworks.com/matlabcentral/fileexchange/27390-mesh-voxelisation
%
% INPUT:
% nii           - filename of NIfTI file or structure with NIfTI data incl 
%                 header information
% surface       - filename of STL-surface or a structure containing
%                 faces/vertices of the surface model
% mask_filename - filename of the mask (with extension .nii.gz) that will
%                 be created


if ischar(surface)
    % Load STL file
    FV = stlread(surface);
else
    FV = surface;
end
clear surface

if ischar(nii)
    % Load image file
    img = load_untouch_nii(nii);
else
    img = nii;
end
clear nii

% Convert all vertices to voxel dimensions
vertices = FV.vertices;
T = [img.hdr.hist.srow_x;...
     img.hdr.hist.srow_y;...
     img.hdr.hist.srow_z;...
     0 0 0 1];
voxel_coords = [vertices ones(size(vertices,1),1)] / (T');

% Create a patch object in voxel dimensions.
FV_voxel.vertices = voxel_coords(:,1:3);
FV_voxel.faces = FV.faces;

% Use VOXELISE to make a binary mask
OUTPUTgrid = VOXELISE(0:img.hdr.dime.dim(2)-1,...
                      0:img.hdr.dime.dim(3)-1,...
                      0:img.hdr.dime.dim(4)-1,...
                      FV_voxel,'xyz');

% Save the mask
mask = img;
mask.img = cast(OUTPUTgrid,'like',mask.img);
mask.hdr.dime.scl_slope = 1;
mask.hdr.dime.scl_inter = 0;
save_untouch_nii(mask,mask_filename)
fprintf('binary mask written to %s\n',mask_filename)
