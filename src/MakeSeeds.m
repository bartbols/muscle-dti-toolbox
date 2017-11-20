function seeds = MakeSeeds( mask,img,spacing_mm )
%MAKESEEDS Creates seeds within a mask with regular spacing and returns the
%seed location in global coordinates. These seed locations can be used in
%dti_tractography
%
% Bart Bolsterlee, Neuroscience Research Australia (NeuRA)
% February 2017
%
% ----------------- USAGE -----------------
% seeds = MakeSeeds( mask,DTI,spacing_mm )
%
% ----------------- INPUT -----------------
% mask : NIfTI structure or filename of NIfTI file with the seed mask.
%        (Usually this will be a mask containing the whole muscle.)
% img  : NIfTI structure or filename of NIfTI file with the corresponding
%        DTI/anatomical data. The header information is required to convert
%        the seeds to the correct global coordinate system. The dimensions
%        of 'img' and 'mask' should be the same, i.e. the mask was created
%        on this image.
% spacing_mm : spacing of the regular grid of seed points in mm

% Load the mask
if ~isstruct(mask)
    mask = load_untouch_nii(mask);
end
M = mask.img;

% Load the DTI data
if ~isstruct(img)
    img = load_untouch_nii(img);
end

% Get all voxels locations with the mask
% Calculate point cloud with all voxel locations
[I,J,K] = ind2sub(mask.hdr.dime.dim(2:4),find(mask.img~=0));

% Transform to global coordinates using the transformation in the header of
% img. This information is not always stored correctly in
% the mask file (e.g. when the mask is created in 3D Slicer).

T = [img.hdr.hist.srow_x;...
     img.hdr.hist.srow_y;...
     img.hdr.hist.srow_z;...
     0 0 0 1];

coords_glob = [I-1 J-1 K-1 ones(size(I))] * T';
coords_glob = coords_glob(:,1:3)';

% Make a regular grid between the boundaries of the mask
[Sx, Sy, Sz] = meshgrid(min(coords_glob(1,:)) : spacing_mm : max(coords_glob(1,:)),...
                        min(coords_glob(2,:)) : spacing_mm : max(coords_glob(2,:)),...
                        min(coords_glob(3,:)) : spacing_mm : max(coords_glob(3,:)));

seeds_glob = [Sx(:) Sy(:) Sz(:)]';

% Determine which seeds are within the mask
seeds_vox = T \ [seeds_glob;ones(1,size(seeds_glob,2))];

in_mask = interp3(mask.img,...
    seeds_vox(2,:)+1,...
    seeds_vox(1,:)+1,...
    seeds_vox(3,:)+1,...
    'nearest')~=0;
seeds = seeds_glob(:,in_mask);

end

