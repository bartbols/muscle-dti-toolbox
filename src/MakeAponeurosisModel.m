function varargout = MakeAponeurosisModel( mask,anat,model_name )
%MAKEAPONEUROSISMODEL This function makes a surface model of the
%aponeurosis.
%
% Bart Bolsterlee
% Neuroscience Research Australia
% November 2017
%
% ----------------- USAGE -----------------
% model = MakeAponeurosisModel( mask,anat,model_name );
% 
% ----------------- INPUT -----------------
% - mask    :  NIfTI structure containing the mask of the aponeurosis as
%              non-zero entries. A filename of a NIfTI file may also be
%              provided.
% - anat    :  NIfTI structure containing the anatomical data on whcih the 
%              aponeurosis was outlined. A filename of a NIfTI file may 
%              also be provided.
% - model_name : filename (with extension .stl) to which the STL-model will
%                be saved
%
%         
% ----------------- OUTPUT -----------------
% - model: structure with the vertices and faces of the model


% Read inputs
p = inputParser;
addRequired(p,'mask')
addRequired(p,'anat')
addRequired(p,'model_name',@(x) endsWith(x,'.stl','IgnoreCase',true))
parse(p, mask,anat,model_name)

%% Check if structure or filename is provided
if ~isstruct(mask)
    mask = load_untouch_nii(mask);
end
if ~isstruct(anat)
    anat = load_untouch_nii(anat);
end

% Calculate point cloud with all voxel locations
[I,J,K] = ind2sub(mask.hdr.dime.dim(2:4),find(mask.img~=0));

% Transform to global coordinates using the transformation in the header of
% the anatomical image. This information is not always stored correctly in
% the mask file (e.g. when the mask is created in 3D Slicer), so that's why
% the anatomical image is used.

T = [anat.hdr.hist.srow_x;anat.hdr.hist.srow_y;anat.hdr.hist.srow_z;0 0 0 1];
coords = [I-1 J-1 K-1 ones(size(I))] * T';
model.vertices = coords(:,1:3);

% Use MATLAB's built-in 'boundary' function to create a closed surface
% around the point cloud.
model.faces = boundary(model.vertices);

% Calculate normal vectors per face
model.normals = facenormals(model);

% Write the model
stlwrite(model_name,model)
fprintf('Aponeurosis model written to %s.\n',model_name)
if nargout == 1
    varargout{1} = model;
end


end

