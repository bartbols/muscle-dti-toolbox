function varargout = MakeAponeurosisModel( mask,anat,model_name,varargin )
%MAKEAPONEUROSISMODEL This function makes a surface model of the
%aponeurosis.
%
% Bart Bolsterlee
% Neuroscience Research Australia
% November 2017
%
% ----------------- USAGE -----------------
% model = MakeAponeurosisModel( mask,anat,model_name );
% model = MakeAponeurosisModel( mask,anat,model_name,'label_number',1);
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
% Optional input, provided as 'parameter',<value> pairs:
% - label_number : number of the label in the mask from which the
%    aponeurosis will be created. If not provided all non-zero voxels are
%    used.
%         
% ----------------- OUTPUT -----------------
% - model: structure with the vertices and faces of the model


% Read inputs
p = inputParser;
addRequired(p,'mask')
addRequired(p,'anat')
addRequired(p,'model_name',@(x) endsWith(x,'.stl','IgnoreCase',true))
addParameter(p,'label_number',[],@(x) isscalar(x))
parse(p, mask,anat,model_name,varargin{:})
label_number = p.Results.label_number;

%% Check if structure or filename is provided
if ~isstruct(mask)
    mask = load_untouch_nii(mask);
end
if ~isstruct(anat)
    anat = load_untouch_nii(anat);
end

% Calculate point cloud with all voxel locations
if isempty(label_number)
    [I,J,K] = ind2sub(mask.hdr.dime.dim(2:4),find(mask.img~=0));
else
    [I,J,K] = ind2sub(mask.hdr.dime.dim(2:4),find(mask.img==label_number));
end

% Transform to global coordinates using the transformation in the header of
% the anatomical image. This information is not always stored correctly in
% the mask file (e.g. when the mask is created in 3D Slicer), so that's why
% the anatomical image is used.

T = [anat.hdr.hist.srow_x;anat.hdr.hist.srow_y;anat.hdr.hist.srow_z;0 0 0 1];
coords = [I-1 J-1 K-1 ones(size(I))] * T';
model.vertices = coords(:,1:3);

% Use MATLAB's built-in 'boundary' function to create a closed surface
% around the point cloud.
model.faces = boundary(model.vertices,1);

% Calculate normal vectors per face
model.normals = facenormals(model);

% Write the model
stlwrite(model_name,model)
fprintf('Aponeurosis model written to %s.\n',model_name)
if nargout == 1
    varargout{1} = model;
end


end

