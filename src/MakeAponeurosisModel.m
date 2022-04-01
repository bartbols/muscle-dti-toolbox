function varargout = MakeAponeurosisModel( mask,model_name,varargin )
%MAKEAPONEUROSISMODEL This function makes a surface model of the
%aponeurosis.
%
% Bart Bolsterlee
% Neuroscience Research Australia
% November 2017
%
% ----------------- USAGE -----------------
% model = MakeAponeurosisModel( mask,model_name );
% model = MakeAponeurosisModel( mask,model_name,'label_number',1);
% 
% ----------------- INPUT -----------------
% - mask       :  Filename of a NIfTI file containing the mask.
% - model_name : filename (with extension .stl) to which the STL-model will
%                be saved
%
% Optional input, provided as 'parameter',<value> pairs:
% - shrink       : shrink factor (see documentation of function 'boundary'
%                  for details. Default: 1
% - label_number : number of the label in the mask from which the
%                  aponeurosis will be created. If not provided all 
%                  non-zero voxels are used.
%         
% ----------------- OUTPUT -----------------
% - model: structure with the vertices and faces of the model

% Read inputs
p = inputParser;
addRequired(p,'mask')
addRequired(p,'model_name',@(x) endsWith(x,'.stl','IgnoreCase',true))
addParameter(p,'label_number',[],@(x) isscalar(x))
addParameter(p,'shrink',1,@(x) isscalar(x) && x>=0 && x<=1)
parse(p, mask,model_name,varargin{:})
label_number = p.Results.label_number;

%%
% If multiple labels are present in the mask, make a new mask with only the
% selected label in it.
tic
M = load_untouch_nii(mask);
if isempty(label_number)
    fprintf('Creating aponeurosis model of all non-zero labels in %s...\n',mask)
    % Set all non-zero voxels to 1.
    M.img = cast(M.img ~= 0,'like',M.img);
else 
    fprintf('Creating aponeurosis model of label %d in %s...\n',label_number,mask)
    % Select voxels within the selected label.
    M.img = cast(M.img == label_number,'like',M.img);
end
% save_untouch_nii(M,fullfile(tempdir,'mask.nii.gz'))    

% % Resample to isotropic dimensions to create a smooth surface with
% % faces that have more or less equally-sized edges (aspect ratio of 
% % triangles close to 1).
% commandTxt = sprintf('c3d -int 3 %s -resample-mm 0.5x0.5x0.5mm -o %s',...
%         fullfile(tempdir,'mask.nii.gz'),fullfile(tempdir,'mask_iso.nii.gz'));
% system(commandTxt);
% mask = load_untouch_nii(fullfile(tempdir,'mask_iso.nii.gz'));

% % Clean up files.
% delete(fullfile(tempdir,'mask_iso.nii.gz'))
% delete(fullfile(tempdir,'mask.nii.gz'))

% Calculate point cloud with all voxel locations
% if isempty(label_number)
%     [I,J,K] = ind2sub(mask.hdr.dime.dim(2:4),find(mask.img~=0));
%     [I,J,K] = ind2sub(mask.hdr.dime.dim(2:4),find(mask.img>=0.5));
[I,J,K] = ind2sub(M.hdr.dime.dim(2:4),find(M.img~=0));

% Transform to global coordinates using the transformation in the header of
% the mask image.

T = [M.hdr.hist.srow_x;M.hdr.hist.srow_y;M.hdr.hist.srow_z;0 0 0 1];
if all(all(T(1:3,1:3)==0)) 
    % The srow information is missing from the header.
    % Calculate the spatial transformation matrix from the
    % quaternion parameters.
    T = makeT_from_quat( M );
end

coords = [I-1 J-1 K-1 ones(size(I))] * T';
model.vertices = coords(:,1:3);

% Use MATLAB's built-in 'boundary' function to create a closed surface
% around the point cloud.
model.faces = boundary(model.vertices,p.Results.shrink);

% Remesh to a more regular distribution of triangles.
opt.gridsize = 0.25;opt.closesize = 0.5;opt.elemsize = 0.5;
[model.vertices,model.faces] = remeshsurf(model.vertices,model.faces,opt);
[model.vertices,model.faces] = surfreorient(model.vertices,model.faces);

% Apply some smoothing
% unsmoothed = model;
conn = meshconn(model.faces,size(model.vertices,1));
model.vertices = smoothsurf(model.vertices,[],conn,5,0.7,'lowpass');

% Calculate normal vectors per face
model.normals = facenormals(model);

fprintf(' completed.\n')
% Write the model
% stlwrite(strrep(model_name,'.stl','_unsmoothed.stl'),unsmoothed)
if exist(fileparts(model_name),'dir') ~=7;mkdir(fileparts(model_name));end
stlwrite2(model_name,model)
fprintf('Aponeurosis model written to %s.\n',model_name)

% plot_surface(model,'ShowNormals',true,'FaceColor','y')'
if nargout == 1
    varargout{1} = model;
end
t_elapsed = toc;
fprintf('It took %.2f seconds to create all models.\n',t_elapsed)

end

