function dti_reg = register_DTI_to_anat( dti,anat,parfile,dti_reg,varargin)
%REGISTER_DTI_TO_ANAT uses elastix to register the DTI image to the
%anatomical scan.
%
% Bart Bolsterlee, Neuroscience Research Australia (NeuRA)
% September 2017
%
% Note: this function  requires Convert3D and elastix to be installed on
% the pc and added to the path so that 'c3d' and 'elastix' are recognised
% as an external command.
% Convert3D can be downloaded here:
% http://www.itksnap.org/pmwiki/pmwiki.php?n=Downloads.C3D
%
% Elastix can be downloaded here:
% http://elastix.isi.uu.nl/download.php
%
% ----------------- USAGE -----------------
% filename = register_DTI_to_anat(dti,anat,parfile,dti_reg,varargin)
%
% ----------------- INPUT -----------------
% ----- REQUIRED -----
% dti:      filename of the DTI data file
% anat:     filename of the anatomical MRI data
% parfile:  filename (char) or a cell string of filenames of the elastix 
%           parameter file(s) with registration parameters. For example, to
%           perform a rigid and then a bspline registration: parfile =
%           {'rigid.txt','bspline.txt'}. Or to only perform a bspline, use
%           'bspline.txt'.
% dti_reg:  filename of the registered DTI file (will be created)
%
% ----- OPTIONAL -----
% Optional inputs are provided as pairs of 'ParameterName',<value> (e.g.
% 'mask','mask.nii.gz')
%
% - mask                 : filename of the mask file used for registration
% - foreground_threshold : threshold intensity for foreground. A foreground
%                          mask will be created using this threshold.
% - stack                : if anatomical scan is 4D, choose which stack is
%                          used for registration.
% - b0_stack             : stack numbers in DTI file to which the
%                          anatomical scan will be registered. Default = 1
%                          (usually the b0-image is the first image in the stack)
% - InspectRegistration  : opens ITK-snap with the anatomical scan and the
%                          DTI scans before and after registration.
%                          This requires ITK-SNAP to be added to the path.
%
% ----------------- OUTPUT -----------------
% dti_reg: filename of the registered DTI file
%
% Example without using a mask:
% dti_reg = register_DTI_to_anat('dti.nii.gz','anat.nii.gz','parameter.txt','dti_reg.nii.gz')
%
% Example with mask:
% dti_reg = register_DTI_to_anat('dti.nii.gz','anat.nii.gz',...
%                 'parameter.txt','dti_reg.nii.gz',...
%                 'mask','mask.nii.gz')
%
% Or to automatically create a foreground mask from the anatomical scan
% including all voxels with intensity larger than 20:
%
% dti_reg = register_DTI_to_anat('dti.nii.gz','anat.nii.gz',...
%                 'parameter.txt','dti_reg.nii.gz',...
%                 'foreground_threshold',20)

% Read the inputs
p = inputParser;
addRequired(p,'dti',@(x) contains(x,'.nii.gz'))
addRequired(p,'anat',@(x) contains(x,'.nii.gz'))
addRequired(p,'parfile',@(x) ischar(x) || iscell(x))
addRequired(p,'dti_reg',@(x) contains(x,'.nii.gz'))
addParameter(p,'mask',[],@(x) isempty(x) || contains(x,'.nii.gz'))
addParameter(p,'foreground_threshold',10,@(x) isscalar(x) || isempty(x))
addParameter(p,'stack',[],@(x) isscalar(x) || isempty(x))
addParameter(p,'b0_stack',1,@(x) isscalar(x) || isempty(x))
addParameter(p,'InspectRegistration',false,@(x) islogical(x) || x==1 || x==0);
parse(p,dti,anat,parfile,dti_reg,varargin{:});

% Assign parameter values
foreground_threshold = p.Results.foreground_threshold;
mask                 = p.Results.mask;
stack                = p.Results.stack;
b0_stack             = p.Results.b0_stack;
InspectFlag          = p.Results.InspectRegistration;

% Create temporary working directory.
char_list = char(['a':'z' '0':'9']) ;
tmpdir = [];
while exist(tmpdir,'dir') == 7 || isempty(tmpdir)
    tmpdir = fullfile(tempdir,char_list(ceil(length(char_list)*rand(1,8))));
end
mkdir(tmpdir)

try
    % Set some filenames
    b0_map         = fullfile(tmpdir,'b0_map.nii.gz');
    
    if ~isempty(stack)
        % Extract one stack from the anatomical scan if the data is 4D
        anat_3D = fullfile(tmpdir,'anat_3D.nii.gz');
        extract_3Dfrom4D(anat,anat_3D,stack);
    else
        anat_3D = anat;
    end
    
    % Create mask based on foreground_threshold if no mask is provided
    if isempty(mask)
        if isempty(foreground_threshold)
            foreground_threshold = -Inf;
        end
        anat_mask = load_untouch_nii(anat_3D);
        anat_mask.img = cast(anat_mask.img * anat_mask.hdr.dime.scl_slope>foreground_threshold,'like',anat_mask.img);
        anat_mask.hdr.dime.scl_slope = 1;
        mask = fullfile(tmpdir,'mask_3D.nii.gz');
        save_untouch_nii(anat_mask,mask);
    end
    
    % Load the DTI data
    dti_data = load_untouch_nii(dti);
    
    % Create an empty nifti-structure that will be filled with the
    % registered data.
    dti_reg_data = dti_data;
    dti_reg_data.img = cast(zeros(size(dti_data.img)),'like',dti_data.img);
    nStacks = dti_data.hdr.dime.dim(5);
    
    % Open a progress bar
    h = waitbar(0,'Registering DTI data...',...
        'Name','Progress bar DTI registration');
    waitbar(1 / (nStacks+2),h,'Resampling anatomical data and mask...')
    
    % Extract B0-map from 4D dti data. Elastix does not accept 4D
    % data, so a 3D image should be created before registration can be done.
    extract_3Dfrom4D(dti,b0_map,b0_stack);
    waitbar(2 / (nStacks+2),h,...
        'Registering the B0-map to the anatomical scan.')
    
    %% Register B0 map to the anatomical scan
    if ischar(parfile)
        % Only one parameter file is provided.
        parfile = {parfile};
    end
    
    % Check how many steps there are in the registration
    nSteps = length(parfile);
    % Now, iterate through the steps, using the transformation from the 
    % previous step as the initial transformation.
    for i = 1 : nSteps
        mkdir(fullfile(tmpdir,sprintf('step%02d',i)))
        if i == 1
            % First step, no initial transform.
            elastix_cmd = sprintf('elastix -f %s -m %s -fMask %s -out %s -p %s',...
                anat_3D,...
                b0_map,...
                mask,...
                fullfile(tmpdir,sprintf('step%02d',i)),...
                parfile{i});
        else
            % after first step, use previous parameter file as initial
            % transform
            elastix_cmd = sprintf('elastix -f %s -m %s -fMask %s -out %s -t0 %s -p %s',...
                anat_3D,...
                b0_map,...
                mask,...
                fullfile(tmpdir,sprintf('step%02d',i)),...
                fullfile(tmpdir,sprintf('step%02d',i-1),'TransformParameters.0.txt'),...
                parfile{i});
        
        end
        % Run registration with Elastix
        system(elastix_cmd)
    end
    
    
    %% APPLY TRANSFORMATION
    
    % First, modify the transform file so that it resamples the DTI image
    % to the dimensions of the original DTI image (and not of the
    % anatomical image to which it defaults)
    transform_file     = fullfile(tmpdir,sprintf('step%02d',i),...
        'TransformParameters.0.txt');
    transform_file_mod = fullfile(tmpdir,'TransformParameters_mod.txt');    
    ModifyTransformFile(transform_file,b0_map,transform_file_mod);
    
    % Now apply this transformation to all stacks in the DTI image, and
    % save again as a new 4D stack   
    for k = 1 : nStacks
        waitbar((k+2) / (nStacks+2),h,sprintf('Transforming stack %d of %d',k,nStacks))
        fname = fullfile(tmpdir,sprintf('image%03d.nii.gz',k));
        
        % Extract a 3D stack
        extract_3Dfrom4D(dti,fname,k);
        
        % Transform the stack
        transformix_cmd = sprintf('transformix -in %s -out %s -tp %s',...
            fname,...
            tmpdir,...
            transform_file_mod);
        system(transformix_cmd);
        
        % Load the result and put into dti_reg_data
        result = load_untouch_nii(fullfile(tmpdir,'result.nii.gz'));
        % Correct the intensity scale of the result image because
        % elastix/transformix set the scl_slope property to 1.
        int_scale = (result.hdr.dime.scl_slope / dti_reg_data.hdr.dime.scl_slope);
        dti_reg_data.img(:,:,:,k) = cast(single(result.img) * int_scale,...
            'like',dti_reg_data.img);
        
    end
    % Save the registered DTI data
    save_untouch_nii(dti_reg_data,dti_reg)
    fprintf('Registered DTI data saved as %s\n',dti_reg)
    
    % Save the elastix transformation file as well.
    movefile(transform_file_mod,strrep(dti_reg,'.nii.gz','_transform.txt'));
    % Close the progress bar.
    close(h)
    
    % Inspect the results in ITK-snap
    if InspectFlag == true
        system(sprintf('itk-snap -g %s -o %s %s &',anat,dti,dti_reg))
    end
    
    % Delete the temporary working directory.
    rmdir(tmpdir,'s')
    
catch ME
    % remove temporary working directory, then throw error message
    rmdir(tmpdir,'s')
    error(ME.message)
end


end

