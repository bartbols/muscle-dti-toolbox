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
% parfile:  fliename of the elastix parameter file with registration parameters
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
addRequired(p,'parfile')
addRequired(p,'dti_reg',@(x) contains(x,'.nii.gz'))
addParameter(p,'mask',[],@(x) isempty(x) || contains(x,'.nii.gz'))
addParameter(p,'foreground_threshold',10,@(x) assert(isscalar(x)))
addParameter(p,'stack',[],@(x) assert(isscalar(x)))
addParameter(p,'b0_stack',1,@(x) assert(isscalar(x)))
parse(p,dti,anat,parfile,dti_reg,varargin{:});

% Assign parameter values
foreground_threshold = p.Results.foreground_threshold;
mask                 = p.Results.mask;
stack                = p.Results.stack;
b0_stack             = p.Results.b0_stack;

% Create temporary working directory.
char_list = char(['a':'z' '0':'9']) ;
tmpdir = [];
while exist(tmpdir,'dir') == 7 || isempty(tmpdir)
    tmpdir = fullfile(pwd,char_list(ceil(length(char_list)*rand(1,8))));
end
mkdir(tmpdir)

try
    % Set some filenames
    b0_map         = fullfile(tmpdir,'b0_map.nii.gz');
    anat_resampled = fullfile(tmpdir,'anat_resampled.nii.gz');
    mask_resampled = fullfile(tmpdir,'mask_resampled.nii.gz');
    
    if ~isempty(stack)
        % Extract one stack from the anatomical scan if the data is 4D
        anat_tmp = fullfile(tmpdir,'anat_tmp.nii.gz');
        extract_3Dfrom4D(anat,anat_tmp,stack);
    else
        anat_tmp = anat;
    end
    
    % Create mask if foreground_threshold is provided 
    
    if ~isempty(foreground_threshold)
        if ~isempty(mask)
            error('mask and foreground_threhsold cannot both be defined. Remove either mask or foreground_threshold from the inputs.')
        else
            mask = fullfile(tmpdir,'mask_tmp.nii.gz');
        end
        anat_data = load_untouch_nii(anat_tmp);
        anat_data.img = cast(anat_data.img>foreground_threshold,'like',anat_data.img);
        anat_data.hdr.dime.scl_slope = 1;
        save_untouch_nii(anat_data,mask);
    end
    
    % Load the DTI data
    dti_data = load_untouch_nii(dti);
    
    % Create an empty nifti-structure that will be filled with the
    % registered data.
    dti_reg_data = dti_data;
    dti_reg_data.img = cast(zeros(size(dti_data.img)),'like',dti_data.img);
    nStacks = dti_data.hdr.dime.dim(5);

    h = waitbar(0,'Registering DTI data...',...
        'Name','Progress bar DTI registration');
    waitbar(1 / (nStacks+2),h,'Resampling anatomical data and mask...')
    
    % Step 1: Extract B0-map from 4D dti data. Elastix does not take in 4D
    % data, so a 3D image should be created before registration can be done.
    % Assume that B0 map is the first image in the 4D DTI data.
    extract_3Dfrom4D(dti,b0_map,b0_stack);
    
    % Resample anatomical data and mask to the resolution of the DTI scan.
    new_dim = dti_data.hdr.dime.dim(2:4);
    c3d_cmd = sprintf('c3d %s -interpolation 1 -resample %dx%dx%dvox -o %s',...
        anat_tmp,...
        new_dim(1),new_dim(2),new_dim(3),....
        anat_resampled);
    fprintf('Resampling anatomical scan...')
    system(c3d_cmd);
    fprintf(' completed.\n')
    
    % Resample mask    
    c3d_cmd = sprintf('c3d %s -interpolation 0 -resample %dx%dx%d -o %s',...
        mask,...
        new_dim(1),new_dim(2),new_dim(3),....
        mask_resampled);
    fprintf('Resampling mask ...')
    system(c3d_cmd);
    fprintf(' completed.\n')
    
    % Make sure that the downsampled anatomical scan and mask have the same 
    % header info as the DTI data. Even though they have the same number of
    % voxels and bounding box, the header is sometimes different, which
    % causes (subtle) misinterpretation of the registration results.
    
    mask_resampled_data = load_untouch_nii(mask_resampled);
    anat_resampled_data = load_untouch_nii(anat_resampled);
    mask_resampled_data.hdr.hist = dti_data.hdr.hist;
    anat_resampled_data.hdr.hist = dti_data.hdr.hist;
    save_untouch_nii(mask_resampled_data,mask_resampled);
    save_untouch_nii(anat_resampled_data,anat_resampled);
    
    waitbar(2 / (nStacks+2),h,...
        'Registering the B0-map to the anatomical scan.')
    
    % Register B0 map to the anatomical scan
    elastix_cmd = sprintf('elastix -f %s -m %s -fMask %s -p %s -out %s',...
        anat_resampled,...
        b0_map,...
        mask_resampled,...
        parfile,...
        tmpdir);
    system(elastix_cmd)
    
    % Now apply this transformation to all stacks in the DTI image, and
    % save again as a new 4D stack
    transform_file = fullfile(tmpdir,'TransformParameters.0.txt');
    for k = 1 : nStacks
        waitbar((k+2) / (nStacks+2),h,sprintf('Transforming stack %d of %d',k,nStacks))
        fname = fullfile(tmpdir,sprintf('image%03d.nii.gz',k));
        % Extract a stack
        extract_3Dfrom4D(dti,fname,k);
        % Transform the stack
        transformix_cmd = sprintf('transformix -in %s -out %s -tp %s',...
            fname,...
            tmpdir,...
            transform_file);
        system(transformix_cmd);
        
        % Downsample to the original DTI dimensions again.
        c3d_cmd = sprintf('c3d %s -interpolation 1 -resample %dx%dx%dvox -o %s',...
            fullfile(tmpdir,'result.nii.gz'),...
            dti_data.hdr.dime.dim(2),...
            dti_data.hdr.dime.dim(3),...
            dti_data.hdr.dime.dim(4),....
            fullfile(tmpdir,'result_dti_res.nii.gz'));
        fprintf('Resampling registered DTI stack...')
        system(c3d_cmd);
        fprintf(' completed.\n')
        
        % Load the result and put into dti_reg_data
        result = load_untouch_nii(fullfile(tmpdir,'result_dti_res.nii.gz'));
        dti_reg_data.img(:,:,:,k) = result.img;
    
    end
    % Save the registered DTI data
    save_untouch_nii(dti_reg_data,dti_reg)
    fprintf('Registered DTI data saved as %s\n',dti_reg)
    % Delete the temporary working directory
    rmdir(tmpdir,'s')
    close(h)
catch ME
    % remove temporary working directory, then throw error message
    rmdir(tmpdir,'s')
    error(ME.message)
end


end

