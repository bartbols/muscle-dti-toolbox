function predict_mask_per_slice( fixed,moving,mask_ref,mask_pred,parfile,varargin)
%PREDICT_MASK_PER_SLICE Registers the mask from the fixed to the moving
%image on a slice by slice basis using Elastix. The fixed and moving images
%are expected to have the same dimensions and expected to be aligned per
%slice.
%
% Bart Bolsterlee, Neuroscience Research Australia
% December 2017
%
% INPUT:
% - fixed          : filename of the fixed image
% - moving         : filename of the moving image
% - mask_ref       : filename of the reference mask
% - mask_pred      : filename of the mask that will be predicted
% - parfile        : Elastix parameter file containing the
%                    transformation parameters (2D registration).
%
% Optional inputs, provided as 'parameter',<value> pairs:
% stack_f          : if fixed image is 4D, the channel number used for
%                    registration.
% stack_m          : if moving image is 4D, the channel number used for
%                    registration.
% slice            : vector of slice numbers used for registration. If empty,
%                    all slices will be predicted.
% label_number     : Label number in the mask that will be predicted.
%                    Only required if mask contains multiple labels.
%                    Default is to use all labels in the mask.
%
% OUTPUT: none (mask-image <mask_pref> will be created)
p = inputParser;
addRequired(p,'fixed')
addRequired(p,'moving')
addRequired(p,'mask_ref')
addRequired(p,'mask_pred')
addRequired(p,'parfile')
addParameter(p,'stack_f',[],@(x) assert(isscalar(x)))
addParameter(p,'stack_m',[],@(x) assert(isscalar(x)))
addParameter(p,'slice',[],@(x) assert(isnumeric(x)));
addParameter(p,'label_number',[],@(x) assert(isnumeric(x)))
parse(p, fixed,moving,mask_ref,mask_pred,parfile,varargin{:})

stack_f = p.Results.stack_f;
stack_m = p.Results.stack_m;
slice   = p.Results.slice;
label_number = p.Results.label_number;
%
% Create temporary working directory.
char_list = char(['a':'z' '0':'9']) ;
tmpdir = [];
while exist(tmpdir,'dir') == 7 || isempty(tmpdir)
    tmpdir = fullfile(pwd,char_list(ceil(length(char_list)*rand(1,8))));
end
mkdir(tmpdir)

try
    % Load the mask and the image data
    fprintf('Loading mask image...')
    mask   = load_untouch_nii(mask_ref);
    fprintf(' completed.\n')
    
    fprintf('Loading fixed image...')
    fixed  = load_untouch_nii(fixed);
    fprintf(' completed.\n')
    
    fprintf('Loading moving image...')
    moving = load_untouch_nii(moving);
    fprintf(' completed.\n')

    % Extract slice numbers in which the mask is present and save as 2D
    % images.
    fixed_2D = fixed;
    fixed_2D.hdr.dime.dim(1) = 2;
    fixed_2D.hdr.dime.dim(4:5) = 1;
    fixed_2D.hdr.dime.pixdim(4) = 0;
    
    moving_2D = moving;
    moving_2D.hdr.dime.dim(1) = 2;
    moving_2D.hdr.dime.dim(4:5) = 1;
    moving_2D.hdr.dime.pixdim(4) = 0;
    
    mask_2D = mask;
    mask_2D.hdr.dime.dim(1) = 2;
    mask_2D.hdr.dime.dim(4:5) = 1;
    mask_2D.hdr.dime.pixdim(4) = 0;
        
    if isempty(slice)
        slice_numbers = 1 : mask.hdr.dime.dim(4);
    else
        slice_numbers = slice;
    end
    fprintf('Extracting 2D images...')
    slice_numbers2 = slice_numbers;
    for slice_nr = slice_numbers
        if all(all(mask.img(:,:,slice_nr)==0))
            slice_numbers2(slice_numbers2==slice_nr)=[];
            continue
        end
        
        if exist(fullfile(tmpdir,'fixed'),'dir') ~= 7;mkdir(fullfile(tmpdir,'fixed'));end
        if exist(fullfile(tmpdir,'moving'),'dir') ~= 7;mkdir(fullfile(tmpdir,'moving'));end
        if exist(fullfile(tmpdir,'mask'),'dir') ~= 7;mkdir(fullfile(tmpdir,'mask'));end
        
        % Extract 2D slice from the fixed image
        if isempty(stack_f)
            fixed_2D.img = fixed.img(:,:,slice_nr);
        else
            fixed_2D.img = fixed.img(:,:,slice_nr,stack_f);
        end
        save_untouch_nii(fixed_2D,fullfile(tmpdir,'fixed',sprintf('slice%03d.nii.gz',slice_nr)))
        
        % Extract 2D slice from the moving image
        if isempty(stack_f)
            moving_2D.img = moving.img(:,:,slice_nr);
        else
            moving_2D.img = moving.img(:,:,slice_nr,stack_m);
        end
        save_untouch_nii(moving_2D,fullfile(tmpdir,'moving',sprintf('slice%03d.nii.gz',slice_nr)))
        
        % Extract 2D slice from the mask
        mask_2D.img = mask.img(:,:,slice_nr);
        save_untouch_nii(mask_2D,fullfile(tmpdir,'mask',sprintf('slice%03d.nii.gz',slice_nr)))
    end
    fprintf(' completed.\n')
    
    % Make empty 3D mask with spatial information of moving image
    mask_tf = mask;
    mask_tf.hdr.hist = mask.hdr.hist;
    mask_tf.img = cast(zeros(size(mask.img)),'like',mask.img);
    
    hwait = waitbar(0,'',...
        'Name','Progress bar predict_mask_per_slice');
    tic
    counter = 0;
    for slice_nr = slice_numbers2
        counter = counter + 1;
        waitbar((counter-1) / length(slice_numbers2),hwait,sprintf('Slice %d of %d',counter,length(slice_numbers2)))
        try
            filename.fixed_2D  = fullfile(tmpdir,'fixed',sprintf('slice%03d.nii.gz',slice_nr));
            filename.moving_2D = fullfile(tmpdir,'moving',sprintf('slice%03d.nii.gz',slice_nr));
            filename.mask_2D   = fullfile(tmpdir,'mask',sprintf('slice%03d.nii.gz',slice_nr));
            
%             if exist(filename.fixed_2D,'file') ~= 2;continue;end
            
            if exist(fullfile(tmpdir,'results'),'dir') == 7
                rmdir(fullfile(tmpdir,'results'),'s')
            end
            mkdir(fullfile(tmpdir,'results'))
            
            reg_elastix(filename.fixed_2D,filename.moving_2D,parfile,...
                'mask',filename.mask_2D,...
                'dilate_mask',2,...
                'transform_inv',fullfile(tmpdir,'results','transform_inv.txt'),...
                'label_number',label_number)
            
            apply_transform_elastix(filename.mask_2D,...
                fullfile(tmpdir,'results','mask_2D_tf.nii.gz'),...
                fullfile(tmpdir,'results','transform_inv.txt'),...
                'ref_image',filename.moving_2D,...
                'mask',true)
            
            mask_2D_tf = load_untouch_nii(fullfile(tmpdir,'results','mask_2D_tf.nii.gz'));
            
            if isempty(label_number)
                mask_tf.img(:,:,slice_nr) = cast(mask_2D_tf.img,'like',mask.img);
            else
                mask_tf.img(:,:,slice_nr) = cast((mask_2D_tf.img == label_number)*label_number,'like',mask.img);
            end
        catch ME
            error(ME.message)
            continue
        end
        % Create the output folder, if it doesn't exist already.
        if exist(fileparts(mask_pred),'dir') ~= 7
            mkdir(fileparts(mask_pred))
        end
        save_untouch_nii(mask_tf,mask_pred);        
    end
    close(hwait)
    fprintf('Predicted mask saved as %s.\n',mask_pred)
    t_elapsed = toc;
    fprintf('Registration took %.1f seconds.\n',t_elapsed)
    rmdir(tmpdir,'s')
catch ME
    % remove temporary working directory, then throw error message
    rmdir(tmpdir,'s')
    error(ME.message)
end




end

