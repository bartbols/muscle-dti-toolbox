function varargout = MakeFA_map( fib_filename,DTI_filename, fa_threshold,varargin )
%MAKEFA_MAP Reads the fractional anisotropy data from the DSI-reconstructed
% .fib-file and saves it as a nifti-file with the same metadata as the
% original DTI data.
%
% Bart Bolsterlee, Neuroscience Research Australia (NeuRA)
% February 2017
%
% ----------------- USAGE ----------------- 
% FA_filename = MakeFA_map( fib_filename,DTI_filename, fa_threshold)
% or 
% FA_filename = MakeFA_map( fib_filename,DTI_filename, fa_threshold,FA_filename)
%
% ----------------- INPUT ----------------- 
% - fib_filename     : filename of .fib file as reconstructed with DSI studio
% - DTI_filename     : filename of .nii file with DTI data
% - fa_threshold     : 1x2 vector with lower and upper limit for fractional 
%                      anisotropy values. Values outside of this range are 
%                      made 0 (black in FA image).
%
% Optional 4th input argument:
% - FA_filename : filename of new file with eigenvalue data
%
% ----------------- OUTPUT ----------------- 
%-  FA_filename   : filename of new file with eigenvalue data
%-  FA            : NIfTI structure with eigenvalue data

%% Check inputs
p = inputParser;
addRequired(p,'fib_filename',@(x) contains(x,'.fib'))
addRequired(p,'DTI_filename',@(x) contains(x,'.nii'))
addRequired(p,'fa_threshold',@(x) validateattributes(x,{'numeric'},{'size',[1 2]}))
addOptional(p,'FA_filename',[DTI_filename(1:end-7) '_FA.nii.gz'],@(x) contains(x,'.nii'))
parse(p,fib_filename,DTI_filename,fa_threshold,varargin{:});
FA_filename = p.Results.FA_filename;
%%

% Load the fib file
fprintf('Loading fib-file...')
[path,name,ext] = fileparts(fib_filename);
if strcmp(ext,'.gz')
    gunzip(fib_filename)
    fib_data = load(fib_filename(1:end-3),'-mat');
    delete(fib_filename(1:end-3))
elseif strcmp(ext,'.fib')
    fib_data = load(fib_filename,'-mat');
else
error('Unknown file format %s. The fib-file should have extension .fib or .fib.gz',ext)
end
fprintf(' completed.\n')


% Load the nifti file with DTI data
DTI_nii = load_untouch_nii(DTI_filename);

% Create new nifti-file with same metadata as original DTI data but with
% the primary eigenvector data (3-channels) as image data
% perm_dim = [2 3 4 1];
% FA_map = permute(reshape(fib_data.dir0,[3 fib_data.dimension]),perm_dim);

% Exclude voxels outside the FA range
fa_map = reshape(fib_data.fa0,fib_data.dimension);
fa_map(fa_map < fa_threshold(1) | fa_map > fa_threshold(2)) = 0;
% L1 = reshape(fib_data.ad,fib_data.dimension);
% L2 = reshape(fib_data.rd1,fib_data.dimension);
% L3 = reshape(fib_data.rd2,fib_data.dimension);
% L1(fa_map < fa_threshold(1) | fa_map > fa_threshold(2)) = 0;
% L2(fa_map < fa_threshold(1) | fa_map > fa_threshold(2)) = 0;
% L3(fa_map < fa_threshold(1) | fa_map > fa_threshold(2)) = 0;

scaling = 1;
FA_map(:,:,:,1) = fa_map * scaling;
% FA_map(:,:,:,2) = L2 * scaling;
% FA_map(:,:,:,3) = L3 * scaling;

% Use all information from the DTI nifti file but overwrite the image data
% with the eigenvalue maps.
% The second dimension needs to be flipped to have a good match in ITK-snap
% between the anatomical data and the eigenvalue map. I'm not sure if this
% also means that the y-direction (or maybe x-direction?) needs to be
% flipped. So it could be that the directions are not stored correctly. For
% segmentation in ITK-snap this doesn't really matter because there will be
% good contrast in colour between tissues with different eigenvectors
% anyway.
% Note 05/02/2018: The eigenvectors are now saved in the ITK coordinate
% system, which is flipped in the x and y coordinates compared to the NIfTI
% coordinate system.

FA_nii = DTI_nii;
if FA_nii.hdr.hist.srow_x(1) > 0
    FA_map = flip(FA_map,1);
%     FA_map(:,:,:,1) = FA_map(:,:,:,1);
end
if FA_nii.hdr.hist.srow_y(2) > 0
    FA_map = flip(FA_map,2);
%     FA_map(:,:,:,2) = FA_map(:,:,:,2);
end
  
FA_nii.img = single(FA_map);
FA_nii.hdr.dime.dim(5) = 1;
FA_nii.hdr.dime.scl_slope = 1; %/scaling;
FA_nii.hdr.dime.scl_inter = 0;
FA_nii.hdr.dime.bitpix = 32;
FA_nii.hdr.dime.datatype = 16;
save_untouch_nii(FA_nii,FA_filename)

if nargout > 0
    varargout{1} = FA_filename;
    if nargout > 1
        varargout{2} = FA_nii;
    end
end


end % of function

