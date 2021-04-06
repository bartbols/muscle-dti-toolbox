function varargout = MakeLAMBDA_map( fib_filename,DTI_filename, fa_threshold,varargin )
%MAKELAMBDA_MAP Reads the eigenvalue data from the DSI-reconstructed
% .fib-file and saves it as a nifti-file with the same metadata as the
% original DTI data.
%
% Bart Bolsterlee, Neuroscience Research Australia (NeuRA)
% February 2017
%
% ----------------- USAGE ----------------- 
% LAMBDA_filename = MakeLAMBDA_map( fib_filename,DTI_filename, fa_threshold)
% or 
% LAMBDA_filename = MakeLAMBDA_map( fib_filename,DTI_filename, fa_threshold,LAMBDA_filename)
%
% ----------------- INPUT ----------------- 
% - fib_filename     : filename of .fib file as reconstructed with DSI studio
% - DTI_filename     : filename of .nii file with DTI data
% - fa_threshold     : 1x2 vector with lower and upper limit for fractional 
%                      anisotropy values. Values outside of this range are 
%                      made 0 (black in LAMBDA image).
%
% Optional 4th input argument:
% - LAMBDA_filename : filename of new file with eigenvalue data
%
% ----------------- OUTPUT ----------------- 
%-  LAMBDA_filename   : filename of new file with eigenvalue data
%-  LAMBDA            : NIfTI structure with eigenvalue data

%% Check inputs
p = inputParser;
addRequired(p,'fib_filename',@(x) contains(x,'.fib'))
addRequired(p,'DTI_filename',@(x) contains(x,'.nii'))
addRequired(p,'fa_threshold',@(x) validateattributes(x,{'numeric'},{'size',[1 2]}))
addOptional(p,'LAMBDA_filename',[DTI_filename(1:end-7) '_LAMBDA.nii.gz'],@(x) contains(x,'.nii'))
parse(p,fib_filename,DTI_filename,fa_threshold,varargin{:});
LAMBDA_filename = p.Results.LAMBDA_filename;
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

fa_map = reshape(fib_data.fa0,fib_data.dimension);

% Load the nifti file with DTI data
DTI_nii = load_untouch_nii(DTI_filename);

% Create new nifti-file with same metadata as original DTI data but with
% the primary eigenvector data (3-channels) as image data
% perm_dim = [2 3 4 1];
% LAMBDA_map = permute(reshape(fib_data.dir0,[3 fib_data.dimension]),perm_dim);

% Filter the eigenvector map based on FA
% L1 = LAMBDA_map(:,:,:,1);
% L2 = LAMBDA_map(:,:,:,2);
% L3 = LAMBDA_map(:,:,:,3);
L1 = reshape(fib_data.ad,fib_data.dimension);
L2 = reshape(fib_data.rd1,fib_data.dimension);
L3 = reshape(fib_data.rd2,fib_data.dimension);
L1(fa_map < fa_threshold(1) | fa_map > fa_threshold(2)) = 0;
L2(fa_map < fa_threshold(1) | fa_map > fa_threshold(2)) = 0;
L3(fa_map < fa_threshold(1) | fa_map > fa_threshold(2)) = 0;

scaling = 1;
LAMBDA_map(:,:,:,1) = L1 * scaling;
LAMBDA_map(:,:,:,2) = L2 * scaling;
LAMBDA_map(:,:,:,3) = L3 * scaling;

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

if DTI_nii.hdr.hist.srow_x(1) > 0
    LAMBDA_map = flip(LAMBDA_map,1);
%     LAMBDA_map(:,:,:,1) = LAMBDA_map(:,:,:,1);
end
if DTI_nii.hdr.hist.srow_y(2) > 0
    LAMBDA_map = flip(LAMBDA_map,2);
%     LAMBDA_map(:,:,:,2) = LAMBDA_map(:,:,:,2);
end
    
DTI_nii.img = single(LAMBDA_map);
DTI_nii.hdr.dime.dim(5) = 3;
DTI_nii.hdr.dime.scl_slope = 1; %/scaling;
DTI_nii.hdr.dime.scl_inter = 0;
DTI_nii.hdr.dime.bitpix = 32;
DTI_nii.hdr.dime.datatype = 16;
save_untouch_nii(DTI_nii,LAMBDA_filename)

if nargout > 0
    varargout{1} = LAMBDA_filename;
    if nargout > 1
        varargout{2} = DTI_nii;
    end
end


end % of function

