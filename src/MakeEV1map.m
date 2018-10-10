function varargout = MakeEV1map( fib_filename,DTI_filename, fa_threshold,varargin )
%FIB_TO_NII Reads the primary eigenvector data from the DSI-reconstructed
% .fib-file and saves it as a nifti-file with the same metadata as the
% original DTI data.
%
% Bart Bolsterlee, Neuroscience Research Australia (NeuRA)
% February 2017
%
% ----------------- USAGE ----------------- 
% EV1_filename = MakeEV1map( fib_filename,DTI_filename, fa_threshold)
% or 
% EV1_filename = MakeEV1map( fib_filename,DTI_filename, fa_threshold,EV1_filename)
%
% ----------------- INPUT ----------------- 
% - fib_filename     : filename of .fib file as reconstructed with DSI studio
% - DTI_filename     : filename of .nii file with DTI data
% - fa_threshold     : 1x2 vector with lower and upper limit for fractional 
%                      anisotropy values. Values outside of this range are 
%                      made 0 (black in EV1 image).
%
% Optional 4th input argument:
% - EV1_filename : filename of new file with primary eigenvector data
%
% ----------------- OUTPUT ----------------- 
%-  EV1_filename      : filename of new file with primary eigenvector data
%-  EV1               : NIfTI structure with primary eigenvector data

%% Check inputs
p = inputParser;
addRequired(p,'fib_filename',@(x) contains(x,'.fib'))
addRequired(p,'DTI_filename',@(x) contains(x,'.nii'))
addRequired(p,'fa_threshold',@(x) validateattributes(x,{'numeric'},{'size',[1 2]}))
addOptional(p,'EV1_filename',[DTI_filename(1:end-7) '_EV1.nii.gz'],@(x) contains(x,'.nii'))
parse(p,fib_filename,DTI_filename,fa_threshold,varargin{:});
EV1_filename = p.Results.EV1_filename;
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
perm_dim = [2 3 4 1];
% perm_dim = [2 3 4 1];
eigvec_map = permute(reshape(fib_data.dir0,[3 fib_data.dimension]),perm_dim);

% Filter the eigenvector map based on FA
R = eigvec_map(:,:,:,1);
G = eigvec_map(:,:,:,2);
B = eigvec_map(:,:,:,3);
R(fa_map < fa_threshold(1) | fa_map > fa_threshold(2)) = 0;
G(fa_map < fa_threshold(1) | fa_map > fa_threshold(2)) = 0;
B(fa_map < fa_threshold(1) | fa_map > fa_threshold(2)) = 0;
% eigvec_map(:,:,:,1) = abs(flip(R,2)) * 100;
% eigvec_map(:,:,:,2) = abs(flip(G,2)) * 100;
% eigvec_map(:,:,:,3) = abs(flip(B,2)) * 100;

% Make all eigenvectors point in the positive z-direction.
R(B<0) = -R(B<0);
G(B<0) = -G(B<0);
B(B<0) = -B(B<0);

scaling = 1;
eigvec_map(:,:,:,1) = R * scaling;
eigvec_map(:,:,:,2) = -G * scaling; % negative to correct for DSI studio's coordinate system
eigvec_map(:,:,:,3) = B * scaling;

% Use all information from the DTI nifti file but overwrite the image data
% with the eigenvector map.
% The second dimension needs to be flipped to have a good match in ITK-snap
% between the anatomical data and the eigenvector map. I'm not sure if this
% also means that the y-direction (or maybe x-direction?) needs to be
% flipped. So it could be that the directions are not stored correctly. For
% segmentation in ITK-snap this doesn't really matter because there will be
% good contrast in colour between tissues with different eigenvectors
% anyway.
% Note 05/02/2018: The eigenvectors are now saved in the ITK coordinate
% system, which is flipped in the x and y coordinates compared to the NIfTI
% coordinate system.

if DTI_nii.hdr.hist.srow_x(1) > 0
    eigvec_map = flip(eigvec_map,1);
    eigvec_map(:,:,:,1) = -eigvec_map(:,:,:,1);
end
if DTI_nii.hdr.hist.srow_y(2) > 0
    eigvec_map = flip(eigvec_map,2);
    eigvec_map(:,:,:,2) = -eigvec_map(:,:,:,2);
end
    
DTI_nii.img = single(eigvec_map);
DTI_nii.hdr.dime.dim(5) = 3;
DTI_nii.hdr.dime.scl_slope = 1; %/scaling;
DTI_nii.hdr.dime.scl_inter = 0;
DTI_nii.hdr.dime.bitpix = 32;
DTI_nii.hdr.dime.datatype = 16;
save_untouch_nii(DTI_nii,EV1_filename)

if nargout > 0
    varargout{1} = EV1_filename;
    if nargout > 1
        varargout{2} = DTI_nii;
    end
end


end % of function

