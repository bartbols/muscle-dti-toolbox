function [T,bvec_corr] = correct_bvec(DTI_filename,old_bvec_filename,new_bvec_filename)
% This file corrects the bvec-file, which contains the gradient directions,
% for oblique image acquisition. A new bvec-file (txt-file) is created with
% the corrected gradient directions.
%
% Bart Bolsterlee, Neuroscience Research Australia (NeuRA)
% February 2017
%
% ----------------- USAGE ----------------- 
% [T,bvec_corr] = correct_bvec(DTI_filename,old_bvec_filename,new_bvec_filename)
%
% ----------------- INPUT -----------------
% DTI_filename      : the filename of the nifti file containing the DTI data
% old_bvec_filename : the filename of the bvec-file as obtained from the
%                     scanner
% new_bvec_filename : the filename of the new, corrected bvec-file
%
% ----------------- OUTPUT -----------------
% T                 : the transformation matrix stored in the header of the 
%                     nifti-file
% bvec_corr         : corrected bvec table as written to the corrected
%                     bvec-file
%
% This script is based on the following email conversation with Frank Yeh,
% developer of DSI studio on February 1st, 2016.
% 
%     You can correct the bvec using the following steps:
% 
% 1 First get the sform matrix (I assume this is the oblique matrix).
% 
% -1.875 0.000379747 -0.00985536
% 9.66389e-05 1.86953 0.38152
% -0.00371396 -0.143069 4.98541
% 
% 2. Correct for voxel size by multiple it with the following matrix
% 1/1.875 0 0
% 0 1/1.875 0
% 0 0 1/5
% 
% 3. multiply bvec with the corrected sform matrix
% 4. scale bvec, making them unique vectors.
% 
%    This should fix the oblique problem. However, a 4 degrees of the
% rotation angle is too small to be detectable in the final result. I
% would recommend you run a test scan using a large oblique angle to see
% if the correction is correct.

%% Check inputs
if nargin ~= 3
    error('Wrong number of input arguments. Usage is correct_bvec(DTI_filename,old_bvec_filename,new_bvec_filename)')
end
p = inputParser;
addRequired(p,'DTI_filename',@(x) ~isempty(strfind(x,'.nii.gz')))
addRequired(p,'old_bvec_filename',@(x) ~isempty(strfind(x,'.bvec')))
addRequired(p,'new_bvec_filename',@(x) ~isempty(strfind(x,'.bvec')))
parse(p,DTI_filename,old_bvec_filename,new_bvec_filename)


%%

% Load the data
bvec = load(old_bvec_filename,'-ascii');
img = load_untouch_nii(DTI_filename);

% Get transformation from nifti-header
T = [img.hdr.hist.srow_x;img.hdr.hist.srow_y;img.hdr.hist.srow_z];

% Get rotation matrix and voxel size
R = T(1:3,1:3);
voxelsize = sqrt(sum(R .^2));
S = diag(1./voxelsize);

% Correct the transformation matrix for the voxel size
Rnew = R*S;

% Two options: pre or post-multiply. Which one is correct?
% Option 1: pre-multiply
% !!! According to Frank Yeh (DSI Studio developer) this option is correct !!!
bvec_corr = Rnew * bvec; % (they are already unit vectors so they don't need to be normalized)

% % Option 2: post-multiply
% bvec_corr2 = (bvec' * Rnew)';

% Flip the directions according to the sign of the transformation matrix
% NOTE: I am not 100% certain this is correct, so always check fibre tracts 
% for plausibility of the orientation.

bvec_corr = bvec_corr .* repmat(sign(diag(Rnew)),1,size(bvec,2));

save(new_bvec_filename,'bvec_corr','-ascii')
fprintf('Corrected bvec-file saved as %s\n',new_bvec_filename)
end