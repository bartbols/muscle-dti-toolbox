clear
addpath(genpath('../matlab_functions'))

%% Set filenames for the example data
% Always use full path names (not relative to current directory). You can
% use MATLAB's pwd command to get a string with the current directory, from
% which you can build up your full filename using fullfile. Type help
% fullfile for more info.
% The extension for image data has to be .nii.gz and should be included in
% the filename.
datapath       = fullfile(pwd,'data','001');   % raw data path for subject 1
results_path   = fullfile(pwd,'results','001'); % results path for subject 1
DTI_filename   = fullfile(datapath,'CALF001_DTI.nii.gz'); 
bval_filename  = fullfile(datapath,'CALF001_DTI.bval');
bvec_filename  = fullfile(datapath,'CALF001_DTI.bvec');

%% Preprocess the DTI data and create the src/fib file
% By calling Preprocessing_and_DTI_recon, the b-vectors are corrected, the
% raw DTI data is filtered and DSI Studio is called to reconstruct the
% diffusion tensor. This will create the src.gz and fib.gz files that are
% required for fibre tracking.
% Type help Preprocessing_and_DTI_recon for more information on how to use
% this function.
FilterFlag = true; 
filename = Preprocessing_and_DTI_recon('DTI',DTI_filename,...
    'bval',bval_filename,'bvec',bvec_filename,...
    'filter',FilterFlag,...
    'ResultsPath',results_path);

%% Create masks
% The next lines of code will create masks at the resolution of the DTI scan that 
% can be used for tractography. The tractography masks will be created 
% from the muscle segmentation of the anatomical scan, which will be
% downsampled to match the DTI resolution.
label_filename = fullfile(datapath,'masks','CALF001_T1_anat_labels.nii.gz');

% Make masks/surfaces for all labels present in the label file
% Type help MakeSurfaceAndMasks for more information on how to use this
% function.
mask_filenames = MakeSurfaceAndMasks( label_filename,DTI_filename,...
    'ResultsPath',fullfile(results_path,'DTI_masks'));

% Make a mask to exclude regions with FA values below/above the provided
% threshold.
% Type help MakeFAmask for more information on how to use this function.
FA_threshold = [0.1 0.5];
filename.FA_mask = MakeFAmask(filename.FIB,FA_threshold,...
    fullfile(results_path,'DTI_masks','CALF001_DTI_LPCA_FA.nii.gz'));

