clear
addpath(genpath('../bin'))
addpath(genpath('../src'))

%% Set filenames for the example data
% Always use full path names (not relative to current directory). You can
% use MATLAB's pwd command to get a string with the current directory, from
% which you can build up your full filename using fullfile. Type help
% fullfile for more info.
% The extension for image data has to be .nii.gz and should be included in
% the filename.
datapath       = fullfile(pwd,'data','001');   % raw data path for subject 1
results_path   = fullfile(pwd,'results','001'); % results path for subject 1
filename.DTI   = fullfile(datapath,'CALF001_DTI.nii.gz'); 
filename.bval  = fullfile(datapath,'CALF001_DTI.bval');
filename.bvec  = fullfile(datapath,'CALF001_DTI.bvec');

%% Preprocess the DTI data and create the src/fib file
% By calling Preprocessing_and_DTI_recon, the b-vectors are corrected, the
% raw DTI data is filtered and DSI Studio is called to reconstruct the
% diffusion tensor. This will create the src.gz and fib.gz files that are
% required for fibre tracking.
% Type help Preprocessing_and_DTI_recon for more information on how to use
% this function.
FilterFlag = true; 
filename = Preprocessing_and_DTI_recon('DTI',filename.DTI,...
    'bval',filename.bval,'bvec',filename.bvec,...
    'filter',FilterFlag,...
    'ResultsPath',results_path);
