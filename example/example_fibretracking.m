% Clear the workspace
clear

%% Set some filenames
% This assumes that the data have been preprocessed with 'example_preprocessing.m' 
% and that masks were created with 'example_make_masks.m'.

results_path   = fullfile(pwd,'results','001');

% Filenames for tractography
% Name of the the fibre file created with DSI studio (the .fib.gz file that
% is the result of Preprocessing_and_DTI_recon)
TrackFileNames.FIB     = fullfile(results_path,'CALF001_DTI_LPCA.fib.gz');
TrackFileNames.DTI     = fullfile(results_path,'CALF001_DTI_LPCA.nii.gz');

% Filename to which the tracts will be saved. Always include the file
% extension .mat.
TrackFileNames.Tracts  = fullfile(results_path,'tracts','CALF001_DTI_MG.mat');

% Set the seed, boundary (TER=terminative region) and ROA (region of avoidance)
% masks for fibre tractography.
TrackFileNames.Seed       = fullfile(results_path,'DTI_masks','CALF001_DTI_MG_seed.nii.gz');
TrackFileNames.TER        = fullfile(results_path,'DTI_masks','CALF001_DTI_MG_boundary.nii.gz');
TrackFileNames.ROA        = fullfile(results_path,'DTI_masks','CALF001_DTI_LPCA_FA.nii.gz');

% Load the surface model
SurfFilename     = fullfile(results_path,'DTI_masks','CALF001_DTI_MG.stl');
SurfModel        = stlread(SurfFilename);

%% Fibre track settings
% type 'help TrackFibres' for the meaning of each of parameters
TrackSettings.MinLength    = 20;
TrackSettings.MaxLength    = 200;
TrackSettings.FA_threshold = 0.10;
TrackSettings.Stepsize     = 1;
TrackSettings.FiberCount   = 250;
TrackSettings.SeedCount    = [];
TrackSettings.MaxAngle     = 10;
TrackSettings.Smoothing    = 0;
TrackSettings.MaxTime      = 30;


%% Perform fibre tracking, truncation and extrapolation
% Type help TrackFibres for more information on how to use this function.
[DTItracts,StopFlag] = TrackFibres(TrackFileNames,TrackSettings);

% Reconstruct facicles and measure architecture
% Type help CalcArchitecture for more information on how to use this function.
DTItracts = CalcArchitecture(DTItracts,SurfModel);

% Add the order of the polynomial fit as an optional 3rd input argument if
% you don't want to use the default order of 3. For example, to fit a
% second-order polynomial:
% DTItracts = CalcArchitecture(DTItracts,surf_model, 2);

% Calculate DTI measures (fractional anisotropy, mean diffusivity and
% eigenvalues) for each tract. This will add the fields fa, md, lambda1,
% lambda2 and lambda3 to DTI tracts.
DTItracts = CalcDTI_indices(DTItracts,TrackFileNames.FIB);

% Save the fields in DTItracts as individual variables using the -struct option
% for saving.
save(TrackFileNames.Tracts,'-struct','DTItracts')

% You can load the data into the workspace again as follows:
% DTItracts = load(TrackFileNames.Tracts);

%% Plot the results
% The function InspectTracts can be used to visually inspect the results of
% fibre tracking and subsequent analyses. It will plot the fibre tracts in
% 3D (by default, the raw and the polynomial fitted tracts will be
% displayed) and the distribution of some architectural parameters and DTI
% indices. You can also call InspectTracts without any input arguments; the
% user is then requested to select a tract file from a file selection
% dialog box.
%
% Type 'help InspectTracts' to explore the options for usage.

InspectTracts('Tracts',DTItracts,...
    'SurfModel',SurfModel);

% Alternatively, you can give the filenames of the DTItracts and surface 
% model (instead of the structure arrays)

% handles = InspectTracts('Tracts',TrackFileNames.Tracts,...
%     'SurfModel',SurfFilename);
