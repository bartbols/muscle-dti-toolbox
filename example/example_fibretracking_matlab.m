% This script is an example of fibre tracking using the MATLAB-based
% tracking algorithms dti_tractography.m
% Clear the workspace
clear

%% Set some filenames
% This assumes that the data have been preprocessed with 'example_preprocessing.m' 
% and that masks were created with 'example_make_masks.m'.

results_path   = fullfile(pwd,'results','001');

% Filenames for tractography
% Name of the the fibre file created with DSI studio (the .fib.gz file that
% is the result of Preprocessing_and_DTI_recon)
filename.FIB     = fullfile(results_path,'CALF001_DTI_LPCA.fib.gz');
filename.DTI     = fullfile(results_path,'CALF001_DTI_LPCA.nii.gz');

% Filename to which the tracts will be saved. Always include the file
% extension .mat.
filename.Tracts  = fullfile(results_path,'tracts','CALF001_DTI_MG_matlab.mat');

% Set the seed, boundary (TER=terminative region) and ROA (region of avoidance)
% masks for fibre tractography.
% filename.Seed       = fullfile(results_path,'DTI_masks','CALF001_DTI_MG_seed.nii.gz');
filename.TER        = fullfile(results_path,'DTI_masks','CALF001_DTI_MG_boundary.nii.gz');
filename.FULL       = fullfile(results_path,'DTI_masks','CALF001_DTI_MG_full.nii.gz');

% Load the surface model
SurfFilename     = fullfile(results_path,'DTI_masks','CALF001_DTI_MG.stl');
SurfModel        = stlread(SurfFilename);

%% Fibre track settings
% type 'help dti_tractography' for the meaning the parameters
TrackSettings.MinLength    = 20;
TrackSettings.MaxLength    = 200;
TrackSettings.FA_threshold = [0.1 0.5]; % minimum and maximum FA
TrackSettings.StepSize     = 1;
TrackSettings.MaxAngle     = 10;

%% Make seeds
% The MATLAB algorithms requires the user to provide a list of seed
% locations (in global coordinates). The function MakeSeeds can be used to
% make a regular grid of seed points with a given spacing within a seed
% mask.
spacing_mm = 2.5; % spacing of the grid in mm
seeds = MakeSeeds(filename.FULL,filename.DTI,spacing_mm);

% %% Plot the seeds and the surface model
% figure
% h1 = patch('Vertices',SurfModel.vertices,...
%     'Faces',SurfModel.faces,...
%     'FaceColor','y',...
%     'FaceAlpha',0.2,...
%     'EdgeColor','none',...
%     'EdgeAlpha',0.5);
% axis equal off
% hold on
% h2 = plot3(seeds(1,:),...
%      seeds(2,:),...
%      seeds(3,:),...
%      'ro',...
%      'MarkerFaceColor','r',...
%      'MarkerSize',2);
%  
%  legend([h1 h2],{'muscle model','seeds'})
%  view(35,50)
 
%% Perform fibre tracking, truncation and extrapolation
% Type 'help dti_tractography' for more information on how to use this function.
DTItracts = dti_tractography(filename.FIB,seeds,TrackSettings,...
    'DTI',filename.DTI,...
    'TER',filename.TER);

% Reconstruct facicles and measure architecture
% Type 'help CalcArchitecture' for more information on how to use this function.
DTItracts = CalcArchitecture(DTItracts,SurfModel);

% Add the order of the polynomial fit as an optional 3rd input argument if
% you don't want to use the default order of 3. For example, to fit a
% second-order polynomial:
% DTItracts = CalcArchitecture(DTItracts,surf_model, 2);

% Calculate DTI measures (fractional anisotropy, mean diffusivity and
% eigenvalues) for each tract. This will add the fields fa, md, lambda1,
% lambda2 and lambda3 to DTI tracts.
DTItracts = CalcDTI_indices(DTItracts,filename.FIB);

% Save the fields in DTItracts as individual variables using the -struct option
% for saving.
save(filename.Tracts,'-struct','DTItracts')

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
