function MD_mask_filename = MakeMDmask( fib_filename,md_threshold,varargin )
%MAKEMDMASK Reads the mean diffusivity (MD) data from the DSI-reconstructed
% .fib.gz-file, creates a binary mask between the md_threshold values and 
% saves the mask as a NIFTI file that can be used as a mask in DSI studio.
%
% Bart Bolsterlee, Neuroscience Research Australia (NeuRA)
% August 2017
%
% ----------------- USAGE ----------------- 
% MD_mask_filename = MakeMDmask( fib,md_threshold )
% or
% MD_mask_filename = MakeMDmask( fib,md_threshold,...
%                                'ResultsPath',<results_path>,...
%                                'filename',<md_mask_filename>)
%
% ----------------- INPUT -----------------
% - fib_filename : filename of .fib.gz file as reconstructed with DSI studio.
% - md_threshold : 1x2 vector with lower and upper limit for fractional 
%                  anisotropy values. Values outside of this range are made
%                  0.
% Optional inputs, provided as 'ParameterName',<value> inputs:
%
% - ResultsPath:  folder where the MD mask will be saved to. If not
%                 provided, the folder in which the fib-file lives is used.
% - filename    : filename of the MD mask, without the path. If not provided, the
%                 fib_filename is used with _MD.nii.gz appended.
%                              
%
% ----------------- OUTPUT -----------------
% - MD_mask_filename (optional): filename of the MD mask.

%% Check inputs
p = inputParser;
addRequired(p,'fib_filename',@(x) ~isempty(strfind(x,'.fib.gz')))
addRequired(p,'md_threshold',@(x) validateattributes(x,{'numeric'},{'size',[1 2]}))
addParameter(p,'filename',[],@(x) ~isempty(strfind(x,'.nii.gz')))
addParameter(p,'ResultsPath',fileparts(fib_filename),@ischar)
parse(p, fib_filename,md_threshold,varargin{:})
MD_mask_filename = p.Results.filename;
ResultsPath      = p.Results.ResultsPath;

% Create results folder if it doesnt' exist
if exist(p.Results.ResultsPath,'dir') ~= 7
    mkdir(p.Results.ResultsPath)
    fprintf('Results directory created: %s\n',p.Results.ResultsPath)
end

%% If filename is not provided use the name of the .fib.gz file 
if isempty(MD_mask_filename)
    [pathstr,name,ext] = fileparts(fib_filename) ;
    MD_mask_filename = [name(1:end-4) '_MDmask.nii.gz'];
end
full_MD_mask_name = fullfile(ResultsPath,MD_mask_filename);

%%
% Call DSI Studio to export an MD map from the fib-file
commandTxt2 = sprintf('dsi_studio --action=exp --source=%s --export=md',...
    fib_filename);
system(commandTxt2);

% Create a binary mask with voxels in between the MD thresholds 1 and all
% other voxels 0. Save as a nifti file.
md_map = load_untouch_nii([fib_filename '.md.nii.gz']);
delete([fib_filename '.md.nii.gz'])
md_mask = md_map;
md_mask.img = cast((md_map.img < md_threshold(1) | md_map.img > md_threshold(2)),'like',md_map.img);
clear md_map

% Save the MD mask as a .nii.gz file
save_untouch_nii(md_mask,full_MD_mask_name)

fprintf('MD mask created with filename %s\n',full_MD_mask_name)
end % of function

