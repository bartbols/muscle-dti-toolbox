function FA_mask_filename = MakeFAmask( fib_filename,fa_threshold,varargin )
%MAKEFAMASK Reads the FA data from the DSI-reconstructed .fib.gz-file, creates
% a binary mask between the fa_threshold values and saves the mask as a 
% .nii.gz file that can be used as a mask in DSI studio.
%
% Bart Bolsterlee, Neuroscience Research Australia (NeuRA)
% February 2017
%
% ----------------- USAGE ----------------- 
% FA_mask_filename = MakeFAmask( fib,fa_threshold )
% or
% FA_mask_filename = MakeFAmask( fib,fa_threshold,...
%                                'ResultsPath',<results_path>,...
%                                'filename',<fa_mask_filename>)
%
% ----------------- INPUT -----------------
% - fib_filename : filename of .fib.gz file as reconstructed with DSI studio.
% - fa_threshold : 1x2 vector with lower and upper limit for fractional 
%                  anisotropy values. Values outside of this range are made
%                  0.
% Optional inputs, provided as 'ParameterName',<value> inputs:
%
% - ResultsPath:  folder where the FA mask will be saved to. If not
%                 provided, the folder in which the fib-file lives is used.
% - filename    : filename of the FA mask, without the path. If not provided, the
%                 fib_filename is used with _FA.nii.gz appended.
%                              
%
% ----------------- OUTPUT -----------------
% - FA_mask_filename (optional): filename of the FA mask.

%% Check inputs
p = inputParser;
addRequired(p,'fib_filename',@(x) ~isempty(strfind(x,'.fib.gz')))
addRequired(p,'fa_threshold',@(x) validateattributes(x,{'numeric'},{'size',[1 2]}))
addParameter(p,'filename',[],@(x) ~isempty(strfind(x,'.nii.gz')))
addParameter(p,'ResultsPath',fileparts(fib_filename),@ischar)
parse(p, fib_filename,fa_threshold,varargin{:})
FA_mask_filename = p.Results.filename;
ResultsPath      = p.Results.ResultsPath;

% Create results folder if it doesnt' exist
if exist(p.Results.ResultsPath,'dir') ~= 7
    mkdir(p.Results.ResultsPath)
    fprintf('Results directory created: %s\n',p.Results.ResultsPath)
end

%% If filename is not provided use the name of the .fib.gz file 
if isempty(FA_mask_filename)
    [pathstr,name,ext] = fileparts(fib_filename) ;
    FA_mask_filename = [name(1:end-4) '_FAmask.nii.gz'];
end
full_FA_mask_name = fullfile(ResultsPath,FA_mask_filename);

%%
% Call DSI Studio to export an FA map from the fib-file
commandTxt2 = sprintf('dsi_studio --action=exp --source=%s --export=fa0',...
    fib_filename);
system(commandTxt2);

% Create a binary mask with voxels in between the FA thresholds 1 and all
% other voxels 0. Save as a nifti file.
fa_map = load_untouch_nii([fib_filename '.fa0.nii.gz']);
delete([fib_filename '.fa0.nii.gz'])
fa_mask = fa_map;
fa_mask.img = cast((fa_map.img < fa_threshold(1) | fa_map.img > fa_threshold(2)),'like',fa_map.img);
clear fa_map

% Save the FA mask as a .nii.gz file
save_untouch_nii(fa_mask,full_FA_mask_name)

fprintf('FA mask created with filename %s\n',full_FA_mask_name)
end % of function

