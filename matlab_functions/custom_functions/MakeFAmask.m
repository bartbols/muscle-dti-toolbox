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
% FA_mask_filename = MakeFAmask( fib,fa_threshold, FA_mask_filename)
%
% ----------------- INPUT -----------------
% - fib_filename : filename of .fib.gz file as reconstructed with DSI studio.
% - fa_threshold : 1x2 vector with lower and upper limit for fractional 
%                  anisotropy values. Values outside of this range are made
%                  0.
%
% Optional 3rd input argument
% - FA_mask_filename: filename of the FA mask. If not provided, the
%                    fib_filename is used with _FA.nii.gz appended.
%                              
%
% ----------------- OUTPUT -----------------
% - FA_mask_filename (optional): filename of the FA mask.

%% Check inputs
p = inputParser;
addRequired(p,'fib_filename',@(x) ~isempty(strfind(x,'.fib.gz')))
addRequired(p,'fa_threshold',@(x) validateattributes(x,{'numeric'},{'size',[1 2]}))
addOptional(p,'FA_mask_filename',[],@(x) ~isempty(strfind(x,'.nii.gz')))
parse(p, fib_filename,fa_threshold,varargin{:})
FA_mask_filename = p.Results.FA_mask_filename;

%% If filename is not provided, create the map in the subfolder 'masks' of 
% the folder in which the fibre file lives.
if isempty(FA_mask_filename)
    [pathname,fname,ext] = fileparts(fib_filename);
    FA_mask_filename = fullfile(pathname,'masks',[fname(1:end-4) '_FA.nii.gz']);
    if exist(fullfile(pathname,'masks'),'dir') ~= 7
        mkdir(pathname,'masks');
        fprintf('Created folder %s\n',fullfile(pathname,'masks'))
    end
end
%%
% Check extension
if ~strcmp(fib_filename(end-6:end) ,'.fib.gz')
    error('Extension of fib-file needs to be ''.fib.gz''')
end
% Check dimensions of fa_threshold (has to be 1x2 array);
validateattributes(fa_threshold,{'numeric'},{'ncols',2,'nrows',1})


%%
% Call DSI Studio to export an FA map from the fib-file
commandTxt2 = sprintf('dsi_studio --action=exp --source=%s --export=fa0',...
    fib_filename);
system(commandTxt2);

% Create a binary mask with voxels in between the FA thresholds 1 and all
% other voxels 0. Save as a nifti file.
fa_map = load_untouch_nii([fib_filename '.fa0.nii.gz']);
fa_mask = fa_map;
fa_mask.img = int16((fa_map.img < fa_threshold(1) | fa_map.img > fa_threshold(2)));
delete([fib_filename '.fa0.nii.gz'])
clear fa_map

% Save the FA mask as a .nii.gz file
save_untouch_nii(fa_mask,FA_mask_filename)

fprintf('FA mask created with filename %s\n',FA_mask_filename)
end % of function

