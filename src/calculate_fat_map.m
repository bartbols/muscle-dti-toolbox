function varargout = calculate_fat_map( mDixon,fat_map,varargin )
%CALCULATE_FAT_MAP Calculates a fat fraction map form the mDixon scan and
%saves the the map as a NIfTI file.
%
% ----------------- USAGE -----------------
% nii = calculate_fat_map(mDixon,fat_map,varargin)
%
% ----------------- INPUT -----------------
% ----- REQUIRED -----
% mDixon    : NIfTI filename or structure containing the mDixon data.
% fat_map   : filename of the fat_map (this file will be created).
% 
% ----- OPTIONAL -----
% Optional inputs are provided as pairs of 'ParameterName',<value>
% - 'water' : stack number of the water image. Default: 4
% - 'fat'   : stack number of the fat image.   Default: 1
%
% ----------------- OUTPUT -----------------
% Only one (optional) output can be provided, which will be the NIfTI
% structure containing the fat map.

p = inputParser;
addRequired(p,'mDixon')
addRequired(p,'fat_map',@(x) contains(x,'.nii'))
addParameter(p,'water',4,@isscalar)
addParameter(p,'fat',1,@isscalar)
parse(p,mDixon,fat_map,varargin{:});

water = p.Results.water;
fat   = p.Results.fat;
    
if ~isstruct(mDixon)
    % filename is provided
    % check if filename of the mDixon scan is not the same as the fat map
    if strcmp(mDixon,fat_map)
        error('Filename of mDixon file and fat map cannot be the same.')
    end
    mDixon = load_untouch_nii(mDixon);
end


F = mDixon;
F.hdr.dime.dim(1) = 3;
F.hdr.dime.dim(5) = 1;
F.img = single(mDixon.img(:,:,:,fat)) ./ (single(mDixon.img(:,:,:,fat)) + single(mDixon.img(:,:,:,water)));
F.hdr.dime.datatype = 16;
F.hdr.dime.bitpix = 32;
F.hdr.dime.scl_slope = 1;
F.hdr.dime.scl_inter = 0;
F.hdr.dime.glmin = 0;
F.hdr.dime.glmax = 1;

save_untouch_nii(F,fat_map);
fprintf('Fat map saved as %s.\n',fat_map)

if nargout == 1
    varargout{1} = F;
end







end

