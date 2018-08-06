function extract_2Dfrom3D( filename_in,output_folder,varargin )
%EXTRACT_2DTO3D extracts 2D images from a 3D/4D image input. The slices are
%saved as separate NIfTI files in the designated output folder.
%
% INPUT
% filename_in: filename of 3D/4D nifti image data
% output_folder: folder in which all the 2D slices will be saved.
% Optional inputs, provided as <parameter>,'value' pairs:
% slice    : list of slices to extract. Default: extract all slices
% stack    : stack number to extract 2D images from in case the input file
%            is 4D. Default: 1
%
p = inputParser;
addRequired(p,'filename_in',@(x) contains(x,'.nii'))
addRequired(p,'output_folder')
addParameter(p,'slice',[],@isnumeric)
addParameter(p,'stack',[],@isscalar)
parse(p,filename_in,output_folder,varargin{:});
stack = p.Results.stack;
slice = p.Results.slice;


if exist(output_folder,'dir') ~= 7
    mkdir(output_folder);
end

% Load the image data
if ~isstruct(filename_in)
    % Filename of NIfTI file is provided as input. Load the file.
    nii = load_untouch_nii(filename_in);   % extract image & header
else
    % NIfTI structure is provided as input.
    nii = filename_in;
    filename_in = [nii.fileprefix '.nii.gz'];
end

% Check dimensionality of the data

switch nii.hdr.dime.dim(1)
    case 3
        image3D = nii.img;
    case 4
        if isempty(stack)
            error('A 4D image is provided, but no stack number to extract the slices from is provide. Add the stack number as an input argument')
        end
        image3D = nii.img(:,:,:,stack);
            
    otherwise
        
    error('Input image should be 3D or 4D.')    
end

if isempty(slice)
    % Extract all slices
    slice = 1 : size(image3D,3);
end

for slice_nr = slice
    nii2D = nii;
    % Extract the i-th stack
    nii2D.img = image3D(:,:,slice_nr);

    % Edit the header to conform with 3D data
    nii2D.hdr.dime.dim(1) = 2;                      % (make it 2D only)
    nii2D.hdr.dime.dim(4:8) = 1;                    % (make it 2D only)

    % Save the 2D data
    filename_2D = fullfile(output_folder,sprintf('slice%03d.nii.gz',slice_nr));
    save_untouch_nii(nii2D,filename_2D)
    fprintf('Slice %d saved as %s\n',slice_nr,filename_2D)
end
end % of function

