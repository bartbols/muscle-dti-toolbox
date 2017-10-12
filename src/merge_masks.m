function multi_label = merge_masks(multi_label_filename,varargin)
%MERGE_MASKS merges a series of masks containing only one label to one
%multi-label mask.
%
% Bart Bolsterlee, Neuroscience Research Australia (NeuRA)
% October 2017
%
% INPUT:
% - multi_label_filename: filename of the multi-label image
%
% All subsequent input arguments should come in pairs with the first
% argument the filename of the single-label image and the second argument
% the label number it should get in the multi-label image.
% 
% For example, to combine two images with only one label to a multi-label
% image containing two labels:
% 
% merged_filename =
%     merge_mask('multi_label_img.nii.gz','label1.nii.gz',1,'label2.nii.gz',2)
%
% All non-zero values in label 1 will be assigned number 1 in the
% multi-label images. All non-zero values in label 2 will be assigned
% number 2 in the multi-label images, etc.
%
% And for three labels:
% merged_filename =
%     merge_mask('multi_label_img.nii.gz','label1.nii.gz',1,'label2.nii.gz',2,...
%                'label3.nii.gz',3)
%
% OUTPUT:
% multi_label : nifti-structure with the multi-label image as created by
% this function.
%


if ~endsWith(multi_label_filename,'.nii.gz')
    error('output filename should end with .nii.gz')
end

inp = varargin;
n = length(inp);

% Load all label files
label_files = inp(1:2:end);
label_numbers = inp(2:2:end);

for i = 1 : n/2
    % Load the label files
    if ~endsWith(label_files{i},'.nii.gz')
        error('Label files should end with .nii.gz')
    end
    current_label = load_untouch_nii(label_files{i});
    
    if i == 1
        % Make empty label with same header information as the first header
        multi_label = current_label;
        multi_label.img = cast(zeros(size(current_label.img)),'like',current_label.img);
    end
    
    % Add the current label to the multi-label image with the assigned
    % label number
    multi_label.img(current_label.img ~= 0) = cast(label_numbers{i},'like',multi_label.img);
    
end

% Save the multi-label image
save_untouch_nii(multi_label,multi_label_filename);
fprintf('Multi-label image saved as %s\n',multi_label_filename)
for i = 1 : n/2
    fprintf('Label %d: %s\n',label_numbers{i},char(label_files{i}))
    
end
end

