function nii = extract_3Dfrom4D( filename_4D,filename_3D,i )
%EXTRACT_3DTO4D extracts the i-th 3D stack from the 4D image input. The
%3D image is saved.
%
% INPUT
% filename_4D: filename of 4D nifti image data
% filename_3D: filename of 3D nifti image data
% i          : index of the 3D stack in the 4D data
%
% OUTPUT      
% nii        : structure array with nifti image information (both data and
%              header)


% Load the 4D data
    nii = load_untouch_nii(filename_4D);   % extract image & header
    
% Extract the i-th stack
    nii.img = nii.img(:,:,:,i);
    
% Edit the header to conform with 3D data
    nii.hdr.dime.dim(1) = 3;                      % (make it 3D only)
    nii.hdr.dime.dim(5:8) = 1;                    % (make it 3D only)

% Save the 3D data

save_untouch_nii(nii,filename_3D)
fprintf('Stack %d saved as %s\n',i,filename_3D)





end

