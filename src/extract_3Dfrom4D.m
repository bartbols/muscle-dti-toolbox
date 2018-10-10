function varargout = extract_3Dfrom4D( filename_4D,filename_3D,i )
%EXTRACT_3DTO4D extracts the i-th 3D stack from the 4D image input. The
%3D image is saved.
%
% INPUT
% filename_4D: filename of 4D nifti image data
% filename_3D: filename of 3D nifti image data or, if i=0, the results
%              folder where all the 3D stacks will be saved to.
% i          : index of the 3D stack in the 4D data. If 0, all stacks will
%              be extracted and stored in the folder 'filename_3D' with the
%              same filename as the 4D data but with _1, _2, _3 etc.
%              appended (e.g. if example.nii.gz contains 3 stacks the
%              following files will be created in the folder filename_3D:
%              example_1.nii.gz,example_2.nii.gz, example_3.nii.gz)
%
% OUTPUT
% nii        : structure array with nifti image information (both data and
%              header)


% Load the 4D data
if ~isstruct(filename_4D)
    % Filename of NIfTI file is provided as input. Load the file.
    nii4D = load_untouch_nii(filename_4D);   % extract image & header
else
    % NIfTI structure is provided as input.
    nii4D = filename_4D;
    filename_4D = [nii4D.fileprefix '.nii.gz'];
end

if nii4D.hdr.dime.dim(1) ~= 4
    error('Input image is not 4D.')
end

if i ~= 0
    % Extract the i-th stack
    nii3D = nii4D;
    nii3D.img = nii4D.img(:,:,:,i);

    % Edit the header to conform with 3D data
    nii3D.hdr.dime.dim(1) = 3;                      % (make it 3D only)
    nii3D.hdr.dime.dim(5:8) = 1;                    % (make it 3D only)

    % Save the 3D data
    save_untouch_nii(nii3D,filename_3D)
    fprintf('Stack %d saved as %s\n',i,filename_3D)
    filename_out = filename_3D;
elseif i == 0
    filename_out = cell(1,nii4D.hdr.dime.dim(5));
    for i = 1 : nii4D.hdr.dime.dim(5)
        % Extract the i-th stack
        nii3D(i)     = nii4D;
        nii3D(i).img = nii4D.img(:,:,:,i);

        % Edit the header to conform with 3D data
        nii3D(i).hdr.dime.dim(1) = 3;                      % (make it 3D only)
        nii3D(i).hdr.dime.dim(5:8) = 1;                    % (make it 3D only)

        % Save the 3D data
        [~,b,c] = fileparts(filename_4D);
        fname = strrep([b c],'.nii.gz',['_' int2str(i) '.nii.gz']);
        fullfilename = fullfile(filename_3D,fname);
        if exist(filename_3D,'dir')~=7
            mkdir(filename_3D)
            fprintf('Folder ''%s'' created.\n',filename_3D)
        end
        
        save_untouch_nii(nii3D(i),fullfilename)
        fprintf('Stack %d saved as %s\n',i,fullfilename)
        filename_out{i} = fullfilename;

    end
end

if nargout > 0
    varargout{1} = filename_out;
    if nargout > 1
        varargout{2} = nii3D;
    end
end
end

