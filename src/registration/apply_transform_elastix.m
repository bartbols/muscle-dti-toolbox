function apply_transform_elastix( filename_in,filename_out,transform_file,varargin )
%APPLY_TRANSFORM_ELASTIX applies the transform file (output of Elastix) to
%an image. Handles both 3D and 4D data
%
% Bart Bolsterlee, Neuroscience Research Australia
% December 2017
%
% INPUT:
% - filename_in    : filename of NIfTI/STL file to be transformed
% - filename_out   : filename of NIfTI/STL file after transformation
% - transform_file : Elastix transformation file containing the
%                   transformation parameters
% Optional inputs, provided as 'parameter',<value> pairs:
% mask : false/true. If true, the image to be transformed is a binary mask.
%                    The BsplineOrder will be set to 0 to avoid 'folding'
%                    of the transformation.
% ref_image         : filename of NIfTI file of which the spatial
%                    information will be used. If not provided, the spatial
%                    information of 'filename_in' is used.
%
% OUTPUT: none (files will be created)
%
% Example:
%
% To transform image1.nii.gz with transformation file 'transform.txt' and
% save the results as image2.nii.gz:
%
% apply_transform_elastix('image1.nii.gz','image2.nii.gz','transform.txt')
%
p = inputParser;
addRequired(p,'filename_in')
addRequired(p,'filename_out')
addRequired(p,'transform_file')
addParameter(p,'mask',false,@(x) islogical(x) || x==0 || x == 1)
addParameter(p,'ref_image',[])
parse(p,filename_in,filename_out,transform_file,varargin{:})
mask         = p.Results.mask;
ref_image    = p.Results.ref_image;

% Detect what type of data (surface/image) is provided.
if endsWith(filename_in,'stl')
    input_type = 'surface';
elseif endsWith(filename_in,{'.nii','.nii.gz'})
    input_type = 'image';
else
    error('Unknown data type. Input should be a NIfTI-image or an STL-surface.')
end

% Create temporary working directory.
char_list = char(['a':'z' '0':'9']) ;
tmpdir = [];
while exist(tmpdir,'dir') == 7 || isempty(tmpdir)
    tmpdir = fullfile(tempdir,char_list(ceil(length(char_list)*rand(1,8))));
end
mkdir(tmpdir)


try
    switch input_type
        case 'image'
            % Check if image is 2D, 3D or 4D by loading the data
            img_in  = load_untouch_nii(filename_in);
            tf_file = fullfile(tmpdir,'transform.txt');
            
            % Set the spatial domain of the transform file to the domain of the
            % input image or the reference image.
            if isempty(ref_image)
                ModifyTransformFile(transform_file,img_in,tf_file)
            else
                ModifyTransformFile(transform_file,ref_image,tf_file)
            end
%             tf_file = transform_file;
            if img_in.hdr.dime.dim(1) == 3 || img_in.hdr.dime.dim(1) == 2
                % Image is 2D or 3D and can be transformed directly with transformix
                
                if mask == true
                    % If a binary mask is provided as input image, the final
                    % B-spline order should be set to 0 to avoid artefacts in the
                    % transformed image.
                    set_ix(tf_file,'FinalBSplineInterpolationOrder',0)
                end
                transformix_cmd = sprintf('transformix -in %s -out %s -tp %s',...
                    filename_in,tmpdir,tf_file);
                system(transformix_cmd)
                % Create output directory, if it doesn't exist already.
                if exist(fileparts(filename_out),'dir') ~= 7
                    mkdir(fileparts(filename_out))
                end
                
                if ~isempty(ref_image)
                    % Reslice to the reference image dimensions
                    if mask == true
                        ip = 0; % nearest neighbor interpolation
                    else
                        ip = 1; % linear interpolation
                    end
                    I = load_untouch_nii(ref_image);
                    if I.hdr.dime.dim(1) > 3
                        % Image has more than 3 dimensions, which Convert3D
                        % cannot handle. Extract a 3D stack and use that as
                        % the reference image.
                        extract_3Dfrom4D(I,fullfile(tmpdir,'ref_image.nii.gz'),1);
                        ref_image = fullfile(tmpdir,'ref_image.nii.gz');
                    end
                    
                    c3d_cmd = sprintf('c3d -int %d %s %s -reslice-identity -o %s',...
                        ip,ref_image,fullfile(tmpdir,'result.nii.gz'),filename_out);
                    system(c3d_cmd);
                else
                    % Check if result is zipped or not
                    tmpname = ls(fullfile(tmpdir,'result.*'));
                    if endsWith(filename_out,'.nii.gz') && endsWith(tmpname,'.nii')
                        % Result image is currently not zipped, while a zipped file is
                        % requested as output.
                        gzip(fullfile(tmpdir,tmpname))
                        movefile(fullfile(tmpdir,'result.nii.gz'),ref_image)
                    elseif endsWith(filename_out,'.nii') && endsWith(tmpname,'.nii.gz')
                        % Result image is currently zipped, while a unzipped file
                        % is requested.
                        gunzip(fullfile(tmpdir,tmpname))
                        movefile(fullfile(tmpdir,'result.nii'),filename_out)
                    elseif endsWith(filename_out,'.nii.gz') && endsWith(tmpname,'.nii.gz')
                        movefile(fullfile(tmpdir,'result.nii.gz'),filename_out)
                    elseif endsWith(filename_out,'.nii') && endsWith(tmpname,'.nii')
                        movefile(fullfile(tmpdir,'result.nii'),filename_out)
                    end    
                    
%                     movefile(fullfile(tmpdir,'result.nii.gz'),filename_out)
                end
                
                fprintf('Transformed file saved as %s.\n',filename_out)
            elseif img_in.hdr.dime.dim(1) == 4
                % Image is 4D, which transformix cannot handle. Transform each 3D
                % stack independently, then combine again into a new 4D image.
                % Extract all 3D images
                fname = extract_3Dfrom4D(filename_in,tmpdir,0);
                
                if mask == true
                    % If a binary mask is provided as input image, the final
                    % B-spline order should be set to 0 to avoid artefacts in the
                    % transformed image.
                    set_ix(tf_file,'FinalBSplineInterpolationOrder',0)
                end
                
                for i = 1 : img_in.hdr.dime.dim(5)
                    % Transform the stack with transformix
                    transformix_cmd = sprintf('transformix -in %s -out %s -tp %s',...
                        fname{i},tmpdir,tf_file);
                    system(transformix_cmd)
                    
                    % Read the results
                    tf = load_untouch_nii(fullfile(tmpdir,'result.nii.gz'));
                    if i == 1
                        img_out = tf;
                        img_out.hdr.dime.dim(1) = 4;
                        img_out.hdr.dime.dim(5) = img_in.hdr.dime.dim(5);
                    end
                    %  Add stack to 4D image.
                    img_out.img(:,:,:,i) = single(tf.img);
                    
                end
                
                % Create output directory, if it doesn't exist already.
                if exist(fileparts(filename_out),'dir') ~= 7
                    mkdir(fileparts(filename_out))
                end
                save_untouch_nii(img_out,filename_out);
                fprintf('Transformed file saved as %s.\n',filename_out)
            else
                error('Only 3D and 4D data is supported, but the input image is %dD.',img_in.hdr.dime.dim(1))
            end
            
            rmdir(tmpdir,'s')
        case 'surface'
            % Read the surface model
            FV = stlread(filename_in);
            
            % Transform the vertices with transformix
            FV.vertices = transformix_points(FV.vertices, transform_file);
            
            % Write as stl-file.
            stlwrite(filename_out,FV)
            fprintf('In reg_elastix: Transformed surface saved as %s\n',filename_out)
    end
catch ME
    % remove temporary working directory, then throw error message
    rmdir(tmpdir,'s')
    getReport(ME)
    error(ME.message)
end


end

