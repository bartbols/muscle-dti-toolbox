
function filename = MakeSurfaceAndMasks(segm_filename, DTI_filename, varargin )
%MAKESURFACEANDMASKS Creates a surface mesh (.stl file) and masks (.nii.gz
% or .nii files) at the resolution of the DTI scan for tractography.
%
% Bart Bolsterlee, Neuroscience Research Australia (NeuRA)
% February 2017
%
% Note: this function  requires Convert3D to be installed on the pc and
% added to the path so that 'c3d' is recognised as an external command.
% Convert3D can be downloaded here:
% http://www.itksnap.org/pmwiki/pmwiki.php?n=Downloads.C3D
%
% ----------------- USAGE -----------------
% filename = MakeSurfaceAndMasks(segm_filename,DTI_filename,varargin)
%
% ----------------- INPUT -----------------
% ----- REQUIRED -----
% - segm_filename: filename of the image file containing the labels at the
% resolution of the anatomical scan. Must include extension .nii or .nii.gz. In the
% same folder as the label file, there must be a .lab file with label
% information. The label info file should contain three columns with the
% label number (the number in segmentation data) in the first column, the
% short name of the muscle in the second and the full name of the muscle in
% the third. Example of a .lab-file:
% 3, MG, medial gastrocnemius
% 5, LG, lateral gastrocnemius
% This means that label 3 contains the medial gastrocnemius (short name:
% MG) and label 5 contains the lateral gastrocnemius (short name: LG)
%
% - DTI_filename: filename of the image file containing the DTI data to which
% the labels will be resampled to. Must include extension .nii.gz.
%
% ----- OPTIONAL -----
% Optional inputs are provided as pairs of 'ParameterName',<value> (e.g.
% 'LabelNumbers',[1 3 5])
%
% - LabelNumbers: list of label numbers in the segmentation file to create
% mask files for.
% - LabelNames: cell string with list of label names to create mask files for.
% - MaskPrefix: prefix for the mask and surface model filename. The short
% name of the label (e.g. MG) will be added to the prefix. The default is
% the DTI filename.
% - ResultsPath: folder in which the mask and surface file will be saved. The
%   default is the folder in which the segmentation file lives.
% - NumberOfVoxelsToRemove: number of layers of voxels to remove from the
%   full mask to generate the seed mask. Default = 2
% - MakeMasks  : if true, the tractography masks are created. Default = true
% - MakeSurface: if true, surfaces (.stl) are created. Default = true.
% - FillHoles  : if true, all holes in the binary mask will be filled prior
%                to surface reconstruction. Default = true.
%
%  MakeMasks and MakeSurface can be set to false if only surfaces or
%  only masks are required. If both set to false, the function does nothing.
%
% Note:
% If LabelNumbers and LabelNames are not provided, masks and surfaces
% will be created for all labels in the .lab file.
%
% If both LabelNumbers and LabelNames are provided, the masks and surfaces
% are created for the labels in LabelNumbers.

% ----------------- OUTPUT -----------------
% - filename: a n x 1 struct with fields boundary_mask, seed_mask, full_mask
% and surface containing the filenames of the mask files and surface model
% files that were created. n is the number of labels for which files were
% created.
%
% Example 1 (using LabelNames):
% filename = MakeSurfaceAndMasks('segmentation.nii.gz', 'DTI_data.nii.gz', ...
% 'LabelNames',{'MG,'LG'})
%
% Example 2 (using LabelNumbers):
% filename = MakeSurfaceAndMasks('segmentation.nii.gz', 'DTI_data.nii.gz', ...
% 'LabelNumbers',[2 4],'ResultsPath','MyMasks')


%% Check input arguments
p = inputParser;
addRequired(p,'segm_filename',@(x) contains(x,'.nii'))
addRequired(p,'DTI_filename',@(x) contains(x,'.nii'))
addParameter(p,'LabelNumbers',[],@isnumeric)
addParameter(p,'LabelNames',[],@(x) assert(ischar(x) || iscell(x)))
addParameter(p,'MaskPrefix',[],@ischar)
addParameter(p,'ResultsPath',fileparts(segm_filename),@ischar)
addParameter(p,'NumberOfVoxelsToRemove',2,@(x)(isnumeric(x) && mod(x,1)==0 && x>0))
addParameter(p,'MakeSurface',true,@(x) x==0 || x==1 || islogical(x) )
addParameter(p,'MakeMasks',true,@(x) x==0 || x==1 || islogical(x) )
addParameter(p,'FillHoles',true,@(x) x==0 || x==1 || islogical(x) )
parse(p,segm_filename,DTI_filename,varargin{:});

MakeSurface = p.Results.MakeSurface; % if true, surfaces are made
MakeMasks   = p.Results.MakeMasks; % if true, masks are made
FillHoles   = p.Results.FillHoles;

%% Create results folder if it doesnt' exist
if exist(p.Results.ResultsPath,'dir') ~= 7
    mkdir(p.Results.ResultsPath)
    fprintf('Results directory created: %s\n',p.Results.ResultsPath)
end

% Create temporary working directory.
char_list = char(['a':'z' '0':'9']) ;
tmpdir = [];
while exist(tmpdir,'dir') == 7 || isempty(tmpdir)
    tmpdir = fullfile(tempdir,char_list(ceil(length(char_list)*rand(1,8))));
end
mkdir(tmpdir)


try
    %%
    % Read in the label information file
    if endsWith(segm_filename,'.gz')
        LabelInfoFilename = [segm_filename(1:end-7) '.lab'];
    elseif endsWith(segm_filename,'.nii')
        LabelInfoFilename = [segm_filename(1:end-4) '.lab'];
    else
        error('Unknown file type for segmentation file. The file should have extension .nii or .nii.gz.')
    end
    assert(exist(LabelInfoFilename,'file') == 2,...
        'Label information file %s does not exist. Create a valid label information file and store in the same folder as the segmentation file.',LabelInfoFilename)
    
    fid = fopen(LabelInfoFilename,'r');
    LabelInfo = textscan(fid,'%d %s %s','Delimiter',',');
    fclose(fid);
    
    % Select which labels to process
    if isempty(p.Results.LabelNumbers) && isempty(p.Results.LabelNames)
        fprintf('Label numbers or names not defined. Process all labels.\n')
        LabelNumToProcess = cell2mat(LabelInfo(1));
        idxToProcess = (1 : length(LabelNumToProcess))';
    elseif ~isempty(p.Results.LabelNumbers)
        LabelNumToProcess = p.Results.LabelNumbers;
        if ~isempty(p.Results.LabelNames)
            warning('Both LabelNames and LabelNumbers are defined. The labels given by LabelNumbers will be processed.')
        end
        
        % Check if the label numbers that are provided for processing exist in
        % the label file
        idxNotPresent = find(~ismember(LabelNumToProcess,LabelInfo{1}));
        if ~isempty(idxNotPresent)
            fprintf('Label not present in label info file\n')
            fprintf('The following labels are provided in %s:\n',LabelInfoFilename)
            for i = 1 : size(LabelInfo{2},1)
                fprintf('%-3d %-8s %s\n',LabelInfo{1}(i),LabelInfo{2}{i},LabelInfo{3}{i})
            end
            error('Label %d not found in label info file <a href="matlab: open(''%s'')">%s</a>.',...
                LabelNumToProcess(idxNotPresent),LabelInfoFilename,LabelInfoFilename)
        end
        [~,idxToProcess] = ismember(LabelInfo{1},LabelNumToProcess);
        idxToProcess = find(idxToProcess);
    elseif ~isempty(p.Results.LabelNames)
        % Check if the label names that are provided for processing exist in
        % the label file
        LabelNames = p.Results.LabelNames; % names provided for processing
        idxNotPresent = find(~ismember(LabelNames,LabelInfo{2}));
        if ~isempty(idxNotPresent)
            fprintf('Label not present in label info file\n')
            fprintf('The following labels are provided in %s:\n',LabelInfoFilename)
            for i = 1 : size(LabelInfo{2},1)
                fprintf('%-3d %-8s %s\n',LabelInfo{1}(i),LabelInfo{2}{i},LabelInfo{3}{i})
            end
            error('Label %s not found in label info file <a href="matlab: open(''%s'')">%s</a>.',...
                LabelNames{idxNotPresent(1)},LabelInfoFilename,LabelInfoFilename)
            
        end
        [~,idx,~] = intersect(LabelInfo{2},LabelNames);
        idxToProcess = LabelInfo{1}(idx);
        LabelNumToProcess = LabelInfo{1}(idxToProcess);
    end
    
    if isempty(p.Results.MaskPrefix)
        [~,tmp] = fileparts(DTI_filename);
        MaskPrefix = tmp(1:end-4);
    else
        MaskPrefix = p.Results.MaskPrefix;
    end
    
    % Check if the label numbers are actually present in the segmentation data
    % file
    label_img = load_untouch_nii(segm_filename);
    ExistingLabelNumbers = double(unique(label_img.img(:)));
    idx = find(~ismember(LabelNumToProcess,ExistingLabelNumbers));
    if ~isempty(idx)
        fprintf('Label %d (%s) does not exist in %s\n',...
            LabelNumToProcess(idx(1)),LabelInfo{2}{idx(1)},segm_filename)
        fprintf('The following label numbers are available: ')
        fprintf('%d ',ExistingLabelNumbers)
        fprintf('\n')
        error('Label does not exist. Modify the label info file <a href="matlab: open(''%s'')">%s</a>',...
            LabelInfoFilename,LabelInfoFilename)
        
    end
    %%
    % Now all labels have been confirmed to exist, create mask files and
    % a surface file for each selected label.
    
    % Extract a 3D stack from the 4D DTI image to be used by convert3D
    % (which can't handle 4D date) to resample the mask from anatomical
    % image dimensions to DTI image dimensions. This only needs to be done
    % once for all labels, so that's why it's placed outside the loop.
    extract_3Dfrom4D(DTI_filename,fullfile(tmpdir,'DTI_3D.nii.gz'),1);
    c = 0;
    for CurrentIdx = idxToProcess'
        c = c + 1;
        CurrentLabelName = LabelInfo{2}{CurrentIdx};
        CurrentLabelNr   = LabelInfo{1}(CurrentIdx);
        
        fprintf('Processing label %d (%s)...\n\n',CurrentLabelNr,CurrentLabelName);
        
        % ------------- MASK RESAMPLING --------
        % Load the DTI data. The masks will be resampled to the resolution of the
        % DTI scan using c3d.
        DTI_data = load_untouch_nii(DTI_filename);
        siz = DTI_data.hdr.dime.dim(2:4);
        
        % Create a new, temporary file with only the selected label.
        mask     = label_img;
        mask.img = cast(label_img.img == CurrentLabelNr,'like',mask.img);
        mask.hdr.dime.scl_slope = 1;
        mask.hdr.dime.scl_inter = 0;
        save_untouch_nii(mask,fullfile(tmpdir,'mask.nii.gz')); % has dimensions of anatomical scan
        
        if MakeSurface == true
            % Make a file with isotropic dimensions for surface model generation.
            commandTxt = sprintf('c3d -int 3 %s -resample-mm 1.5x1.5x1.5mm -o %s',...
                fullfile(tmpdir,'mask.nii.gz'),fullfile(tmpdir,'mask_iso.nii.gz'));
            [~,cmdout] = system(commandTxt);
            % Check if Convert 3D is installed.
            if contains(cmdout,'not recognized as an internal or external command')
                error('Convert3D not found. Please install Convert3D and add to the path environment.')
            end
            % Read into the workspace
            mask_iso = load_untouch_nii(fullfile(tmpdir,'mask_iso.nii.gz'));
            
            % -------------------- SURFACE MODELS ----------------------
            % Create a surface model from the binary mask using the MATLAB toolbox
            % iso2mesh
            filename(c).surface = fullfile(p.Results.ResultsPath, sprintf('%s_%s.stl',MaskPrefix,CurrentLabelName));
            opt.radbound = 1;
            opt.maxsurf = 1;
            method = 'cgalsurf';
            
            binary_mask = mask_iso.img;
            
            if FillHoles == true
                for slice_nr = 1 : size(binary_mask,3)
                    binary_mask(:,:,slice_nr) = imfill(binary_mask(:,:,slice_nr),'holes');
                end
            end
            
            [FV.vertices,FV.faces]   = v2s(binary_mask,0.5,opt,method);
            [FV.vertices,FV.faces]   = surfreorient(FV.vertices,FV.faces);
            
            % Apply some smoothing
            [conn,connnum,count] = meshconn(FV.faces,size(FV.vertices,1));
            %             FV.vertices = smoothsurf(FV.vertices,[],conn,10,0.1,'laplacian');
            FV.vertices = smoothsurf(FV.vertices,[],conn,50,0.7,'lowpass');
            
            % Transform to global coordinates
            % An offset of 1/2 voxel should be applied to have the surface
            % model in the correct coordinate system after transformation
            % to global coordinates (because voxel coordinates present the
            % centre of the voxel, not the edge like v2s assumes).
            FV.vertices = FV.vertices - 0.5;
            
            T = [mask_iso.hdr.hist.srow_x;...
                mask_iso.hdr.hist.srow_y;...
                mask_iso.hdr.hist.srow_z;...
                0 0 0 1];
            if all(all(T(1:3,1:3)==0))
                % The srow information is missing from the header.
                % Calculate the spatial transformation matrix from the
                % quaternion parameters.
                T = makeT_from_quat( mask_iso );
            end
            
            tf = T * [FV.vertices ones(size(FV.vertices,1),1)]';
            FV.vertices = tf(1:3,:)';
            
            % Flip the normals
            FV.faces(:,[1 2]) = FV.faces(:,[2 1]);
            
            stlwrite(filename(c).surface,FV);
            fprintf('Surface saved as %s\n',filename(c).surface)
            clear FV
        end
        
        % Check if bounding box of both images are the same. If not, the mask
        % that was outlined on the anatomical images needs to be cropped to the
        % dimensions of the DTI scans.
        
%         BB_anat = mask.hdr.dime.pixdim(2:4) .* mask.hdr.dime.dim(2:4);
%         BB_DTI  = DTI_data.hdr.dime.pixdim(2:4) .* DTI_data.hdr.dime.dim(2:4);
        
        if MakeMasks == true
% -------------------------------------------------------------------------
%             This is the old code that was used before I found out
%             that c3d can resample to the dimensions of a reference image
%             using the reslice-identity command
%             if abs(BB_anat(3) - BB_DTI(3)) > 0.1
%                 % Bounding box of DTI and anatomical mask are different.
%                 % Crop the mask to the dimensions of the DTI scan, assuming
%                 % that the centre of the bounding boxes are aligned.
%                 
%                 nSlices = round(BB_DTI(3) ./ mask.hdr.dime.pixdim(4));
%                 first_slice = round((BB_anat(3) - BB_DTI(3)) / (2*mask.hdr.dime.pixdim(4)))+1;
%                 offset = (BB_anat(3) - BB_DTI(3)) /  (2*mask.hdr.dime.pixdim(4));
%                 
%                 mask_cropped = mask;
%                 mask_cropped.img = mask.img(:,:,first_slice : 1 : (first_slice+nSlices-1));
%                 
%                 mask_cropped.hdr.dime.dim(4) = nSlices;
%                 mask_cropped.hdr.hist.qoffset_z = mask_cropped.hdr.hist.qoffset_z + offset;
%                 mask_cropped.hdr.hist.srow_z(4) = mask_cropped.hdr.hist.srow_z(4) + offset;
%                 
%                 % Overwrite the mask with the cropped version.
%                 save_untouch_nii(mask_cropped,fullfile(tmpdir,'mask.nii.gz'));
%             end
%             
%             commandTxt = sprintf('c3d -int 0 %s -resample %dx%dx%d -o %s',...
%                 fullfile(tmpdir,'mask.nii.gz'),...
%                 siz(1),siz(2),siz(3),...
%                 fullfile(tmpdir,'mask_resampled.nii.gz'));
%             [status,cmdout] = system(commandTxt);
% -------------------------------------------------------------------------
            % Resample mask to DTI dimensions using the -reslice-identity
            % command in convert3D, which resamples an image to the
            % dimensions of a reference image.
            commandTxt = sprintf('c3d %s %s -reslice-identity -o %s',...
                fullfile(tmpdir,'DTI_3D.nii.gz'),...
                fullfile(tmpdir,'mask.nii.gz'),...
                fullfile(tmpdir,'mask_resampled.nii.gz'));
            [status,cmdout] = system(commandTxt);
            
            % Read resampled mask into the workspace
            mask_resampled = load_untouch_nii(fullfile(tmpdir,'mask_resampled.nii.gz'));
            
            % -------------------- TRACTOGRAPHY MASK FILES ----------------------
            
            % Create a new structure that will be written to a nifti file with only
            % the selected label, resampled to the DTI resolution.
            % Three masks are created for each label:
            % - the full mask with all voxels of the selected label
            % - the boundary mask with only voxels at the boundary of the selected
            %   label. This mask is the region of termination for muscle tractography.
            % - the seed mask with only voxels 2 voxels away from the boundary. This
            %   mask is the seed region for muscle tractography.
            
            % Create the full mask
            full_mask       = DTI_data;
            full_mask.img   = cast(mask_resampled.img~=0,'like',DTI_data.img);
            full_mask.hdr.dime.dim(1) = 3;
            full_mask.hdr.dime.dim(5) = 1;
            full_mask.hdr.dime.scl_slope = 1;
            
            % Overwrite some of the header information, because DSI studio doesn't deal
            % well with more complicated transform matrices in the header. The
            % following changes have been checked to work for the leg scans.
            % Note: changing the sign in the srow_x/y/z matrices affects the
            % orientation in DSI Studio, but pixdim does not. DSI studio
            % also does not read any translational data (last element of
            % the srow-vectors).
            full_mask.hdr.hist.srow_x = [sign(DTI_data.hdr.hist.srow_x(1)) * DTI_data.hdr.dime.pixdim(2) 0 0 0];
            full_mask.hdr.hist.srow_y = [0 sign(DTI_data.hdr.hist.srow_y(2)) * DTI_data.hdr.dime.pixdim(3) 0 0];
            full_mask.hdr.hist.srow_z = [0 0 sign(DTI_data.hdr.hist.srow_z(3)) * DTI_data.hdr.dime.pixdim(4) 0];
            
            % Create the boundary mask
            boundary_mask = full_mask;
            %             boundary_mask.img = cast(bwperim(boundary_mask.img,4),'like',boundary_mask.img);
            boundary_mask.img = cast(~boundary_mask.img,'like',boundary_mask.img);
            
            % Remove n voxels from the boundary to make a seed mask
            seed_mask = full_mask;
            for k = 1:p.Results.NumberOfVoxelsToRemove
                seed_mask.img = seed_mask.img - cast(bwperim(seed_mask.img,4),'like',seed_mask.img);
            end
            
            % Create nifti-files for all three masks
            % Save the data in the same folder as the segmentation file
            filename(c).full_mask     = fullfile(p.Results.ResultsPath, sprintf('%s_%s_full.nii.gz',MaskPrefix,CurrentLabelName));
            filename(c).boundary_mask = fullfile(p.Results.ResultsPath, sprintf('%s_%s_boundary.nii.gz',MaskPrefix,CurrentLabelName));
            filename(c).seed_mask     = fullfile(p.Results.ResultsPath, sprintf('%s_%s_seed.nii.gz',MaskPrefix,CurrentLabelName));
            
            save_untouch_nii(full_mask,filename(c).full_mask)
            fprintf('Full mask saved as %s\n',filename(c).full_mask)
            
            save_untouch_nii(boundary_mask,filename(c).boundary_mask)
            fprintf('Boundary mask saved as %s\n',filename(c).boundary_mask)
            
            save_untouch_nii(seed_mask,filename(c).seed_mask)
            fprintf('Seed mask saved as %s\n',filename(c).seed_mask)
        end
    end
    
    % Remove temporary working directory.
    rmdir(tmpdir,'s')
    
    
catch ME
    % remove temporary working directory, then throw error message
    rmdir(tmpdir,'s')
    error(ME.message)
end
