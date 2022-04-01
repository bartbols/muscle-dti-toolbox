function filename = segm2stl(segm_filename, varargin )
%SEGM2SURFACE Creates surface meshes (.stl file) from a label map/
% segmentation file.
%
% Bart Bolsterlee, Neuroscience Research Australia (NeuRA)
% April 2022 (largely based on MakeSurfaceAndMasks.m)
%
% Note: this function  requires Convert3D to be installed on the pc and
% added to the path so that 'c3d' is recognised as an external command.
% Convert3D can be downloaded here:
% http://www.itksnap.org/pmwiki/pmwiki.php?n=Downloads.C3D
%
% ----------------- USAGE -----------------
% filename = segm2stl(segm_filename,varargin)
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
% ----- OPTIONAL -----
% Optional inputs are provided as pairs of 'ParameterName',<value> (e.g.
% 'LabelNumbers',[1 3 5])
%
% - LabelNumbers: list of label numbers in the segmentation file to create
% mask files for.
% - LabelNames: cell string with list of label names to create mask files for.
% - MaskPrefix: prefix for the mask and surface model filename. The short
% name of the label (e.g. MG) will be added to the prefix. The default is
% the segmentation filename.
% - ResultsPath: folder in which the mask and surface file will be saved. The
%   default is the folder in which the segmentation file lives.
% - MakeSurface: if true, surfaces (.stl) are created. Default = true.
% - FillHoles  : if true, all holes in the binary mask will be filled prior
%                to surface reconstruction. Default = true.
% - FlipNormals : if true, the normal vectors of the surface models are
%                 flipped. Use this option if the normals are pointing
%                 inside. Default: false
% - ResampleRes : isotropic voxelsize used for resampling the mask to. This
%                 parameters determines the edgelength of the triangles in
%                 the model. Default : 1.5.
%
%
% Note:
% If LabelNumbers and LabelNames are not provided, surfaces
% will be created for all labels in the .lab file.
%
% If both LabelNumbers and LabelNames are provided, the masks and surfaces
% are created for the labels in LabelNumbers.

% ----------------- OUTPUT -----------------
% - filename: a n x 1 struct with field "surface" containing the filenames
% of the surface models that were created. n is the number of labels
% for which files were created.
%
% Example 1 (using LabelNames):
% filename = segm2stl('segmentation.nii.gz', ...
% 'LabelNames',{'MG,'LG'})
%
% Example 2 (using LabelNumbers):
% filename = segm2stl('segmentation.nii.gz', ...
% 'LabelNumbers',[2 4],'ResultsPath','MyMasks')


%% Check input arguments
p = inputParser;
addRequired(p,'segm_filename',@(x) contains(x,'.nii'))
addParameter(p,'LabelNumbers',[],@isnumeric)
addParameter(p,'LabelNames',[],@(x) assert(ischar(x) || iscell(x)))
addParameter(p,'MaskPrefix',[],@ischar)
addParameter(p,'ResultsPath',fileparts(segm_filename),@ischar)
addParameter(p,'FillHoles',true,@(x) x==0 || x==1 || islogical(x) )
addParameter(p,'FlipNormals',false,@(x) x==0 || x==1 || islogical(x) )
addParameter(p,'ResampleRes',1.5,@isscalar)
parse(p,segm_filename,varargin{:});

res         = p.Results.ResampleRes;

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
    % Load the label map.
    label_img = load_untouch_nii(segm_filename);
    %
    % Read in the label information file
    if endsWith(segm_filename,'.gz')
        LabelInfoFilename = [segm_filename(1:end-7) '.lab'];
    elseif endsWith(segm_filename,'.nii')
        LabelInfoFilename = [segm_filename(1:end-4) '.lab'];
    else
        error('Unknown file type for segmentation file. The file should have extension .nii or .nii.gz.')
    end
    
    if exist(LabelInfoFilename,'file') == 2 % Load the .lab file with label name information
        fid = fopen(LabelInfoFilename,'r');
        LabelInfo = textscan(fid,'%d %s %s','Delimiter',',');
        fclose(fid);
    else % Create label info from segmentation file.
        unique_label_nrs = unique(label_img.img(:));
        unique_label_nrs(unique_label_nrs==0) = []; % Remove 0 label.
        LabelInfo = cell(1,3);
        for ii = 1 : length(unique_label_nrs)
            LabelInfo{1,1}(ii,1) = unique_label_nrs(ii);
            LabelInfo{1,2}{ii,1} = sprintf('%02d',unique_label_nrs(ii));
            LabelInfo{1,3}{ii,1} = sprintf('label_%02d',unique_label_nrs(ii));
        end
    end
    
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
        [~,tmp,ext] = fileparts(segm_filename);
        if strcmp(ext,'.nii')
            MaskPrefix = tmp;
        elseif strcmp(ext,'.gz')
            MaskPrefix = tmp(1:end-4);
        else
            error('Unknown file format %s',ext)
        end
    else
        MaskPrefix = p.Results.MaskPrefix;
    end
    
    % Check if the label numbers are actually present in the segmentation data
    % file
    ExistingLabelNumbers = double(unique(label_img.img(:)));
    idx = find(~ismember(LabelNumToProcess,ExistingLabelNumbers));
    if ~isempty(idx)
        fprintf('The following labels exist in the .lab file but not in the segmentation file.:\n')
        for ii = 1 : length(idx)
            fprintf('%10d (%10s)\n',...
                LabelNumToProcess(idx(ii)),LabelInfo{2}{idx(ii)})
            
        end
        fprintf('The following label numbers are available: ')
        fprintf('%d ',ExistingLabelNumbers)
        fprintf('\n')
        %         error('Label does not exist. Modify the label info file <a href="matlab: open(''%s'')">%s</a>',...
        %             LabelInfoFilename,LabelInfoFilename)
        % Remove non-existing labels
        idxToProcess(idx) = [];
        
    end
    
    %%
    % Now all labels have been confirmed to exist, create mask files and
    % a surface file for each selected label.
    
    c = 0;
    for CurrentIdx = idxToProcess'
        c = c + 1;
        CurrentLabelName = LabelInfo{2}{CurrentIdx};
        CurrentLabelNr   = LabelInfo{1}(CurrentIdx);
        
        fprintf('Processing label %d (%s)...\n\n',CurrentLabelNr,CurrentLabelName);
        
        % ------------- MASK RESAMPLING --------
        
        % Create a new, temporary file with only the selected label.
        mask     = label_img;
        mask.img = cast(label_img.img == CurrentLabelNr,'like',mask.img);
        mask.hdr.dime.scl_slope = 1;
        mask.hdr.dime.scl_inter = 0;
        if p.Results.FillHoles == true
            for slice_nr = 1 : size(mask.img,3)
                mask.img(:,:,slice_nr) = imfill(mask.img(:,:,slice_nr),'holes');
            end
        end
        save_untouch_nii(mask,fullfile(tmpdir,'mask.nii.gz')); % has dimensions of anatomical scan
        
        % Make a file with isotropic dimensions for surface model generation.
        commandTxt = sprintf('c3d -int 3 %s -resample-mm %.2fx%.2fx%.2fmm -o %s',...
            fullfile(tmpdir,'mask.nii.gz'),res,res,res,fullfile(tmpdir,'mask_iso.nii.gz'));
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
        
        [FV.vertices,FV.faces]   = v2s(mask_iso.img,0.5,opt,method);
        [FV.vertices,FV.faces]   = surfreorient(FV.vertices,FV.faces);
        
        % Apply some smoothing
        [conn,connnum,count] = meshconn(FV.faces,size(FV.vertices,1));
        %  FV.vertices = smoothsurf(FV.vertices,[],conn,10,0.1,'laplacian');
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
        if p.Results.FlipNormals == true
            FV.faces(:,[1 2]) = FV.faces(:,[2 1]);
        end
        stlwrite2(filename(c).surface,FV);
        fprintf('Surface saved as %s\n',filename(c).surface)
        clear FV
    end
    
    % Remove temporary working directory.
    rmdir(tmpdir,'s')
    
    
catch ME
    % remove temporary working directory, then throw error message
    rmdir(tmpdir,'s')
    error(ME.message)
end
