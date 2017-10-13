function varargout = Preprocessing_and_DTI_recon( varargin )
%PREPROCESSING_AND_DTI_RECON This function filters the raw DTI data,
%corrects the bvec-file and calls DSI studio to reconstruct the diffusion
%tensor. The data is then saved.
%
% Bart Bolsterlee, Neuroscience Research Australia (NeuRA)
% February 2017
%
% Note: this function requires DSI studio to be installed on the pc and 
% added to the path so that 'dsi_studio' is recognised as an external command. 
% DSI Studio can be downloaded here:
% http://dsi-studio.labsolver.org/dsi-studio-download
%
% ----------------- USAGE ----------------- 
% filename = Preprocessing_and_DTI_recon(varargin)
% 
% ----------------- INPUT ----------------- 
% This function only has optional input arguments, provided as
% 'ParameterName',<value> pairs. If no inputs are provided, the user has
% to select filenames for the DTI-, bval- and bvec-file and is asked
% whether to filter the data and whether to correct the BVEC file. The
% DSI Studio SRC and FIB files are created with the default filenames (same
% name as DTI file but with extension .src.gz and .fib.gz, respectively).
%
% ----- OPTIONAL -----
% DTI: name of the 4D .nii.gz file containing the DTI data. If not
% provided, the user is requested to select a file.
%
% bvec: name of the file containing the b-vector data. If not
% provided, the user is requested to select a file.
%
% bval: name of the file containing the b-value data. If not
% provided, the user is requested to select a file.
%
% ResultsPath: path where results (filtered DTI, src and fib-files, etc)
%              will be saved to.
%
% SRC: name of the DSI-studio source file. If not provided, the same name
% as the DTI file is used (but with extension .src.gz). 
%
% FIB: name of the DSI-studio fibre file. If not provided, the same name
% as the DTI file is used (but with extension .fib.gz). 
%
% filter: if true, the DTI data is filtered with a LPCA filter. Otherwise,
% the raw DTI data is used. Default is that data is filtered.
%
% CorrectBVEC: if true, the bvec-file is corrected for oblique image
% acquisition using the transform matrix stored in the NIFTI header of the
% DTI file. If false, the provided .bvec file is not corrected. Default = true.
%
% register: if true, the DTI data will be registered to a reference scan.
% The following parameters are then also required:
%
% anat: filename of anatomical reference scan
% parfile: the elastix parameter file used for registration
%
% And these parameters are optional parameters
%
% - mask                 : filename of the mask file used for registration
% - foreground_threshold : threshold intensity for foreground. A foreground 
%                          mask will be created using this threshold.
% - stack                : if anatomical scan is 4D, choose which stack is
%                          used for registration.
% - b0_stack             : stack numbers in DTI file to which the
%                          anatomical scan will be registered. Default = 1 
%                          (usually the b0-image is the first image in the stack)

%
% ----------------- OUTPUT ----------------- 
% filename: structure array with the filenames of the files used for tensor
%           reconstruction.
% Example:
% filename = Preprocessing_and_DTI_recon('DTI','DTI_data.nii.gz',...
%                'bval',DTI_data.bval,'bvec','DTI_data.bvec',...
%                'filter',true,'CorrectBVEC',true)
%
% Example with registration:
% filename = Preprocessing_and_DTI_recon('DTI','DTI_data.nii.gz',...
%                'bval',DTI_data.bval,'bvec','DTI_data.bvec',...
%                'filter',true,'CorrectBVEC',true,...
%                'register',true,'anat','anatomical.nii.gz',...
%                'parfile','affine_parameters.txt')

% Read input arguments
p = inputParser;
addParameter(p,'DTI',[],@(x) contains(x,'.nii.gz'))
addParameter(p,'bval',[],@(x) contains(x,'.bval'))
addParameter(p,'bvec',[],@(x) contains(x,'.bvec'))
addParameter(p,'SRC',[],@(x) contains(x,'.src.gz'))
addParameter(p,'FIB',[],@(x) contains(x,'.fib.gz'))
addParameter(p,'filter',true,@(x) x==0 || x==1 || islogical(x) )
addParameter(p,'CorrectBVEC',true,@(x) x==0 || x==1 || islogical(x) )
addParameter(p,'ResultsPath',[])

% Registration parameters
addParameter(p,'register',false,@(x) x==0 || x==1 || islogical(x) )
addParameter(p,'anat',[],@(x) contains(x,'.nii.gz'))
addParameter(p,'parfile',[],@(x) exist(x,'file')==2)
addParameter(p,'mask',[],@(x) contains(x,'.nii.gz'))
addParameter(p,'foreground_threshold',10,@(x) assert(isscalar(x)))
addParameter(p,'stack',[],@(x) assert(isscalar(x)))
addParameter(p,'b0_stack',1,@(x) assert(isscalar(x)))

parse(p,varargin{:});

% Make cell structure with input arguments
F(:,1) = fieldnames(p.Results);
F(:,2) = struct2cell(p.Results);


%% Add parameters that are not provided
if isempty(F{strcmp(F(:,1),'DTI'),2})
%    No DTI file is provided. Select a file now.
    [fname,path] = uigetfile({'*.nii.gz'},'Select the unfiltered DTI data');
    F{strcmp(F(:,1),'DTI'),2} = fullfile(path,fname);
end
  
if isempty(F{strcmp(F(:,1),'bval'),2})
%   No bval file is provided. Select a file now.
    if exist('path','var')
        [fname,path] = uigetfile({'*.bval'},'Select a bval file',path);
    else
        [fname,path] = uigetfile({'*.bval'},'Select a bval file');
    end
    F{strcmp(F(:,1),'bval'),2} = fullfile(path,fname);
end  
  
if isempty(F{strcmp(F(:,1),'bvec'),2})
%   No bvec file is provided. Select a file now.
    if exist('path','var')
        [fname,path] = uigetfile({'*.bvec'},'Select a bvec file',path);
    else
        [fname,path] = uigetfile({'*.bvec'},'Select a bvec file');
    end
    F{strcmp(F(:,1),'bvec'),2} = fullfile(path,fname);
end  

if isempty(F{strcmp(F(:,1),'ResultsPath'),2})
%   No ResultsPath is provided. Select a path now.
    if exist('path','var')
        ResultsPath = uigetdir(path,'Select a directory where the results will be saved to');
    else
        ResultsPath = uigetdir(pwd,'Select a directory where the results will be saved to');
    end
    F{strcmp(F(:,1),'ResultsPath'),2} = ResultsPath;
end  

% Create the results path, if it doesn't exist already
ResultsPath = F{strcmp(F(:,1),'ResultsPath'),2};
if exist(ResultsPath,'dir') ~= 2
    mkdir(ResultsPath)
    fprintf('Results directory created: %s\n',ResultsPath)
end% Decide whether to filter or not

FilterFlag      = F{strcmp(F(:,1),'filter'),2};
if nargin == 0
    answer = questdlg('Do you want to filter the data?','To filter or not to filter',...
        'Yes','No','Cancel','Yes');
    switch answer
        case 'Yes'
            FilterFlag = true;
        case 'No'
            FilterFlag = false;
        case 'Cancel'
            fprintf('Preprocessing cancelled\n');
            return;
    end
end
if FilterFlag == true
    app = '_LPCA';
else
    app = '';
end

RegistrationFlag = F{strcmp(F(:,1),'register'),2};
if RegistrationFlag == true
    app2 = '_reg';
else
    app2 = '';
end

% Decide whether to correct the bvec file or not
BVEC_correction = F{strcmp(F(:,1),'CorrectBVEC'),2};
if nargin == 0
    answer = questdlg('Do you want to correct the bvec-file?','BVEC correction?',...
        'Yes','No','Cancel','Yes');
    switch answer
        case 'Yes'
            BVEC_correction = true;
        case 'No'
            BVEC_correction = false;
        case 'Cancel'
            fprintf('Preprocessing cancelled\n');
            return;
    end
end

if isempty(F{strcmp(F(:,1),'SRC'),2})
%    No name for the src-file is provided. Use name of DTI file but change
%    extension to .src.gz
    [path,file,ext] = fileparts(F{strcmp(F(:,1),'DTI'),2});
    F{strcmp(F(:,1),'SRC'),2} = fullfile(F{strcmp(F(:,1),'ResultsPath'),2},[file(1:end-4) app app2 '.src.gz']);
end  
 
if isempty(F{strcmp(F(:,1),'FIB'),2})
%    No name for the fib-file is provided. Use name of DTI file but change
%    extension to .fib.gz
    [path,file,ext] = fileparts(F{strcmp(F(:,1),'DTI'),2});
    F{strcmp(F(:,1),'FIB'),2} = fullfile(F{strcmp(F(:,1),'ResultsPath'),2},[file(1:end-4) app app2 '.fib.gz']);
end  

% Put the filenames in a structure
filename.DTI_raw = F{strcmp(F(:,1),'DTI'),2};
if FilterFlag == true
    [path,file,ext] = fileparts(F{strcmp(F(:,1),'DTI'),2});
    filename.DTI_filt = fullfile(F{strcmp(F(:,1),'ResultsPath'),2},[file(1:end-4) app '.nii.gz']);
else
    filename.DTI_filt = [];
end
filename.bval      = F{strcmp(F(:,1),'bval'),2};
filename.bvec      = F{strcmp(F(:,1),'bvec'),2};
if BVEC_correction == true
    [path,file,ext] = fileparts(filename.bvec);
    filename.bvec_corr = fullfile(F{strcmp(F(:,1),'ResultsPath'),2},[file '_corr' ext]);
end
filename.SRC       = F{strcmp(F(:,1),'SRC'),2};
filename.FIB       = F{strcmp(F(:,1),'FIB'),2};
    
    %% ---- Filter the DTI data ----
    DWI_data = load_untouch_nii(filename.DTI_raw);
    if FilterFlag == true
    
        % **************************************************************************
        % * The LPCA filter is described in:                                       *
        % *                                                                        *
        % * J. V. Manjon, P. Coupe, L. Concha, A. Buades, D. L. Collins, M. Robles *
        % * Diffusion Weighted Image Denoising using overcomplete Local PCA.       *
        % * PLoS ONE 8(9): e73021.                                                 *
        % **************************************************************************

        % Filtering DTI data (LPCA method)
        % --- Notes for using LPCA filtering from personal correspondence with
        % Pierrick Coupe on usage of DWIDenoisingLPCA.p ----- 
        % DWIdenoised = DWIDenoisingLPCA(ima, 1, rician, nbthread, verbose);
        % ima is the 4D DWI to denoise
        % rician = 0 for gaussian or 1 for rician noise.
        % nbthread depends on your computer.
        % verbose = 0;

        % Apply filter with the following settings
        rician   = 1; % if 1, rician noise supprression
        nbthread = 1;
        verbose  = 1;
        DWIdenoised = DWIDenoisingLPCA(DWI_data.img, 1, rician, nbthread, verbose);
        
        % Save the filtered data as a nifti-file again.
        % Make NaN's zeros
        DWIdenoised(isnan(DWIdenoised)) = 0;
        
        DWI_data.img = cast(DWIdenoised,'like',DWI_data.img);
        save_untouch_nii(DWI_data,filename.DTI_filt)
    else
        warning('DTI data is not filtered')
    end
 
    %% ---- Registration ---
    if RegistrationFlag == true
        if FilterFlag == true
            filename.dti_before_reg = filename.DTI_filt;
        else
            filename.dti_before_reg = filename.DTI_raw;
        end
        filename.DTI_reg = [filename.dti_before_reg(1:end-7) '_reg.nii.gz'];
       
        % Register data to the anatomical scan using the provided settings
        anat                 = F{strcmp(F(:,1),'anat'),2};
        parfile              = F{strcmp(F(:,1),'parfile'),2};
        stack                = F{strcmp(F(:,1),'stack'),2};
        b0_stack             = F{strcmp(F(:,1),'b0_stack'),2};
        foreground_threshold = F{strcmp(F(:,1),'foreground_threshold'),2};
        mask                 = F{strcmp(F(:,1),'mask'),2};
        
        register_DTI_to_anat(filename.dti_before_reg,anat,...
            parfile,filename.DTI_reg,...
            'stack',stack,'b0_stack',b0_stack,...
            'foreground_threshold',foreground_threshold,...
            'mask',mask)
    end

    
    %% ---- Correct the bvec-file ----
    if BVEC_correction == true
        [T,bvec_corr] = correct_bvec(filename.DTI_raw,filename.bvec,filename.bvec_corr);
    end   
    %% ---- Create src file  ----
    
%     Make SRC file with DSI Studio
    if BVEC_correction == true
        bvec_file = filename.bvec_corr;
    else
        bvec_file = filename.bvec;
    end
    
    % Decide which DTI data file to use for tensor reconstruction
    if RegistrationFlag == true
        % The registered file. (This can be filtered or not filtered.)
        DTI_fname = filename.DTI_reg;
    else
        if FilterFlag == true
            % The filtered DTI file (not registered)
            DTI_fname = filename.DTI_filt;
        else
            % The raw DTI file (not registered)
            DTI_fname = filename.DTI_raw;
        end
    end
    
    % Sometimes, DSI Studio does not correctly interpret image angulation 
    % from the header, causing dimensions to be switched in the SRC file. 
    % To prevent this from happening, the header information is altered
    % here before the SRC-file and FIB-file are created.
    
    DWI_uncorr = load_untouch_nii(DTI_fname);
    
    % Correct header
    DWI_uncorr.hdr.hist.srow_x = [sign(DWI_uncorr.hdr.hist.srow_x(1))*DWI_uncorr.hdr.dime.pixdim(2) 0 0 0];
    DWI_uncorr.hdr.hist.srow_y = [0 sign(DWI_uncorr.hdr.hist.srow_y(2))*DWI_uncorr.hdr.dime.pixdim(3) 0 0];
    DWI_uncorr.hdr.hist.srow_z = [0 0 sign(DWI_uncorr.hdr.hist.srow_z(3))*DWI_uncorr.hdr.dime.pixdim(4) 0];
    
    try
        % Save as new file
        char_list = char(['a':'z' '0':'9']) ;
        tmp_file = fullfile(pwd,[char_list(ceil(length(char_list)*rand(1,8))) '.nii.gz']);
        save_untouch_nii(DWI_uncorr,tmp_file);
        
        % Make SRC file with DSI Studio
        commandTxt = sprintf('dsi_studio --action=src --source=%s --bval=%s --bvec=%s --output=%s',...
            tmp_file,...
            filename.bval,...
            bvec_file,...
            filename.SRC);
        %
        [status,cmdout] = system(commandTxt,'-echo');
        delete(tmp_file) % delete temporary DTI file again
    catch ME
        % To prevent the temporary file to remain in existence if an
        % error occurs, first remove temporary file, then throw error message.
        delete(tmp_file)
        error(ME.message)
    end
    
%% ---- Create fib file  ----    
    % Before fibre reconstruction, a mask with only 1's with the
    % dimensions of the DTI scans needs to be created to avoid automatic
    % thresholding.
    mask = DWI_data;
    dim = size(DWI_data.img);
    mask.img = cast(ones(dim(1:3)),'like',DWI_data.img);
    voxel_size = DWI_data.hdr.dime.pixdim(2:4);
    mask.hdr.dime.dim(1) = 3;
    mask.hdr.dime.dim(5) = 1;
    mask.hdr.hist.srow_x = [voxel_size(1) 0 0 0];
    mask.hdr.hist.srow_y = [0 voxel_size(2) 0 0];
    mask.hdr.hist.srow_z = [0 0 voxel_size(3) 0];
    tmp_mask_fname = fullfile(pwd,'tmp.nii.gz');
    save_untouch_nii(mask,tmp_mask_fname);
    
    commandTxt = sprintf('dsi_studio --action=rec --source=%s --method=%d --mask=%s --check_btable=0 --output_tensor=1 --output=%s',...
        filename.SRC,1,tmp_mask_fname,filename.FIB);
    
    [status,cmdout] = system(commandTxt,'-echo');    delete(tmp_mask_fname); % delete the mask file again
    
    % Rename the files generated by DSI studio
    movefile([filename.SRC '.dti.fib.gz'],filename.FIB);

    %% EV1 map with primary eigenvector data 
    % Create an image with primary eigenvector data, which can be loaded as a EV1 map in ITK-snap
    FA_threshold = [0.05 0.5];
    [~,b,~] = fileparts(filename.FIB);
    filename.EV1 = fullfile(ResultsPath,[b(1:end-4) '_EV1.nii.gz']);
    MakeEV1map(filename.FIB,...
        DTI_fname,...
        FA_threshold,...
        filename.EV1);
    
    %% Report to command window
    fprintf('DTI preprocessing completed. The following files were created:\n')
    if BVEC_correction == true
        fprintf('%-30s: %s\n','bvec file',filename.bvec_corr)
    else
        fprintf('%-30s: %s\n','bvec file',filename.bvec)
    end
    fprintf('%-30s: %s\n','FIB-file',filename.FIB)
    fprintf('%-30s: %s\n','SRC-file',filename.SRC)
    fprintf('%-30s: %s\n','Eigenvector1-map',filename.EV1)
    
    if nargout == 1
        varargout{1} = filename;
    end
end

