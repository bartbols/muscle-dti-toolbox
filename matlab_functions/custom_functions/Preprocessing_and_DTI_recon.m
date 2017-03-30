function filename = Preprocessing_and_DTI_recon( varargin )
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
% ----------------- OUTPUT ----------------- 
% filename: structure array with the filenames of the files used for tensor
%           reconstruction.
% Example:
% filename = Preprocessing_and_DTI_recon('DTI','DTI_data.nii.gz',...
%                'bval',DTI_data.bval,'bvec','DTI_data.bvec',...
%                'filter',true,'CorrectBVEC',true)


% Read input arguments
p = inputParser;
addParameter(p,'DTI',[],@(x) ~isempty(strfind(x,'.nii.gz')))
addParameter(p,'bval',[],@(x) ~isempty(strfind(x,'.bval')))
addParameter(p,'bvec',[],@(x) ~isempty(strfind(x,'.bvec')))
addParameter(p,'SRC',[],@(x) ~isempty(strfind(x,'.src.gz')))
addParameter(p,'FIB',[],@(x) ~isempty(strfind(x,'.fib.gz')))
addParameter(p,'filter',true,@(x) x==0 || x==1 || islogical(x) )
addParameter(p,'CorrectBVEC',true,@(x) x==0 || x==1 || islogical(x) )
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

% Decide whether to filter or not
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

% Decide whether to correct the bvec file or noet
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
   F{strcmp(F(:,1),'SRC'),2} = [F{strcmp(F(:,1),'DTI'),2}(1:end-7) app '.src.gz'];
end  
 
if isempty(F{strcmp(F(:,1),'FIB'),2})
%    No name for the fib-file is provided. Use name of DTI file but change
%    extension to .fib.gz
   F{strcmp(F(:,1),'FIB'),2} = [F{strcmp(F(:,1),'DTI'),2}(1:end-7) app '.fib.gz'];
end  

% Put the filenames in a structure
filename.DTI_raw = F{strcmp(F(:,1),'DTI'),2};
if FilterFlag == true
    filename.DTI_filt = [F{strcmp(F(:,1),'DTI'),2}(1:end-7) '_LPCA.nii.gz'];
else
    filename.DTI_filt = [];
end
filename.bval      = F{strcmp(F(:,1),'bval'),2};
filename.bvec      = F{strcmp(F(:,1),'bvec'),2};
if BVEC_correction == true
    filename.bvec_corr = [filename.bvec(1:end-5) '_corr.bvec'];
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
        DWI_data.img = int16(DWIdenoised);
        save_untouch_nii(DWI_data,filename.DTI_filt)
    else
        warning('DTI data is not filtered')
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
    if FilterFlag == true
        DTI_fname = filename.DTI_filt;
    else
        DTI_fname = filename.DTI_raw;
    end
    
    commandTxt = sprintf('dsi_studio --action=src --source=%s --bval=%s --bvec=%s --output=%s',...
        DTI_fname,...
        filename.bval,...
        bvec_file,...
        filename.SRC);
%     
    [status,cmdout] = system(commandTxt,'-echo');
    
%     % Make .src.gz file directly from MATLAB
%     src_new.voxel_size = DWI_data.hdr.dime.pixdim(2:4);
%     src_new.dimension  = DWI_data.hdr.dime.dim(2:4);
%     
%     bval = load(filename.bval,'-ascii');
%     if BVEC_correction == true
%         src_new.b_table = [bval;bvec_corr];
%     else
%         bvec = load(filename.bvec);
%         src_new.b_table = [bval;bvec];
%     end
%     
% %     % Decide which dimensions to flip based on the information in the srow_
% %     % header entries.
% %     % !!! the next has not been tested yet!!!
% %     flipx=false;flipy=false;flipz=false;
% %     bvec_new = bvec_corr;
% %     if DWI_data.hdr.hist.srow_x(1) < 0;flipx = true;bvec_new(1,:) = -bvec_new(1,:);end
% %     if DWI_data.hdr.hist.srow_y(2) < 0;flipy = true;bvec_new(2,:) = -bvec_new(2,:);end
% %     if DWI_data.hdr.hist.srow_z(3) < 0;flipz = true;bvec_new(3,:) = -bvec_new(3,:);end
%     
%     for i = 1 : size(DWI_data.img,4)
%         data = double(DWI_data.img(:,:,:,i));
%         
% %     % !!! this code has not been tested yet!!!
% %         if flipx==true;data = flip(data,1);end
% %         if flipy==true;data = flip(data,2);end
% %         if flipz==true;data = flip(data,3);end
%         
%         eval(['src_new.image' int2str(i-1) '=data(:)'';'])
%     end
% 
%     save(filename.SRC(1:end-3),'-struct','src_new','-v4')
%     gzip(filename.SRC(1:end-3))
%     delete(filename.SRC(1:end-3));
    
%% ---- Create fib file  ----    
    % Before fibre reconstruction, a mask with only 1's with the
    % dimensions of the DTI scans needs to be created to avoid automatic
    % thresholding.
    mask = DWI_data;
    dim = size(DWI_data.img);
    mask.img = int16(ones(dim(1:3)));
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
    if FilterFlag == true
        fname_DTI = filename.DTI_filt;
    else
        fname_DTI = filename.DTI_raw;
    end
    filename.EV1 = MakeEV1map(filename.FIB,...
        fname_DTI,...
        FA_threshold);
    
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
    
    
end

