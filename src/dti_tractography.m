function DTItracts = dti_tractography( EV1,seeds,settings,varargin)
%FIBRE_TRACKING_MATLAB performs DTI fibre tractography on the vector field 
% provided in EV1, starting from 'seeds'.
%
% Bart Bolsterlee, Neuroscience Research Australia (NeuRA)
% November 2017
%
% ----------------- USAGE -----------------
% DTItracts = dti_tractography(EV1,seeds,settings )
% or
% DTItracts = dti_tractography(EV1,seeds,settings,...
%                                    'parameter',<value>)
%
% ----------------- INPUT -----------------
% Required inputs:
% EV1      : A NIfTI structure or filename of a NIfTI file containing
%            primary eigenvector data (4D). Alternatively, a DSI
%            Studio-generated FIB file (or filename of the FIB-file) can be
%            provided.
% seeds    : 3 x n array of seed locations in global coordinates
% settings : structure with tractography settings containing the fields:
%            - StepSize: stepsize in mm
%            - MinLength : minimum fibre tract length in mm
%            - MaxLength : maximum fibre tract length in mm
%            - FA_threshold : 2-element vector containing the minimum and
%                            maximum FA value 
%            - MaxAngle  : maximum turning angle in degrees
%
% Optional inputs, provided as 'parameter',<value> pairs:
% DTI      : A NIfTI structure or filename of a NIfTI file containing the
%            DTI data. This field is mandatory when EV1 is a FIB-file from
%            DSI Studio.
% FA       : A NIfTI structure of filename of a NIfTI file containing the
%            fractional anisotropy data. This field is mandatory when EV1
%            is not a FIB-file from DSI Studio.
% TER      : Region of termination. Tracts are terminated whenever they enter
%            this region. The last point on the tract will be inside the TER.
% ROI      : Region of interest. Only tracts that pass through this region
%            are included.
% ROA      : Region of avoidance. Tracts that pass through this region
%            are excluded.
% Save     : MAT-filename to which the structure DTItracts will be saved.
%
% ----------------- OUTPUT -----------------
% DTItracts : structure containing the fields
%             - tracts: 3 x np fibre tract points in voxel coordinates
%             - tracts_xyz : 3 x nF fibre tract points in global coordinates
%             - fibindex : nF x 2 array with the first and last index in
%              'tracts' or 'tracts_xyz' of all nF fibre tracts
%             - StepSize : stepsize used for tracking in mm
%             - TrackSettings: array with the settings used for tracking
%             - seed_idx : index of seed points that resulted in successful tracts.

% Check inputs
p = inputParser;
addRequired(p,'EV1')
addRequired(p,'seeds',@isnumeric)
addRequired(p,'settings')
addParameter(p,'DTI',[])
addParameter(p,'FA',[])
addParameter(p,'TER',[])
addParameter(p,'ROI',[])
addParameter(p,'ROA',[])
addParameter(p,'Save',[])
parse(p,EV1,seeds,settings,varargin{:})

DTI  = p.Results.DTI;
FA   = p.Results.FA;
TER  = p.Results.TER;
ROI  = p.Results.ROI;
ROA  = p.Results.ROA;
savename = p.Results.Save;

if settings.FA_threshold(1) > settings.FA_threshold(2)
    error('Minimum FA threshold (%.2f) is larger than the maximum FA threshold (%.2f). Please change the ''FA_threshold'' property and try again.',...
        settings.FA_threshold(1), settings.FA_threshold(2))
end

%% Check if mex-file for vectorfield interpolation is available.
if exist('interp3vec_mex','file') == 3
    use_mex = true;
else
    warning('The MEX-file interp3vec_mex is not available on the MATLAB path. Fibre tracking may be slow.')
    use_mex = false;
end

%% Load the EV1 data
if isstruct(EV1)
    if isfield(EV1,'img')
        % NIfTI structure is provided.
        type = 'fsl';
    elseif isfield(EV1,'dir0')
        % Structure with DSI Studio fib-file is provided.
        type = 'dsi';
    end
else
    fprintf('Loading %s...',EV1)
    if contains(EV1,'.nii')
        % Read the NIfTI file with the primary eigenvector data
        type = 'fsl';
        EV1 = load_untouch_nii(EV1);
        V = EV1.img;
        voxelsize = EV1.hdr.dime.pixdim(2:4);
    elseif contains(EV1,'.fib')
        % Filename of DSI Studio fib-file is provided.
        type = 'dsi';
        gunzip(EV1);
        FIB = load(EV1(1:end-2),'-mat');
        delete(EV1(1:end-2))
        EV1 = FIB;
        clear FIB
    end
    fprintf(' completed.\n')
end

%% Load the DTI data
if ~isstruct(DTI) && ~isempty(DTI)
    DTI = load_untouch_nii(DTI);
end

%% Load the eigenvector data and store in the appropriate format used for fibre tractography
if strcmp(type,'dsi')
    if isempty(DTI)
        error('If DSI Studio fib is provided, the corresponding DTI data should be provided. Add ''DTI'',<dti_filename> to the input.')
    end
    % Get the primary eigenvector data from the FIB file
    V = permute(reshape(EV1.dir0,[3 EV1.dimension]),[2 3 4 1]);
    voxelsize = EV1.voxel_size;
    
    % Get the transformation from the header of the DTI file
    T = [DTI.hdr.hist.srow_x;...
        DTI.hdr.hist.srow_y;...
        DTI.hdr.hist.srow_z;...
        0 0 0 1];
    
    % Get FA map from the FIB file
    FA = reshape(EV1.fa0,EV1.dimension);
    
    % Flip the primary eigenvector and FA map according to the transformation so
    % that the map is compatible with the expected format in interp3vec.m.
    % (Note: Regardless of the header information in the DTI NIfTI file,
    % the second dimensions needs to be flipped and the sign of the 
    % y-direction of the direction vectors needs to be changed.
    V  = flip(V,2);
    FA = flip(FA,2);
    V(:,:,:,2) = -V(:,:,:,2);
    
elseif strcmp(type,'fsl')
    % Get the primary eigenvector data and transformation from the
    % NIfTI data
    V = EV1.img;
    T = [EV1.hdr.hist.srow_x;...
         EV1.hdr.hist.srow_y;...
         EV1.hdr.hist.srow_z;...
         0 0 0 1];
     % Flip direction and signs according to the header information. This
     % has been checked to be correct for all combinations of negative and
     % positive x- and y-components in T.
     sgn = sign(diag(T(1:3,1:3)));
     if sgn(1) == 1
         V = flip(V,1);
         V(:,:,:,2) = -V(:,:,:,2);
     end
     if sgn(2) == -1
         V = flip(V,2);
     end
     if sgn(1) == -1 && sgn(2) == 1
         V(:,:,:,2) = -V(:,:,:,2);
     end
         
    voxelsize = EV1.hdr.dime.pixdim(2:4);
    % Get FA map
    if isempty(FA)
        error('If fsl data is provided, the FA map should be provided as well. Add ''FA'',<fa_filename> to the input.')
    elseif ~isstruct(FA)
        FA = load_untouch_nii(FA);
        FA = FA.img;
    end
end
clear EV1

%% Load region of termination
if ~isempty(TER)
    if ~isstruct(TER)
        % If the input is not a structure, it should be the filename of a
        % NIfTI mask.
        TER = load_untouch_nii(TER);
    end
    settings.TER = true;
end

%% Load region of interest
if ~isempty(ROI)
    if ~isstruct(ROI)
        % If the input is not a structure, it should be the filename of a
        % NIfTI mask.
        ROI = load_untouch_nii(ROI);
    end
    settings.ROI = true;
end

%% Load region of avoidance
if ~isempty(ROA)
    if ~isstruct(ROA)
        % If the input is not a structure, it should be the filename of a
        % NIfTI mask.
        ROA = load_untouch_nii(ROA);
    end
    settings.ROA = true;
end

%% Transform seeds to voxel space
nSeeds = size(seeds,2);
seeds_vox = T \ [seeds;ones(1,nSeeds)];
seeds_vox = seeds_vox(1:3,:);
maxSteps = round(settings.MaxLength / settings.StepSize); % maximum number of steps in one direction

%% ========================== TRACKING ===============================
tic
fprintf('Performing fibre tractography with %d seeds...',nSeeds)
for start_dir   = [-1 1] % bi-directional tracking

    % Start tracking form the seed
    tracts = NaN(3,nSeeds,maxSteps+1); % uni-directional array with tract points in voxel coordinates
    stepnr = 0;
    seed_incl = 1 : nSeeds; % vector with indices of seeds that are still included in fibre tracking.
    tracts(:,:,1) = seeds_vox;
    while stepnr <= maxSteps && ~isempty(seed_incl)
        stepnr = stepnr + 1;
%         fprintf('Step number %d in direction %d\n',stepnr,start_dir);
        
        % !! the voxel coordinates in tracts are indexed from 0 (not 1) !!
%         use_mex = false;
        if use_mex == true
            d = interp3vec_mex(double(V),double(tracts(:,seed_incl,stepnr)));
        else            
            d = interp3vec(V,tracts(:,seed_incl,stepnr));
        end
        
        if stepnr == 1
            s = ones(1,length(seed_incl)) * start_dir;
        else
            % Compare step direction with previous direction to ensure that the
            % tracts propagate in the same direction as the previous step.
            s = sign(sum(d .* d_prev,1));
        end
        s(s==0) = 1;
        
        % Flip sign of step direction if it is different from previous step
        d = d .* (ones(3,1)*s);
        
        % Calculate direction in voxel dimensions
        d_vox = d ./ (voxelsize' * ones(1,length(seed_incl)));
        
        % Update tract points for next step
        next_step = tracts(:,seed_incl,stepnr) + d_vox * settings.StepSize;
        
        % =============== Check stopping criteria ==================
        % Check if turning angle exceeds the maximum angle
        if stepnr > 1
            turning_angle = acosd(sum(d_prev .* d));
            crit1 = turning_angle < settings.MaxAngle; % true if smaller, false if larger
        else
            crit1 = true(1,nSeeds);
        end
        
        % ============== FA ==================
        % Check if next step is outside FA threshold. If so, this will be
        % the last step for the tract.
        
        % Note that this step also checks whether the next step is outside
        % the image domain. If outside, a NaN will be returned for the FA
        % value causing the tract to be terminated.
        
        fa_value = interp3(FA,...
            next_step(2,:)+1,...
            next_step(1,:)+1,...
            next_step(3,:)+1);
        crit2 = fa_value > settings.FA_threshold(1) & ...
                fa_value < settings.FA_threshold(2);
        
         % ============= region of termination =====================
         % Check if current step is inside region of termination. If so,
         % terminate tracking.
         if ~isempty(TER)
             crit3 = interp3(TER.img,...
                        tracts(2,seed_incl,stepnr)+1,...
                        tracts(1,seed_incl,stepnr)+1,...
                        tracts(3,seed_incl,stepnr)+1,'nearest') == 0;
         else
             crit3 = true(size(crit1));
         end
                 
                 %% Exclude fibres from next step based on stopping criteria
         excluded = ~crit1 | ~crit2 | ~crit3;
         
         % Update seed indices of fibres that are still included
         seed_incl(excluded) = [];
        
        % Add the next step to 'tracts'. Exclude fibres based on the
        % stopping criteria.
        tracts(:,seed_incl,stepnr+1) = next_step(:,~excluded);
        
        % Store current direction so it can be used in the next step
        d_prev = d(:,~excluded); 
    end  
    if start_dir == -1
        tracts_1 = tracts; % results of tracking in direction -1
    else
        tracts_2 = tracts; % Results of tracking in direction +1 (arbitrary, but opposite to direction -1)
    end    
end
t_elapsed = toc;

%%
% figure
% axis equal
% hold on
% tracts_1 = flip(tracts_1,3);
% for ii = 1 : size(tracts_1,3)
%     plot3(tracts_1(1,:,ii),tracts_1(2,:,ii),tracts_1(3,:,ii),'ro')
%     text(tracts_1(1,:,ii),tracts_1(2,:,ii),tracts_1(3,:,ii),int2str(ii))
% end
% 
% for ii = 1 : size(tracts_2,3)
%     plot3(tracts_2(1,:,ii),tracts_2(2,:,ii),tracts_2(3,:,ii),'bo')
%     text(tracts_2(1,:,ii),tracts_2(2,:,ii),tracts_2(3,:,ii),int2str(ii))
% end
% plot3(seeds(1),seeds(2),seeds(3),'go','MarkerSize',8)

%% Combine the results from forwards and backwards tracking
% Store results in the structure DTItracts
DTItracts.tracts = NaN(3,nSeeds * (2*maxSteps+1));
k1=0;
fibnr = 0;

DTItracts.fibindex  = zeros(nSeeds,2);
DTItracts.length_mm = zeros(nSeeds,1);
DTItracts.seed_idx  = zeros(nSeeds,1);


for i = 1 : nSeeds
    combined = [flip(squeeze(tracts_1(:,i,:)),2) squeeze(tracts_2(:,i,2:end))];
    
    % Remove the NaNs
    current_tract = combined(:,~isnan(combined(1,:)));
    
    % If ROI is provided, include only fibres that pass through the ROI.
     if ~isempty(ROI)
         crit4 = interp3(ROI.img,...
                    current_tract(2,:) + 1,...
                    current_tract(1,:) + 1,...
                    current_tract(3,:) + 1,'nearest');
                
            if  ~any(crit4)
                % Tract does not go through ROI. Exclude from results.
                continue
            end
     end
     
    % If ROA is provided, exlude fibres that pass through the ROA.
     if ~isempty(ROA)
         crit5 = interp3(ROA.img,...
                    current_tract(2,:) + 1,...
                    current_tract(1,:) + 1,...
                    current_tract(3,:) + 1,'nearest');
                
            if  any(crit5)
                % Tract goes through ROA. Exclude from results.
                continue
            end
     end
     
     % Check length constraints
    nP = size(current_tract,2);
    length_mm = (nP-1) * settings.StepSize;
    
    if length_mm < settings.MinLength || ...
       length_mm > settings.MaxLength
       % Fibre is longer than the maximum length or shorter than the
       % minimum length. Exclude from results.
        continue
    end
    
    % Add to DTItracts
    fibnr = fibnr + 1;
    DTItracts.seed_idx(fibnr) = i;
    DTItracts.tracts(:,(1:nP)+k1) = current_tract;
    DTItracts.fibindex(fibnr,1) = k1+1;
    DTItracts.fibindex(fibnr,2) = k1+nP;
    DTItracts.length_mm(fibnr) = length_mm;
    k1 = k1+nP;
end

fprintf(' completed.\nIt took %.2f seconds to find %d fibres (%.1f%% of number of seeds).\n',...
    t_elapsed,fibnr,fibnr/nSeeds*100)

DTItracts.seeds      = seeds;
DTItracts.seeds_vox  = seeds_vox;
DTItracts.seed_idx   = DTItracts.seed_idx(1:fibnr);
DTItracts.fibindex   = DTItracts.fibindex(1:fibnr,:);
DTItracts.length_mm  = DTItracts.length_mm(1:fibnr);
DTItracts.tracts     = DTItracts.tracts(:,1:k1); % remove the NaNs;

% Transform tracts to global coordinates
tracts_glob = T * [DTItracts.tracts;ones(1,size(DTItracts.tracts,2) )];
DTItracts.tracts_xyz  = tracts_glob(1:3,:);
DTItracts.StepSize    = settings.StepSize;
DTItracts.TrackSettings = settings;
DTItracts.TrackSettings.algorithm = 'matlab';

%% Save the results
if ~isempty(savename)
    % Appends '.mat' if the extension is not provided.
    if ~endsWith(savename,'.mat')
        savename = [savename '.mat'];        
    end
    % Create directory, if it doesn't exist yet.
    if exist(fileparts(savename),'dir') ~= 7
        mkdir(fileparts(savename));
    end
    save(savename,'-struct','DTItracts')
    fprintf('Tracts saved as %s\n',savename)
end

end

