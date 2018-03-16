function [ DTItracts, StopFlag] = TrackFibres( filename,TrackSettings,varargin )
%%TRACKFIBERS calls DSI studio to track fibres with the settings given by
%the structure array TrackSettings.
%
% Note: this function requires DSI studio to be installed on the computer and
% added to the path so that 'dsi_studio' is recognised as an external command.
% DSI Studio can be downloaded here:
% http://dsi-studio.labsolver.org/dsi-studio-download
%
% INPUT:
% 1) The structure array 'filename' that at least contains the fields:
%    - FIB        : name of DSI studio fibre file (.fib or .fib.gz)
%    - DTI        : name of DTI data (with original header)(.nii or .nii.gz)
%    - Tracts     : full filename of the file to which tracts are saved
%    - Seed       : filename of the seed file
%    Optional fields:
%    - ROI1 : filename of the first region of interest file
%    - ROI2 : filename of a second region of interest file
%    - ROA  : filename of the region of avoidance file
%    - ROA2 : filename of second region of avoidance file
%    - TER  : filename of terminative region file
%
%  2) The struct 'TrackSettings' that contains the following fields:
%    Required fields:
%    - FiberCount  : maximum number of fibers.
%
%    Optional fields:
%    - MinLength   : mininum length (in mm) of the tracts. Default: 20
%    - MaxLength   : maximum length (in mm) of the tracts. Default: 200
%    - FA_threshold: minimum threshold for FA; tracts can only pass through
%                    voxels above this threshold. Default: 0.1
%    - StepSize    : step size for propagating tracts (in mm). Default: 1
%    - SeedCount   : maximum number of seeds (doesn't need to be set if
%                    FiberCount is already defined)
%    - MaxAngle    : maximum angle in degrees between subsequent tract segments.
%                    Default = 10
%    - Smoothing   : smoothing factor (between 0-1 where 0 = no smoothing
%                    and 1 = heavy smoothing. Default = 0
%    - MaxTime     : maximum time (in sec) that DSI studio is allowed to
%                    search for fibers. After MaxTime seconds, the search
%                    is terminated: no tracts will be saved and the
%                    StopFlag is set to 0. Default 10 sec.
%
% Optional inputs, provided as parameter_name,<value> pairs.
%    - seed_file    : filename of seed locations in DSI Studio format. See
%                     MakeSeeds.m for further details.
%    - seed_spacing : if provided, non-random seeding will be applied
%                     within the seed region. A grid with spacing in mm
%                     indicated by the value will be created within the
%                     seed region. For example, adding 'seed_spacing', 2
%                     will generate a regular seed grid with spacing 2mm
%                     within the specified seed region.
%                     E.g. 'seed_spacing',2 will generate a regular grid
%                     of seed points with 2 mm spacing within the seed
%                     region.
%
%    For fields that are not present, default values are used. See default
%    values below.
%
% OUTPUT: 1) A struct 'DTItracts' containing the fields 'length', 'tracts'
%            and 't_elapsed':
%            - tracts: points on the tract segments in voxel coordinates
%            - length: the indices of the points in 'tracts' that are
%            connected. See DSI studio documentation for a description of
%            how 'length' and 'tracts' are formatted.
%            - t_elapsed: time it took to find the tracts (in seconds).
%            - tracts_xyz: points on the tract segments in global
%            coordinates
%            - length_mm: tract lengths in millimetres
%
%          2) A logical 'StopFlag': false (0) means no fibres were found
%             in the maximum time and true (1) means fibre tracking was
%             successful.
%
% More information can be found in the online DSI Studio documentation:
% http://dsi-studio.labsolver.org/Manual
%
% -------------------------------------------------------------------------

% Version 1.0
% Author: Bart Bolsterlee, Neuroscience Research Australia
% Date: June 25th, 2015
%
% Change log:
% - BB, 09/02/2017: added conversion from voxel coordinates to global
% coordinates. DTI voxelsize can either be provided as field 'VoxelSize' in
% the input structure TrackSettings, or it will be read from the source
% file (fibre file). The fields tracts_xyz and length_mm are
% added to DTItracts.
% - BB 27/02/2017:
%     - changed input field in filename from .Source to .fib
%     - Replaced fields OutputDir and OutputName in filename by the full
%       tract filename in field 'tracts'
% BB 31/03/2017: changed fieldname of the fibre-filename in structure
% 'filename' from 'fib' to 'FIB'. The old usage with fib in small letters
% is still compatible as well.
% BB 14/08/2017: Calculates the tract point in global coordinates by
% reading in the header information from the DTI file. In addition to the
% fibre-file, the DTI file should now also be provided (as the field 'DTI'
% in the structure 'filename')
% BB 12/12/2017: Added non-random seeding. If an optional extra input
% argument 'seed_spacing' is provided, a regular seed grid with the spacing
% specified by 'seed_spacing' will be generated. Fibre tracking is then
% performed starting from these seed points.
%
% -------------------------------------------------------------------------

% Check the inputs
p = inputParser;
addRequired(p,'filename',@(x) validateattributes(x,{'struct'},{'nonempty'}))
addRequired(p,'TrackSettings',@(x) validateattributes(x,{'struct'},{'nonempty'}))
addParameter(p,'seed_spacing',[],@isscalar)
addParameter(p,'seed_file',[])
parse(p,filename,TrackSettings,varargin{:})
seed_spacing = p.Results.seed_spacing;
seed_file    = p.Results.seed_file;

if ~isempty(seed_spacing) && ~isempty(seed_file)
    error('seed_spacing and seed_file cannot both be defined. Remove one of the inputs.')
end

% Set some default track settings which are used when they are not defined
% in TrackSettings
default.MinLength    = 20;
default.MaxLength    = 200;
default.FA_threshold = 0.1;
default.StepSize     = 1;
default.MaxAngle     = 10;
default.Smoothing    = 0;
default.MaxTime      = 30;

%% ---------------------------------------------------------------------
% ---------------------- DSI STUDIO COMMAND ---------------------------
% ---------------------------------------------------------------------
% Build up the command that will be used to run DSI studio

% Start with the basic command to call DSI studio and ask it to track
% fibres.
CommandTxt = horzcat('dsi_studio --action=trk');

% Change input .fib to .FIB, if only .fib (with small letters) is provided.
if isfield(filename,'fib') && ~isfield(filename,'FIB')
    filename.FIB = filename.fib;
    filename = rmfield(filename,'fib');
end
% add the fibre file name (DSI-studio .fib-file)
if isfield(filename,'FIB')
    if exist(filename.FIB,'file') ~= 2
        error('prog:input',...
            'fib file ''%s'' does not exist. Tracking cannot be started.',...
            filename.FIB)
    else
        % Add fibre filename to command
        fprintf('%-20s: %s\n','fib',filename.FIB)
        CommandTxt = horzcat(CommandTxt,[' --source=' filename.FIB]);
    end
else
    error('fib file not defined. Tracking cannot be started.')
end

% Check tract filename
[OutputDir,fname,ext] = fileparts(filename.Tracts);
% Check extension
assert(strcmp(ext,'.mat'),sprintf('Tract filename must have extension .mat, but has extentension %s.',ext));

% Create the output folder if it doesn't exist yet.
if ~exist(OutputDir,'dir')
    fprintf('Output directory did not exist yet. Creating directory %s\n',OutputDir)
    mkdir(OutputDir)
end

% Add the output filename to the command
if isfield(filename,'Tracts')
    fprintf('%-20s: %s\n','Output',filename.Tracts)
    CommandTxt = horzcat(CommandTxt,[' --output=' filename.Tracts]);
else
    error('Tract filename not defined. Tracking cannot be started.')
end

if exist(filename.Tracts,'file')
    delete(filename.Tracts);
    warning('%s already exists and will be overwritten.',filename.Tracts)
end

% Add the seed filename to the command, if it exists.
if isfield(filename,'Seed') || ~isempty(seed_file)
    
    if isfield(filename,'Seed')
        % Check if seed region file exists
        if exist(filename.Seed,'file') ~= 2
            error('prog:input',...
                'Seed region file ''%s'' does not exist. Tracking cannot be started.',...
                filename.Seed)
        else
            fprintf('%-20s: %s\n', 'Seed region',filename.Seed)
        end
    end
    if ~isempty(seed_file)
        if exist(seed_file,'file') ~= 2
            error('prog:input',...
                'Seed file ''%s'' does not exist. Tracking cannot be started.',...
                seed_file)
        else
            if isfield(filename,'Seed')
                warning('Seed file %s is used for seeding, not region %s',seed_file,filename.Seed)
            end
            fprintf('%-20s: %s\n', 'Seed voxels',seed_file)
        end
    end
    if ~isempty(seed_spacing)
        % Non-random seeding is selected. Generate seed locations
        % as a regular grid with spacing 'seed_spacing' (in mm) within
        % the seed region. Save the seeds with the tracts data
        filename.SeedLocations = strrep(filename.Tracts,'.mat','_seeds.txt');
        fprintf('Creating a regular grid of seeds with spacing %.2fmm within region %s...',...
            seed_spacing,filename.Seed);
        
        seeds  =MakeSeeds(filename.Seed,filename.DTI,seed_spacing,...
            'seed_file',filename.SeedLocations);
        nSeeds = size(seeds,2);
        CommandTxt = horzcat(CommandTxt,[' --seed='    filename.SeedLocations ' --seed_plan=1']);
        fprintf(' completed.\n')
        
        % Remove limit on number of fibres from the track settings.
        if isfield(TrackSettings,'FiberCount');TrackSettings = rmfield(TrackSettings,'FiberCount');end
    elseif ~isempty(seed_file)
        % A file with seed locations is provided.
        CommandTxt = horzcat(CommandTxt,[' --seed='    seed_file ' --seed_plan=1']);
        
        % Remove limit on number of fibres from the track settings.
        if isfield(TrackSettings,'FiberCount');TrackSettings = rmfield(TrackSettings,'FiberCount');end
    else
        % Random seeding within the seed region
        CommandTxt = horzcat(CommandTxt,[' --seed='    filename.Seed ' --seed_plan=0']);
    end
else
    error('Seed region or seed filename not defined. Tracking cannot be started.')
end

% Add the ROI1 filename to the command, if it exists.
fprintf('%-20s: ','ROI1: ')
if isfield(filename,'ROI1')
    fprintf('%s\n', filename.ROI1)
    if exist(filename.ROI1,'file') ~= 2
        fprintf('ROI1 file does not exist. Continue without ROI1.\n')
    else
        CommandTxt = horzcat(CommandTxt,[' --roi='    filename.ROI1]);
    end
else
    fprintf('not defined. Continue without ROI1.\n')
end

% Add the ROI2 filename to the command, if it exists.
fprintf('%-20s: ','ROI2: ')
if isfield(filename,'ROI2')
    fprintf('%s\n', filename.ROI2)
    if exist(filename.ROI2,'file') ~= 2
        fprintf('ROI2 file does not exist. Continue without ROI2.\n')
    else
        CommandTxt = horzcat(CommandTxt,[' --roi2='    filename.ROI2]);
    end
else
    fprintf('not defined. Continue without ROI2.\n')
end

% Add the ROA filename to the command, if it exists.
fprintf('%-20s: ','ROA: ')
if isfield(filename,'ROA')
    fprintf('%s\n', filename.ROA)
    if exist(filename.ROA,'file') ~= 2
        fprintf('ROA file does not exist. Continue without ROA.\n')
    else
        CommandTxt = horzcat(CommandTxt,[' --roa='    filename.ROA]);
    end
else
    fprintf('not defined. Continue without ROA.\n')
end

% Add the second ROA filename to the command, if it exists.
fprintf('%-20s: ','ROA2: ')
if isfield(filename,'ROA2')
    fprintf('%s\n', filename.ROA2)
    if exist(filename.ROA2,'file') ~= 2
        fprintf('ROA2 file does not exist. Continue without ROA2.\n')
    else
        CommandTxt = horzcat(CommandTxt,[' --roa2='    filename.ROA2]);
    end
else
    fprintf('not defined. Continue without ROA2.\n')
end

% Add the TER filename to the command, if it exists.
UseTER2 = false;
fprintf('%-20s: ','TER: ')
if isfield(filename,'TER')
    fprintf('%s\n', filename.TER)
    if exist(filename.TER,'file') ~= 2
        fprintf('TER file does not exist. Continue without TER.\n')
    else
        % Check if a second region of termination is provided. If so,
        % combine into one mask, because DSI studio cannot handle multiple
        % regions of termination.
        if isfield(filename,'TER2')
            fprintf('%-20s: ','TER2: ')
            fprintf('%s\n', filename.TER2)
            if exist(filename.TER,'file') ~= 2
                fprintf('TER2 file does not exist. Only TER is used.\n')
                CommandTxt = horzcat(CommandTxt,[' --ter='    filename.TER]);
            else
                % Combine TER and TER2 into one new mask.
                TER  = load_untouch_nii(filename.TER);
                TER2 = load_untouch_nii(filename.TER2);
                I1 = TER.img;
                I2 = TER2.img;
                if sign(TER2.hdr.hist.srow_x(1)) ~= sign(TER.hdr.hist.srow_x(1))
                    I2 = flip(I2,1);
                end
                if sign(TER2.hdr.hist.srow_y(2)) ~= sign(TER.hdr.hist.srow_y(2))
                    I2 = flip(I2,2);
                end
                TER_comb = TER;
                TER_comb.img = cast(I1 | I2,'like',TER.img);
                char_list = char(['a':'z' '0':'9']) ;
                fname_TER2 = fullfile(tempdir,['TER_combined_' char_list(ceil(length(char_list)*rand(1,8))) '.nii.gz']);
                save_untouch_nii(TER_comb,fname_TER2)
                CommandTxt = horzcat(CommandTxt,[' --ter=' fname_TER2]);
                UseTER2=true;
            end
        else
           % TER2 is not defined. Just use TER as region of termination.
           CommandTxt = horzcat(CommandTxt,[' --ter='    filename.TER]);
        end
    end
else
    fprintf('not defined. Continue without TER.\n')
end

%% ---------------------------------------------------------------------
% ------------------------ TRACK SETTINGS -----------------------------
% ---------------------------------------------------------------------

% Set the minimum track length, if it is defined
fprintf('%-20s: ','Min. length (mm)')
if isfield(TrackSettings,'MinLength')
    fprintf('%.2f\n',TrackSettings.MinLength)
else
    fprintf('not defined. Default of %.2f is used.\n',default.MinLength)
    TrackSettings.MinLength = default.MinLength;
end
CommandTxt = horzcat(CommandTxt,sprintf(' --min_length=%d',TrackSettings.MinLength));

% Set the maximum track length, if it is defined
fprintf('%-20s: ','Max. length (mm)')
if isfield(TrackSettings,'MaxLength')
    fprintf('%.2f\n',TrackSettings.MaxLength)
else
    fprintf('not defined. Default of %.2f is used.\n',default.MaxLength)
    TrackSettings.MaxLength = default.MaxLength;
end
CommandTxt = horzcat(CommandTxt,sprintf(' --max_length=%d',TrackSettings.MaxLength));

% Set the FA threshold, if it is defined
fprintf('%-20s: ','FA threshold (-)')
if isfield(TrackSettings,'FA_threshold')
    fprintf('%.2f\n',TrackSettings.FA_threshold)
else
    fprintf('not defined. Default of %.2f is used.\n',default.FA_threshold)
    TrackSettings.FA_threshold = default.FA_threshold;
end
CommandTxt = horzcat(CommandTxt,sprintf(' --fa_threshold=%-5.3f',TrackSettings.FA_threshold));

% Set the step size, if it is defined
fprintf('%-20s: ','Step size (mm)')
if isfield(TrackSettings,'StepSize')
    fprintf('%.2f\n',TrackSettings.StepSize)
else
    fprintf('not defined. Default of %.2f is used.\n',default.StepSize)
    TrackSettings.StepSize = default.StepSize;
end
CommandTxt = horzcat(CommandTxt,sprintf(' --step_size=%-5.3f',TrackSettings.StepSize));

% Set the number of fibers, if it is defined
fprintf('%-20s: ','Number of fibers')
if isfield(TrackSettings,'FiberCount')
    if ~isempty(TrackSettings.FiberCount)
        fprintf('%d\n',TrackSettings.FiberCount)
        CommandTxt = horzcat(CommandTxt,sprintf(' --fiber_count=%d',TrackSettings.FiberCount));
    else
        fprintf('not defined. There is no limit on the number of fiber tracts.\n')
    end
else
    fprintf('not defined. There is no limit on the number of fiber tracts.\n')
end

% Set the number of seeds, if it is defined
fprintf('%-20s: ','Number of seeds')
if isfield(TrackSettings,'SeedCount')
    if ~isempty(TrackSettings.SeedCount)
        fprintf('%d\n',TrackSettings.SeedCount)
        CommandTxt = horzcat(CommandTxt,sprintf(' --seed_count=%d',TrackSettings.SeedCount));
    else
        fprintf('not defined. There is no limit on the number of seeds.\n')
        
    end
else
    fprintf('not defined. There is no limit on the number of seeds.\n')
end

% Set the maximum angle, if it is defined
fprintf('%-20s: ','Maximum angle (deg)')
if isfield(TrackSettings,'MaxAngle')
    fprintf('%.2f\n',TrackSettings.MaxAngle)
else
    fprintf('not defined. Default of %.2f is used.\n',default.MaxAngle)
    TrackSettings.MaxAngle = default.MaxAngle;
end
CommandTxt = horzcat(CommandTxt,sprintf(' --turning_angle=%d',TrackSettings.MaxAngle));

% Set the smoothing value, if it is not defined
fprintf('%-20s: ','Smoothing')
if isfield(TrackSettings,'Smoothing')
    fprintf('%.2f\n',TrackSettings.Smoothing)
else
    fprintf('not defined. Default of %.2f is used.\n',default.Smoothing)
    TrackSettings.Smoothing = default.Smoothing;
end
CommandTxt = horzcat(CommandTxt,sprintf(' --smoothing=%-5.3f',TrackSettings.Smoothing));

% Set the maximum time, if it is not defined
if isfield(TrackSettings,'MaxTime')
    MaxTime = TrackSettings.MaxTime;
else
    MaxTime = default.MaxTime; % Default maximum time (sec) to find fibres in DSI studio
end
fprintf('Maximum searching time (sec): %.2f\n', MaxTime)

%% ---------------------------------------------------------------------
% ------------------------ CALL DSI STUDIO -----------------------------
% ----------------------------------------------------------------------

% Call DSI studio to track the fibres
fprintf('%s','Calling DSI Studio to track the fibers...')


tic % Start the timer

% Call DSI studio to perform fibre tracking
% The '&' must be added to the command so that DSI studio is running in the
% background.
[status,cmdout] = system(CommandTxt);
disp(cmdout)

% [status,cmdout] = system(CommandTxt);
% if ~isempty(strfind(cmdout,'internal or external command'))
%     error('%s\n%s',cmdout,'Install DSI Studio and add to the path of your operating system. See documentation of the Muscle DTI toolbox for more information.')
%
% end


% The next loop prevents DSI studio from getting stuck when no fibres
% can be found. If the fibre tracking runs for more than 'MaxTime'
% seconds, the program is terminated.
StopFlag = 0;t_elapsed = 0;
while StopFlag == 0 && t_elapsed < MaxTime
    t_elapsed = toc;
    % Check if the tract-file has been created, which is the case if fibre
    % tracking was succesful
    if exist(filename.Tracts,'file') == 2
        % If it has been created, the tracking was successful.
        fprintf(' completed in %6.2f seconds.\n',t_elapsed)
        fprintf('Loading %s... ',filename.Tracts)
        
        % Load the data in the struct 'DTItracts'
        DTItracts                = load(filename.Tracts);
        DTItracts.FibreTrackTime = t_elapsed; % Add the time it took to find the tracts to the struct
        fprintf('completed.\n')
        
        % Set the StopFlag to 1; this will exit the loop.
        StopFlag = 1;
    end
end

if StopFlag == 0
    % Behead DSI Studio!!!
    dos('taskkill /f /im dsi_studio.exe'); % kill DSI studio
    % Create empty fields in 'DTItracts'
    DTItracts.length = [];
    DTItracts.tracts = [];
    DTItracts.FibreTrackTime = t_elapsed;
    warning('Tracts could not be found in %6.2f seconds so the search was terminated.',MaxTime)
else
    
    % Delete temporary TER file
    if UseTER2 == true
        delete(fname_TER2)
    end
    % Add the FileNames and TrackSettings as fields to the DTItracts struct, so
    % that it can always be checked later which settings were used for
    % tractography.
    DTItracts.FileNames     = filename;
    DTItracts.TrackSettings = TrackSettings;
    DTItracts.length        = DTItracts.length';
    DTItracts.TrackSettings.algorithm = 'dsi_studio';
    
    % Add the field 'fibindex' which is a m x 2 array where m is the number
    % of fibers. The first column contains the index of the first point in
    % fiber m, the second column the index of the last point. This
    % information is also stored in 'length', but I find this a more
    % convenient way to store the tracts.
    fibindex(1,1) = 1;
    fibindex(1,2) = DTItracts.length(1);
    for fibnr = 2:length(DTItracts.length)
        fibindex(fibnr,1) = fibindex(fibnr-1,2)+1;
        fibindex(fibnr,2) = sum(DTItracts.length(1:fibnr));
        %     if fibnr == 1
        %         p_start = 1;
        %     else
        %         p_start = DTItracts.length_summed(fibnr-1)+1;
        %     end
        %     p_end = DTItracts.length_summed(fibnr);
    end
    
    DTItracts.fibindex = fibindex;
    
    % Convert tract points from voxel coordinates to global coordinates.
    % Get transformation from voxel to global coordinates from the
    % header of the DTI file.
    DTI = load_untouch_nii(filename.DTI);
    T = [DTI.hdr.hist.srow_x;...
        DTI.hdr.hist.srow_y;...
        DTI.hdr.hist.srow_z;...
        0 0 0 1];
    
    % Flip the second dimension of tracts. I think this is because DSI
    % studio works in a different coordinate system. This may need further
    % checking but works well for the data we have collected so far.
    tracts = DTItracts.tracts;
    tracts(2,:) = (DTI.hdr.dime.dim(3)-1) - tracts(2,:);
    tracts_glob = T * [tracts;ones(1,size(tracts,2))];
    DTItracts.tracts_xyz = tracts_glob(1:3,:);
    
    DTItracts.length_mm   = (DTItracts.length-1) * TrackSettings.StepSize;
    
    save(filename.Tracts,'-struct','DTItracts');
end



end % of the function

