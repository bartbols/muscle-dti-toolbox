function tract2vtk(DTItracts,vtk_filename,varargin)
%%TRACT2VTK writes the data in the structure DTItracts to a .vtk
%polydata file that can be read in by programs like 3D slicer and ParaView
%for nice visualisation of the fibre tracts.
%
% Bart Bolsterlee, Neuroscience Research Australia (NeuRA)
% August 2017
%
% % Acknowledgement: Inspiration taken from the function TransTrk2Vtk.m on
% MATLAB Central by Wyn Mew, Jan 2016.
%
% ----------------- USAGE -----------------
% tract2vtk(DTItracts,vtk_filename,'ParameterName',<value>)
%
% ----------------- INPUT -----------------
%
% - DTItracts : a structure array with the DTItracts. Alternatively, a
%               MATLAB tract-filename (including extension '.mat') can be 
%               provided as well. If a n x 1 structure is provided, the
%               data from all DTItracts are combined in one file.
% - vtk_filename: ASCII vtkPolydata file to which the data will be written.
%                 Include '.vtk' in the filename.
% ----- OPTIONAL -----
% Optional inputs are provided as pairs of 'ParameterName',<value> (e.g.
% 'pct_threshold',40)
%
% - pct_threshold: If provided, all fibres from the DTItracts with percentage
%                  extrapolations larger than this threshold are excluded.
%                  Default = [] (no threshold is applied)
% - mm_threshold:  If provided, all fibres from the DTItracts with
%                  extrapolations in mm larger than this threshold are 
%                  excluded.
%                  Default = [] (no threshold is applied)
% - selection    : List of indices of tracts to be included in the .vtk
%                  file. Default: [].(i.e. all fibres are included, except 
%                  if pct_threshold is defined as well)/
%                  If multiple DTI tracts are provided, selection should be
%                  a cell with the same dimensions as DTItracts in which
%                  each element 
%              
% - ToWrite      : Tract data to be written to vtk-file:
%                  - 'raw': raw tracts
%                  - 'trunc': truncated tracts
%                  - 'poly': polynomial fitted curves incl extensions (i.e.
%                           'fascicles')
%                  Default: 'poly'
% - ColorData    : List of parameter values to be appended as scalars to the
%                  VTK-file. These scalars can be used for color-coding the
%                  fibres in ParaView or 3D slicer. Default: all
%                  architectural and diffusion parameters that are
%                  available in DTItracts.
% - ParaView     : true/false. Open the .vtk file in ParaView. Only works 
%                  if ParaView is installed and added to the path.
%                  Default: false
% - Title        : title of VTK-file. Default: 'Tract data from MATLAB'
%

% Read inputs 

if nargin == 0
    % If no inputs are defined, select a DTI-tracts file and set the
    % VTK-filename interactively.
    
    [FileName,PathName] = uigetfile('*.mat','Select a DTI tracts in MATLAB-format');
    if PathName == 0;return;end
    DTItracts = load(fullfile(PathName,FileName));
    [FileName,PathName] = uiputfile('*.vtk','Set the filename for the VTK-file');
    vtk_filename = fullfile(PathName,FileName);
    if PathName == 0;return;end
end


tic
default_clist = {'fibrelength','penangle','penangle_l','curvature',...
    'pennation1','pennation2','pct_ext','fa','md','lambda1','lambda2','lambda3'};

p = inputParser;
addRequired(p,'DTItracts')
addRequired(p,'vtk_filename',@(x) contains(x,'.vtk'))
addParameter(p,'pct_threshold',[],@isnumeric)
addParameter(p,'mm_threshold',[],@isnumeric)
addParameter(p,'selection',[],@(x) isnumeric(x)||iscell(x))
addParameter(p,'ToWrite','poly',@(x)contains(x,{'poly','raw','trunc'},'IgnoreCase',true))
addParameter(p,'ColorData',default_clist,@(x) iscell(x)||ischar(x))
addParameter(p,'ParaView',false,@(x) islogical(x)||x==0||x==1)
addParameter(p,'Title','Tract data from MATLAB',@(x) ischar(x))
parse(p,DTItracts,vtk_filename,varargin{:});

pct_threshold = p.Results.pct_threshold;
mm_threshold  = p.Results.mm_threshold;
selection     = p.Results.selection;
towrite       = p.Results.ToWrite;
color_list    = p.Results.ColorData;
ParaView      = p.Results.ParaView;
title_txt     = p.Results.Title;

if ischar(color_list)
    color_list = cellstr(color_list);
end

% Check if structure or filename is provided
if ischar(DTItracts)
    % filename is provided
    if exist(DTItracts,'file') == 2
        fprintf('DTItracts loading from: %s ...',DTItracts)
        DTItracts = load(DTItracts);
        fprintf(' completed.\n')
    else
        error('File not found: %s\nNo VTK-file was written.\n', DTItracts)
    end
elseif iscell(DTItracts)
    % list of filenames is provided
    filenames = DTItracts;
    clear DTItracts
    for d = 1 : length(filenames)
        if exist(filenames{d},'file') == 2
            fprintf('DTItracts loading from: %s ...',filenames{d})
            DTItracts(d) = load(filenames{d});
            fprintf(' completed.\n')
        else
            error('File not found: %s\nNo VTK-file was written.\n', filenames{d})
        end
    end
    
end

%% Select data from DTItracts for plotting

fprintf('---------------- tract2vtk --------------------\n')
nFiles = numel(DTItracts); % number of DTItracts
fprintf('Number of tract sets provided: %d\n',nFiles)
SEL = cell(nFiles,1);
for d = 1 : nFiles
    % Check if selection is defined
    if isempty(selection)
        % ... if not defined, include all fibres.
        SEL{d} = 1 : size(DTItracts(d).fibindex,1);
    else
        SEL{d} = selection;
    end
    
    % Exclude fibres based on percentage extrapolations
    if ~isempty(pct_threshold)
        SEL{d}(DTItracts(d).pct_ext(SEL{d}) > pct_threshold | isnan(DTItracts(d).pct_ext(SEL{d})) ) = [];
    end

    % Exclude fibres based on absolute extrapolations
    if ~isempty(mm_threshold)
        SEL{d}(nansum(DTItracts(d).ext(SEL{d},:),2) > mm_threshold) = [];
    end
    
    fprintf('A total of %d fibres (%.1f%% of total) are included for DTItracts(%d).\n',...
        numel(SEL{d}),...
        numel(SEL{d})/size(DTItracts(d).fibindex,1)*100,...
        d)
        
end
%% POINT_DATA
% Make arrays with only the selected datapoints
points = [];
t=0;
fprintf('Fibre type written to file: %s\n',char(towrite))
for d = 1 : nFiles
    nTracts = numel(SEL{d});

    ns = 20; % number of points to sample along the polynomial curve
    for j = 1 : nTracts
        i = SEL{d}(j);
        switch char(towrite)
            case 'raw'
                % get point data from raw fibre tracts
                sel = min(DTItracts(d).fibindex(i,:)):max(DTItracts(d).fibindex(i,:));
                newpoints = DTItracts(d).tracts_xyz(:,sel);
            case 'trunc'
                % get point data from truncated fibre tracts
                sel = min(DTItracts(d).fibindex_trunc(i,:)):max(DTItracts(d).fibindex_trunc(i,:));
                newpoints = DTItracts(d).tracts_xyz(:,sel);
            case 'poly'
                % get point_data from polynomials + extensions
                P = DTItracts(d).PolyCoeff;
                X = [DTItracts(d).endpoints(i,1,1) polyval(P(i).x,linspace(P(i).t0,P(i).t1,ns)) DTItracts(d).endpoints(i,2,1)];
                Y = [DTItracts(d).endpoints(i,1,2) polyval(P(i).y,linspace(P(i).t0,P(i).t1,ns)) DTItracts(d).endpoints(i,2,2)];
                Z = [DTItracts(d).endpoints(i,1,3) polyval(P(i).z,linspace(P(i).t0,P(i).t1,ns)) DTItracts(d).endpoints(i,2,3)];
                newpoints = [X;Y;Z];
        end
        t = t+1;
        len(t) = size(newpoints,2);
        points = [points newpoints];
    end
end
nTracts = t;
nPoints = size(points,2);

%%
% Open a VTK file for writing
fid = fopen(vtk_filename, 'w'); 

% Write the version information 
fprintf(fid, '# vtk DataFile Version 3.0\n');
% Add a title
fprintf(fid, '%s\n',title_txt);
fprintf(fid, 'ASCII\n');
% Define what data this vtk file containts (polydata)
fprintf(fid, 'DATASET POLYDATA\n');

% Set number of points
fprintf(fid, ['POINTS ' num2str(nPoints) ' float\n']);
precision = '3';
% Write the point data
fprintf('Writing point data...')
spec = [repmat(['%0.',precision,'f '],1,3),'\n'];
fprintf(fid,spec, single(points));
fprintf(' completed.\n')

%% Append the connections
fprintf(fid,'\nLINES %d %d\n',nTracts,nPoints + nTracts);
l = -1;
segm_fibre_nr = zeros(nPoints,1);
fprintf('Writing line data...')
for j = 1:nTracts
    f = l+1;
    l = f+len(j)-1;
    spec = [repmat('%d ',1,len(j)+1) '\n'];
    
    % Add the line segments to the vtk file
    fprintf(fid,spec,[len(j) f:l]);
    segm_fibre_nr((f:l)+1,1) = j;
end
fprintf(' completed.\n')

%% Append scalar data for color-coding

% Check if variables in the color list are fields in DTItracts.

% If multiple DTI tracts are combined in one vtk-file, also number the
% fibres with the file number they originate from.
if nFiles > 1
    color_list{end+1} = 'muscle_nr';
end

counter = 0;
for c = 1 : numel(color_list)
    name = color_list{c};    
    if contains(name,{'pennation1','pennation2'})
        if ~isfield(DTItracts(1),'penangle')
            fprintf('%s is not a field in DTItracts. Color data not written.\n',name)
            continue
        end
    else
        if ~isfield(DTItracts(1),name) && ~strcmp(name,'muscle_nr')
            fprintf('%s is not a field in DTItracts. Color data not written.\n',name)
            continue
        end
    end
    fprintf('Writing color data for variable ''%s''... ',name)
    
    V = [];
    for d = 1 : nFiles
        switch name
            case 'fibrelength'
                value_per_fibre = DTItracts(d).fibrelength(SEL{d});
            case 'penangle'
                value_per_fibre = mean(DTItracts(d).penangle(SEL{d},:),2);
            case 'penangle_l'
                value_per_fibre = mean(DTItracts(d).penangle_l(SEL{d},:),2);
            case 'curvature'
                value_per_fibre = DTItracts(d).curvature(SEL{d});
            case 'pennation1'
                value_per_fibre = DTItracts(d).penangle(SEL{d},1);
            case 'pennation2'
                value_per_fibre = DTItracts(d).penangle(SEL{d},2);
            case 'pct_ext'
                value_per_fibre = DTItracts(d).pct_ext(SEL{d});
            case 'fa'
                value_per_fibre = DTItracts(d).fa(SEL{d});
            case 'md'
                value_per_fibre = DTItracts(d).md(SEL{d});
            case 'lambda1'
                value_per_fibre = DTItracts(d).lambda1(SEL{d});
            case 'lambda2'
                value_per_fibre = DTItracts(d).lambda2(SEL{d});
            case 'lambda3'
                value_per_fibre = DTItracts(d).lambda3(SEL{d});
            case 'muscle_nr'
                value_per_fibre = ones(length(SEL{d}),1)*d;
        end
        V = [V;value_per_fibre];
    end
    % Replace NaN's by -999
    V(isnan(V)) = -999;
    color_data = zeros(nPoints,1);
    for j = 1:nTracts
        % Specify the same color for all elements in this fibre
        color_data(segm_fibre_nr == j) = V(j);
    end
    nScalars = 1;
    scalar_value = color_data;
%     spec  = [repmat(['%0.',precision,'f '],1,nScalars) ' \n'];
    if strcmp(name,'muscle_nr')
        spec = '%d\n';
        cl = 'int';
    else
        spec  = '%.3f\n';
        cl = 'float';
    end
    if counter == 0
        % Only write 'POINT_DATA' once
        fprintf(fid, 'POINT_DATA %d\n',nPoints);
    end
%     fprintf(fid, 'SCALARS %s float %d\n',name,nScalars); %ASCII header
    fprintf(fid, 'SCALARS %s %s\n',name,cl); %ASCII header
    fprintf(fid, 'LOOKUP_TABLE default\n'); %ASCII header
    fprintf(fid,spec,scalar_value);
    fprintf(' completed.\n')
    counter = counter + 1;
end
% Close the vtk file
fclose(fid);
fprintf('\nFibres succesfully written to %s\n',vtk_filename)
fprintf('File contains %d fibres, %d points and %d scalar indices\n',nTracts,nPoints,counter)

% Open the file in ParaView, if ParaView is installed and added to the path
if ParaView == true
    try
        cmdtxt = sprintf('paraview --data=%s &',vtk_filename);
        system(cmdtxt)
    catch ME
        error(ME.message)

    end
end

t_elapsed = toc;
fprintf('Total time to write vtk-file: %.1f seconds\n',t_elapsed)
fprintf('------------------------------------------------\n')
end % of function
