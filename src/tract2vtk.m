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
%               provided as well.
% - vtk_filename: ASCII vtkPolydata file to which the data will be written.
%                 Include '.vtk' in the filename.
% ----- OPTIONAL -----
% Optional inputs are provided as pairs of 'ParameterName',<value> (e.g.
% 'ext_threshold',40)
%
% - ext_threshold: If provided, all fibres from the DTItracts with percentage
%                  extensions larger than this threshold are excluded.
%                  Default = [] (all fibres are included)
% - selection    : List of indices of tracts to be included in the .vtk
%                  file. Default: [] (all fibres are included)
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

% Read inputs 
default_clist = {'fibrelength','penangle','curvature',...
    'pennation1','pennation2','pct_ext','fa','md','lambda1','lambda2','lambda3'};
p = inputParser;
addRequired(p,'DTItracts')
addRequired(p,'vtk_filename',@(x) contains(x,'.vtk'))
addParameter(p,'ext_threshold',[],@isnumeric)
addParameter(p,'selection',[],@(x) isnumeric(x))
addParameter(p,'ToWrite','poly',@(x)contains(x,{'poly','raw','trunc'},'IgnoreCase',true))
addParameter(p,'ColorData',default_clist,@(x) iscell(x)||ischar(x))
addParameter(p,'ParaView',false,@(x) islogical(x)||x==0||x==1)
parse(p,DTItracts,vtk_filename,varargin{:});

ext_threshold = p.Results.ext_threshold;
selection     = p.Results.selection;
towrite       = p.Results.ToWrite;
color_list    = p.Results.ColorData;
ParaView      = p.Results.ParaView;

if ischar(color_list)
    color_list = cellstr(color_list);
end

% Check if structure or filename is provided
if ~isstruct(DTItracts)
    if exist(DTItracts,'file') == 2
        fprintf('DTItracts loading from: %s ...',DTItracts)
        DTItracts = load(DTItracts);
        fprintf(' completed.\n')
    else
        error('File not found: %s\nNo VTK-file was written.\n', DTItracts)
    end
      
end

%% Select data from DTItracts for plotting
nFiles = length(DTItracts);
SEL = cell(nFiles,1);
for d = 1 : nFiles
    if isempty(selection)
        SEL{d} = 1 : size(DTItracts(d).fibindex,1);
    else
        SEL{d} = selection;
    end

    if ~isempty(ext_threshold)
        SEL{d} = find(DTItracts(d).pct_ext(SEL{d}) < ext_threshold);
    end
end
%% POINT_DATA
% Make arrays with only the selected datapoints
points = [];
t=0;
fprintf('------------------------------------------------\n')
fprintf('Fibre type written to file: %s\n',char(towrite))
for d = 1 : nFiles
    P = DTItracts(d).PolyCoeff;
    nTracts = length(SEL{d});

    ns = 10; % number of points to sample along the polynomial curve
    for j = 1 : nTracts
        i = SEL{d}(j);
        switch char(towrite)
            case 'raw'
                sel = min(DTItracts(d).fibindex(i,:)):max(DTItracts(d).fibindex(i,:));
                newpoints = DTItracts(d).tracts_xyz(:,sel);
            case 'trunc'
                sel = min(DTItracts(d).fibindex_trunc(i,:)):max(DTItracts(d).fibindex_trunc(i,:));
                newpoints = DTItracts(d).tracts_xyz(:,sel);
            case 'poly'
                % (get point_data from polynomials + extensions)
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
fprintf(fid, 'Tract data from MATLAB\n');
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

counter = 0;
for c = 1 : length(color_list)
    name = color_list{c};    
    if contains(name,{'pennation1','pennation2'})
        if ~isfield(DTItracts(1),'penangle')
            fprintf('%s is not a field in DTItracts. Color data not written.\n',name)
            continue
        end
    else
        if ~isfield(DTItracts(1),name)
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
        end
        V = [V;value_per_fibre];
    end
    
    color_data = zeros(nPoints,1);
    for j = 1:nTracts
        % Specify the same color for all elements in this fibre
        color_data(segm_fibre_nr == j) = V(j);
    end
    nScalars = 1;
    scalar_value = color_data;
    spec  = [repmat(['%0.',precision,'f '],1,nScalars) ' \n'];
    if counter == 0
        % Only write 'POINT_DATA' once
        fprintf(fid, 'POINT_DATA %d\n',nPoints);
        first = false;
    end
%     fprintf(fid, 'SCALARS %s float %d\n',name,nScalars); %ASCII header
    fprintf(fid, 'SCALARS %s float\n',name); %ASCII header
    fprintf(fid, 'LOOKUP_TABLE default\n'); %ASCII header
    fprintf(fid,spec,scalar_value);
    fprintf(' completed.\n')
    counter = counter + 1;
end
% Close the vtk file
fclose(fid);
fprintf('\nFibres succesfully written to %s\n',vtk_filename)
fprintf('File contains %d fibres, %d points and %d scalar indices\n',nTracts,nPoints,counter)
fprintf('------------------------------------------------\n')

% Open the file in ParaView, if ParaView is installed and added to the path
if ParaView == true
    try
        cmdtxt = sprintf('paraview --data=%s &',vtk_filename);
        system(cmdtxt)
    catch ME
        error(ME.message)

    end
end

end % of function
