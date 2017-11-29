function varargout = InspectTracts(varargin)
%%INSPECTTRACTS This function can be used to inspect the results of fiber
%tracking and subsequent tract analyses.
%
% Bart Bolsterlee, Neuroscience Research Australia (NeuRA)
% February 2017
%
% ----------------- USAGE -----------------
% handles = InspectTracts  (user is requested to select a tract file from a
%                        dialog box)
% handles = InspectTracts(TractFilename,'ParameterName',<value>)
%
% ----------------- INPUT -----------------
%
% This function has only optional inputs, provided as 'ParameterName',<value> pairs:
%
% - Tracts:    filename of the tract file (including extension .mat)
%              OR structure array DTItracts
%
% - ToPlot: cell string with data to plot. The options are:
%            - 'raw'   : raw tracts
%            - 'trunc' : truncated tracts
%            - 'cut'   : the parts of the tracts that were truncated
%            - 'poly'  : the polynomial fitted fibres
%       e.g. 'ToPlot', {'raw','poly'} will plot the raw and the polynomial fitted
%       fibres.
%       Default: {'raw','poly'}
% - Color: list of colors for the tracts, in the order of ToPlot. Default:
%          'jet' colormap
% - LineWidth: linewidth for tracts/fibres. Default = 1
% - SurfModel: structure array containing fields 'vertices' and 'faces' OR
%              the full filename of the STL-file of the surface model.
%              If provided, the surface model will be added to the plot.
% - aponeurosis: structure array containing fields 'vertices' and 'faces' OR
%              the full filename of the STL-file of the aponeurosis model.
%              If provided, the aponeurosis model will be added to the plot.
% - SurfModelColor: Color of the surface model. Default: 'y'
% - SurfModelAlpha: alpha level of the surface model patch. Default: 0.2
% - Selection     : vector of indices of fibres to plot. If not 
%                   provided, all fibres are displayed. Both row and column
%                   vectors are accepted.
% - PlotStats     : also plot the histograms of the parameters in separate
%                   subplots. Default: true
% - fraction      : fraction (0 < fraction <= 1) of total tracts to plot. Default: 1 (=100%)
% - pct_threshold : Include only fibres that were extrapolated by less than
%                   this percentage of their total length.  Default = [] (include all fibres)
% - mm_threshold  : Include only fibres that were extrapolated by less than
%                   this length (in mm). Default = [] (include all fibres)
%
% ----------------- OUTPUT -----------------
% handles : handles to the plot objects
%
% Examples:
%
% handles = InspectTracts('Tracts',DTItracts);
%
% handles = InspectTracts('Tracts',tract_filename,'SurfModel,surf_filename);
%
% handles = InspectTracts('Tracts',DTItracts,'ToPlot',{'raw','poly','cut'},...
%                'Color',{'r','g','c'},'SurfModel',surf_model,...
%                'SurfModelColor','r','SurfModelAlpha',0.3)
%
% InspectTracts('Tracts',DTItracts,'mm_threshold',30)

%% Check input
p = inputParser;
addParameter(p,'Tracts',[],@(x) isstruct(x) || exist(x,'file')==2)
addParameter(p,'ToPlot',{'raw','poly'},@(x) iscell(x) || ischar(x))

% Surface model
addParameter(p,'SurfModel',[],@(x) isstruct(x) || endsWith(lower(x),'.stl'))
addParameter(p,'SurfModelColor','y',@(x) ischar(x) || isnumeric(x))
addParameter(p,'SurfModelAlpha',0.25,@(x)validateattributes(x,{'numeric'},{'scalar'}))
addParameter(p,'SurfModelEdgeColor','y',@(x) ischar(x) || isnumeric(x))
addParameter(p,'SurfModelEdgeAlpha',0.1,@(x)validateattributes(x,{'numeric'},{'scalar'}))

% Aponeurosis model
addParameter(p,'aponeurosis',[],@(x) isstruct(x) || endsWith(x,'.stl','IgnoreCase',true))
addParameter(p,'aponeurosisColor','r',@(x) ischar(x) || isnumeric(x))
addParameter(p,'aponeurosisAlpha',0.5,@(x)validateattributes(x,{'numeric'},{'scalar'}))
addParameter(p,'aponeurosisEdgeColor','r',@(x) ischar(x) || isnumeric(x))
addParameter(p,'aponeurosisEdgeAlpha',0.25,@(x)validateattributes(x,{'numeric'},{'scalar'}))

% Some more plot options
addParameter(p,'Color',jet(4),@(x) iscell(x) || isnumeric(x) || ischar(x))
addParameter(p,'LineWidth',1,@(x)validateattributes(x,{'numeric'},{'scalar'}))
addParameter(p,'PlotStats',true,@(x) x==0 || x==1 || islogical(x) )
addParameter(p,'Selection',[],@(x) isnumeric(x))
addParameter(p,'fraction',1,@(x) validateattributes(x,{'numeric'},{'scalar','>',0,'<=',1}))
addParameter(p,'pct_threshold',[],@(x) isscalar(x) || isempty(x))
addParameter(p,'mm_threshold',[],@(x) isscalar(x) || isempty(x))
parse(p,varargin{:})

ToPlot          = p.Results.ToPlot;
Color           = p.Results.Color;
SurfModel       = p.Results.SurfModel;
LineWidth       = p.Results.LineWidth;
PlotStats       = p.Results.PlotStats;
Tracts          = p.Results.Tracts;
Selection       = p.Results.Selection;
fraction        = p.Results.fraction;
mm_threshold    = p.Results.mm_threshold;
pct_threshold   = p.Results.pct_threshold;
aponeurosis     = p.Results.aponeurosis;

SurfModelColor       = p.Results.SurfModelColor;
SurfModelAlpha       = p.Results.SurfModelAlpha;
SurfModelEdgeColor   = p.Results.SurfModelEdgeColor;
SurfModelEdgeAlpha   = p.Results.SurfModelEdgeAlpha;

aponeurosisColor     = p.Results.aponeurosisColor;
aponeurosisAlpha     = p.Results.aponeurosisAlpha;
aponeurosisEdgeColor = p.Results.aponeurosisEdgeColor;
aponeurosisEdgeAlpha = p.Results.aponeurosisEdgeAlpha;

if ischar(ToPlot);ToPlot={ToPlot};end
if ischar(Color);Color={Color};end

% If input1 is not provided, select a file from a dialog box.
if isempty(Tracts)
    [filename,pathstr]= uigetfile('*.mat','Select a DTI tract file');
    fullfilename = fullfile(pathstr,filename);
    Tracts = fullfilename;
end

% Check if a filename or a structure array is provided as input
if isstruct(Tracts)
    DTItracts = Tracts;
    clear Tracts
else
    if exist(Tracts,'file') == 2
        DTItracts      = load(Tracts);
        fullfilename   = Tracts;
        clear Tracts
    else
        error('%s not found.',Tracts)
    end
end

if isempty(Selection)
    % If no selection is provided, plot all fibres.
    Selection = 1:size(DTItracts.fibindex,1);
else
    % If column vector is provided, transpose to row vector.
    [nr,nc] = size(Selection);
    if nc == 1 && nr > 1
        Selection = Selection';
    end
    if islogical(Selection)
        Selection = find(Selection);
    end
end

% Check if extrapolation threshold are provided. If so, exclude fibres 
% above these thresholds
if ~isempty(mm_threshold)
    mm_ext = sum(DTItracts.ext(Selection,:),2);
    Selection(mm_ext > mm_threshold) = [];
end
if ~isempty(pct_threshold)
    Selection(DTItracts.pct_ext(Selection) > pct_threshold) = [];
end

% Check if surface model with fields 'vertices' and 'faces' is provided or
% the filename of the STL file. If an STL-filename is provided, load the
% surface model.
if ~isempty(SurfModel)
    if ~isstruct(SurfModel)
        if exist(SurfModel,'file') == 2
            SurfModel = stlread(SurfModel);
        else
            warning('%s does not exist. Surface model not displayed.',SurfModel)
            SurfModel = [];
        end
    end
end

% Check if aponeurosis model with fields 'vertices' and 'faces' is provided or
% the filename of the STL file. If an STL-filename is provided, load the
% surface model.
if ~isempty(aponeurosis)
    if ~isstruct(aponeurosis)
        if exist(aponeurosis,'file') == 2
            aponeurosis = stlread(aponeurosis);
        else
            warning('%s does not exist. Aponeurosis model not displayed.',SurfModel)
            aponeurosis = [];
        end
    end
end
%%

if PlotStats == true
    if exist('fullfilename','var')==1
        FigureTitle = fullfilename;
    else
%         FigureTitle = DTItracts.FileNames.Tracts;
        FigureTitle = '';
    end
    fig = figure('Name',sprintf('InspectTracts for %s',FigureTitle),'NumberTitle','Off');
    axis tight
    handle_3D = subplot(3,5,[1 6 11]);
    hold on
    
    subplotnr = 1;
    for ii = 1 : 12
        subplotnr = subplotnr + 1;
        if mod(subplotnr-1,5) == 0;subplotnr = subplotnr + 1;end
        subhandles(ii) = subplot(3,5,subplotnr);
        set(subhandles(ii),'Visible','off')
        hold on
    end
else
    handle_3D = gca;
end

% Get the hold state at the start of the function so it can be changed to
% this state at the end.
holdstate = ishold;
hold on

% Convert 'option' to cell array if necessary
if ischar(ToPlot);ToPlot = {ToPlot};end


% Reduce the tracts for plotting to the fraction of the total number 
% of tracts specified by 'fraction'.
if fraction == 1
    SelectionPlot = Selection;
else
    tmp = randperm(length(Selection),round(fraction*length(Selection)));
    SelectionPlot = Selection(tmp);
end
k=0;
for opt = ToPlot
    k = k + 1;
    
    switch char(opt)
        case 'raw'
            % Plot the raw tracts
            PlotX = [];PlotY =[];PlotZ = [];
            for fibnr = SelectionPlot
                PlotX = [PlotX DTItracts.tracts_xyz(1,min(DTItracts.fibindex(fibnr,1:2)):max(DTItracts.fibindex(fibnr,1:2))) NaN];
                PlotY = [PlotY DTItracts.tracts_xyz(2,min(DTItracts.fibindex(fibnr,1:2)):max(DTItracts.fibindex(fibnr,1:2))) NaN];
                PlotZ = [PlotZ DTItracts.tracts_xyz(3,min(DTItracts.fibindex(fibnr,1:2)):max(DTItracts.fibindex(fibnr,1:2))) NaN];
            end
            legendTxt{k} = 'Raw tracts';
            
        case 'trunc'
            % Plot the truncated tracts
            %             nFib = size(DTItracts.fibindex_trunc,1);
            PlotX = [];PlotY =[];PlotZ = [];
            for fibnr = SelectionPlot
                if any(isnan(DTItracts.fibindex_trunc(fibnr,1:2)))
                    continue
                end
                PlotX = [PlotX DTItracts.tracts_xyz(1,min(DTItracts.fibindex_trunc(fibnr,1:2)):max(DTItracts.fibindex_trunc(fibnr,1:2))) NaN];
                PlotY = [PlotY DTItracts.tracts_xyz(2,min(DTItracts.fibindex_trunc(fibnr,1:2)):max(DTItracts.fibindex_trunc(fibnr,1:2))) NaN];
                PlotZ = [PlotZ DTItracts.tracts_xyz(3,min(DTItracts.fibindex_trunc(fibnr,1:2)):max(DTItracts.fibindex_trunc(fibnr,1:2))) NaN];
                
            end
            legendTxt{k} = 'Truncated tracts';
        case 'cut'
            % Plot only the parts of the fascicles that were cut off
            nFib = size(DTItracts.fibindex_trunc,1);
            PlotX = [];PlotY =[];PlotZ = [];
            for fibnr = SelectionPlot
                if any(isnan(DTItracts.fibindex_trunc(fibnr,1:2)))
                    continue
                end

                PlotX = [PlotX DTItracts.tracts_xyz(1,min(DTItracts.fibindex(fibnr,1:2)):(min(DTItracts.fibindex_trunc(fibnr,1:2))-1)) NaN,...
                    DTItracts.tracts_xyz(1,(max(DTItracts.fibindex_trunc(fibnr,1:2))+1):max(DTItracts.fibindex(fibnr,1:2))) NaN];
                PlotY = [PlotY DTItracts.tracts_xyz(2,min(DTItracts.fibindex(fibnr,1:2)):(min(DTItracts.fibindex_trunc(fibnr,1:2))-1)) NaN,...
                    DTItracts.tracts_xyz(2,(max(DTItracts.fibindex_trunc(fibnr,1:2))+1):max(DTItracts.fibindex(fibnr,1:2))) NaN];
                PlotZ = [PlotZ DTItracts.tracts_xyz(3,min(DTItracts.fibindex(fibnr,1:2)):min(DTItracts.fibindex_trunc(fibnr,1:2))-1) NaN,...
                    DTItracts.tracts_xyz(3,(max(DTItracts.fibindex_trunc(fibnr,1:2))+1):max(DTItracts.fibindex(fibnr,1:2))) NaN];
            end
            legendTxt{k} = 'Cut parts of tracts';
        case 'poly'
            % Plot the polynomial fitted tracts, including extrapolations
            PlotX = [];PlotY = [];PlotZ = [];
            if isfield(DTItracts,'PolyCoeff')
                P = DTItracts.PolyCoeff;
            else
                error('''PolyCoeff'' not found as a field in DTITracts. Polynomials fits cannot be plotted.')
            end
            for fibnr = SelectionPlot
                if isempty(P(fibnr).x);continue;end
                if any(isnan(DTItracts.fibindex_trunc(fibnr,1:2)))
                    continue
                end

                tmpX = [DTItracts.endpoints(fibnr,1,1) polyval(P(fibnr).x,linspace(P(fibnr).t0,P(fibnr).t1,100)) DTItracts.endpoints(fibnr,2,1)];
                tmpY = [DTItracts.endpoints(fibnr,1,2) polyval(P(fibnr).y,linspace(P(fibnr).t0,P(fibnr).t1,100)) DTItracts.endpoints(fibnr,2,2)];
                tmpZ = [DTItracts.endpoints(fibnr,1,3) polyval(P(fibnr).z,linspace(P(fibnr).t0,P(fibnr).t1,100)) DTItracts.endpoints(fibnr,2,3)];
                PlotX = [PlotX NaN tmpX];
                PlotY = [PlotY NaN tmpY];
                PlotZ = [PlotZ NaN tmpZ];
            end
            legendTxt{k} = 'Polynomial fits+extrapolations';
            
        otherwise
            warning('unknown option ToPlot %s in PlotTracts.m',char(opt))
    end
    
    % Plot the tracts
    if iscell(Color)
        c = Color{k};
    else
        c = Color(k,:);
    end
    subplot(handle_3D);
    handles(k) = plot3(PlotX(:),PlotY(:),PlotZ(:),'LineWidth',LineWidth,'Color',c);
end

if ~isempty(SurfModel)
    
    % Plot the surface model
    handles(k+1) = patch('Vertices',SurfModel.vertices,...
        'Faces',SurfModel.faces,...
        'FaceColor',SurfModelColor,...
        'FaceAlpha',SurfModelAlpha,...
        'EdgeColor',SurfModelEdgeColor,...
        'EdgeAlpha',SurfModelEdgeAlpha);
    legendTxt{k+1} = 'Muscle surface';
%     set(gca,'XLim',[min(SurfModel.vertices(:,1)),max(SurfModel.vertices(:,1))],...
%         'YLim',[min(SurfModel.vertices(:,2)),max(SurfModel.vertices(:,2))],...
%         'ZLim',[min(SurfModel.vertices(:,3)),max(SurfModel.vertices(:,3))])
end

if ~isempty(aponeurosis)
    
    % Plot the surface model
    handles(k+2) = patch('Vertices',aponeurosis.vertices,...
        'Faces',aponeurosis.faces,...
        'FaceColor',aponeurosisColor,...
        'FaceAlpha',aponeurosisAlpha,...
        'EdgeColor',aponeurosisEdgeColor,...
        'EdgeAlpha',aponeurosisEdgeAlpha);
    legendTxt{k+2} = 'Aponeurosis';
end

if exist('FigureTitle','var')
    [~,titleTxt] = fileparts(FigureTitle);
else
    titleTxt = '';
end
title(titleTxt,'Interpreter','None')
legend(handles,legendTxt)
view(-37.5,30)
axis equal off
% if PlotStats == true
%     set(gca,'Position',[0.13 0.11 0.2 0.8])
% end

% Add some lights
lightangle(-37.5,30)
lightangle(-37.5+180,30)

% Set the 'hold' state back to the same value as at the start of the
% function.
if holdstate == 0
    hold off
end

%% Plot distribution of architecture and DTI indices for all tracts
if PlotStats == true
    VARS = {'Raw tract length','mm',2,'length_mm';...
        'Percentage extension','%',2,'pct_ext';...
        'Abs. extension','mm',2,'abs_ext';...
        'Fascicle length','mm',2,'fibrelength';...
        'Pennation angle','degr',2,'penangle';...
        'Curvature','1/m',2,'curvature';...
        'Fractional anisotropy','-',0.01,'fa';...
        'Mean diffusivity','1e-9 mm^2/s',0.05,'md';...
        '\lambda_1','1e-3 mm^2/s',0.05,'lambda1';...
        '\lambda_2','1e-3 mm^2/s',0.05,'lambda2';...
        '\lambda_3','1e-3 mm^2/s',0.05,'lambda3'};
    
    % Make variable abs_ext (absolute extension = sum of extension at
    % either end)
    if isfield(DTItracts,'ext')
        DTItracts.abs_ext = nansum(DTItracts.ext,2);
    end
    CELLDATA = struct2cell(DTItracts);
    varNames = fieldnames(DTItracts);
%     subplotnr = 1;
    for k = 1 : size(VARS,1)
        % Decide which subplot to display the data in
%         subplotnr = subplotnr + 1;
%         if mod(subplotnr-1,4) == 0;subplotnr = subplotnr + 1;end
%         subplot(3,4,subplotnr)
        subplot(subhandles(k))
        set(subhandles(k),'Visible','on')
        
        % Check if variable exist. If not, continue with next.
        idx = strcmp(varNames,VARS{k,4});
        if ~any(idx)
            title([VARS{k,1} ' data not available'])
            %             subplotnr = subplotnr - 1;
            continue
        end
        
        % Get the data
        data = mean(CELLDATA{idx},2);
        
        % Filter the data to only plot histograms for the selected fibres
        % (if not selection is provided, histograms for all fibres are
        % provided)
        data = data(Selection);
        
        % Plot a histogram
        MyHist( data,VARS{k,3},VARS{k,1},VARS{k,2});
        axis tight
    end
    
    set(gcf,'Position',get(0,'ScreenSize') + [0 100 0 -200])
end
if nargout == 1
    varargout{1} = handles;
end
rotate3d
