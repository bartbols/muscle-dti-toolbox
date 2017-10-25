function CompareTracts(varargin)
%%COMPARETRACTS plots histograms of the parameters in DTItracts, which is a
% n x 1 structure of DTItracts, or a cell structure with DTI filenames.
%
% Bart Bolsterlee, Neuroscience Research Australia (NeuRA)
% October 2017
%
% ----------------- USAGE -----------------
% handles = CompareTracts  (user is requested to select tract files from a
%                           dialog box)
% handles = CompareTracts(TractFilename,names,'ParameterName',<value>)
%
% ----------------- INPUT -----------------
% input1     - n x 1 structure containing data from one or multiple sets of
%             tracts. A histogram curve will be plotted for all parameters
%             present in the structure, with a separate line for each set
%             of tracts. Alternatively, a cell array with tract filenames
%             can be provided.
% input2     - n x 1 cell structure with character arrays indicating the
%             name of each set of tracts. If not provided, the tracts will
%             be called 'tract1','tract2' etc. or (if filenames are
%             provided in DTItracts) the filenames.
%
% Optional parameters, provided as 'ParameterName',<value> pairs:
% - pct_threshold : Include only fibres that were extrapolated by less than
%                   this percentage of their total length.  Default = [] (include all fibres)
% - mm_threshold  : Include only fibres that were extrapolated by less than
%                   this length (in mm). Default = [] (include all fibres)


% Read input arguments
%% Check input
p = inputParser;
% addRequired(p,'DTItracts',@isstruct)
addOptional(p,'input1',[],@(x) isstruct(x) || iscell(x))
addOptional(p,'input2',[],@(x) iscell(x))
addParameter(p,'pct_threshold',Inf,@(x) isscalar(x) || isempty(x))
addParameter(p,'mm_threshold',Inf,@(x) isscalar(x) || isempty(x))
parse(p,varargin{:})

mm_threshold  = p.Results.mm_threshold;
pct_threshold = p.Results.pct_threshold;

if isempty(mm_threshold);mm_threshold = Inf;end
if isempty(pct_threshold);pct_threshold = Inf;end

input1 = p.Results.input1;
names  = p.Results.input2;
if isempty(input1)
    % No input arguments are provided. Select files from a dialog box.
    [filename, pathname] = uigetfile('*.mat','Select multiple DTI tract files','MultiSelect','on');
    if ischar(filename)
        % Only one file is selected. Load this file.
        filenames{1} = fullfile(pathname,filename);
        names{1} = filename;
    else
        % Multiple files are selected. Load all files.
        for i = 1 : length(filename)
            filenames{i} = fullfile(pathname,filename{i});
            names{i} = filename{i};
        end
    end
elseif ischar(input1)
    % one filename is provided
    filenames{1} = input1;
elseif iscell(input1)
    % A list of filenames is provided
    filenames = input1;
end

%% Load the data
if ~isstruct(input1)
    N = length(filenames);
    for i = 1 : N
        DTItracts(i) = load(filenames{i});
    end
else
    DTItracts = input1;
    N = length(DTItracts);
end

if isempty(names)
    % Give the tracts the default names
    for i = 1 : length(filenames)
        names{i} = sprintf('tract%02d',i);
    end
elseif ischar(names)
    % One name is provided. Convert to cell.
    names{1} = input2;
end

%%
figure('Name','Compare tracts')
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
colors = linspecer(N,'qualitative');
for i = 1 : N
    DTItracts(i).abs_ext = nansum(DTItracts(i).ext,2);
    CELLDATA = struct2cell(DTItracts(i));
    varNames = fieldnames(DTItracts(i));
    %     subplotnr = 1;
    for k = 1 : size(VARS,1)
        % Decide which subplot to display the data in
        %         subplotnr = subplotnr + 1;
        %         if mod(subplotnr-1,4) == 0;subplotnr = subplotnr + 1;end
        %         subplot(3,4,subplotnr)
        subplot(3,4,k);hold on
        %         set(subhandles(k),'Visible','on')
        
        % Check if variable exist. If not, continue with next.
        idx = strcmp(varNames,VARS{k,4});
        if ~any(idx)
            title([VARS{k,1} ' data not available'])
            %             subplotnr = subplotnr - 1;
            continue
        end
        
        % Get the data
        data = mean(CELLDATA{idx},2);
        
        % Select fibres within the defined thresholds for extrapolation.
        selection = DTItracts(i).abs_ext < mm_threshold & DTItracts(i).pct_ext < pct_threshold;
        
        % Plot the distribution
        h = MyHist(data(selection),VARS{k,3},VARS{k,1},VARS{k,2},...
            'Type','line',...
            'Color',colors(i,:),...
            'Text','mean',...
            'TextPos',[0.95 0.95-0.1*(i-1)],...
            'BarAlpha',0.5);
        if k == 1
            handle(i) = h;
        end
        axis tight
        
    end
end
% Make all y-axis equal
YLIM = cell2mat(get(findobj(gcf,'type','Axes'),'YLim'));
set(findobj(gcf,'type','Axes'),'YLim',[0 max(YLIM(:,2))])
legend(handle,names,'Interpreter','None')
