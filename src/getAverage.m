function [mean_value,median_value,SD_value,vars,selection] = getAverage( DTItracts,varargin )
%GETAVERAGE Calculates the average value of the variables in DTItracts.
%Several options can be provided to select fibres.
%
%
% Bart Bolsterlee, Neuroscience Research Australia (NeuRA)
% February 2017
%
% ----------------- USAGE ----------------- 
% [mean_value,median_value,SD_value,vars,selection] = getAverage(DTItracts,varargin)
% 
% ----------------- INPUT ----------------- 
% - DTItracts : a structure array (or a MAT-filename) containing the DTI tracts.
%
% Optional inputs, provided as 'parameter',<value> pairs:


p = inputParser;
addRequired(p,'DTItracts')
addParameter(p,'selection',[],@(x) isnumeric(x)||iscell(x))
addParameter(p,'pct_threshold',Inf,@(x) isscalar(x) || isempty(x))
addParameter(p,'mm_threshold',Inf,@(x) isscalar(x) || isempty(x))
addParameter(p,'min_length',0,@(x) isscalar(x) || isempty(x))
addParameter(p,'max_length',Inf,@(x) isscalar(x) || isempty(x))
addParameter(p,'min_pennation',0,@(x) isscalar(x) || isempty(x))
addParameter(p,'max_pennation',90,@(x) isscalar(x) || isempty(x))
addParameter(p,'min_curv',0,@(x) isscalar(x) || isempty(x))
addParameter(p,'max_curv',Inf,@(x) isscalar(x) || isempty(x))
addParameter(p,'min_ang',0,@(x) isscalar(x) || isempty(x))
addParameter(p,'max_ang',180,@(x) isscalar(x) || isempty(x))
addParameter(p,'max_z',Inf,@(x) isscalar(x) || isempty(x))
addParameter(p,'min_z',-Inf,@(x) isscalar(x) || isempty(x))
addParameter(p,'apo_to_mus',true,@(x) islogical(x) || x==0 || x==1)
addParameter(p,'vars',[],@(x) ischar(x) || iscell(x))
parse(p,DTItracts,varargin{:})

pct_threshold = p.Results.pct_threshold;
mm_threshold  = p.Results.mm_threshold;
min_length    = p.Results.min_length;
max_length    = p.Results.max_length;
min_pennation = p.Results.min_pennation;
max_pennation = p.Results.max_pennation;
min_curv      = p.Results.min_curv;
max_curv      = p.Results.max_curv;
min_ang       = p.Results.min_ang;
max_ang       = p.Results.max_ang;
vars          = p.Results.vars;
selection     = p.Results.selection;
apo_to_mus    = p.Results.apo_to_mus;
max_z         = p.Results.max_z;
min_z         = p.Results.min_z;

% Load data if a filename is provided instead of a structure array with the
% DTItract data
if ~isstruct(DTItracts)
    if exist(DTItracts,'file') == 2
        DTItracts = load(DTItracts);
    else
        error('Tract file %s does not exist. Check if extension was included.',DTItracts)
    end
end

% Check if a custom selection was provided
N = size(DTItracts.fibindex,1);
if isempty(selection)
    % Use all fibres
    selection = 1 : N;
end

if isfield(DTItracts,'pct_ext')
    % No custom selection was provided. Check if inclusion thresholds were
    % provided.
    if isempty(pct_threshold);pct_threshold = +Inf;end
    if isempty(mm_threshold); mm_threshold  = +Inf;end
    if isempty(min_length);   min_length    = 0;end
    if isempty(max_length);   max_length    = +Inf;end
    if isempty(min_curv);     min_curv      = 0;end
    if isempty(max_curv);     max_curv      = +Inf;end
    if isempty(min_ang);      min_ang       = 0;end
    if isempty(max_ang);      max_ang       = 180;end
    if isempty(min_pennation); min_pennation  = 0;end
    if isempty(max_pennation); max_pennation  = 90;end
    
    % The angle between endpoint slopes is usually not saved, so needs to
    % be calculated now.
    if ~isfield(DTItracts,'ang') && isfield(DTItracts,'endpoints_dir')
        DTItracts.ang = acosd(sum(squeeze(DTItracts.endpoints_dir(:,1,:)) .*...
                                  squeeze(DTItracts.endpoints_dir(:,2,:)),2));
    end
    % Add pennation angle as single-column vectors (mean pennation and
    % pennation at endpoint 1 and 2).
    if isfield(DTItracts,'penangle')
        DTItracts.pennation  = nanmean(DTItracts.penangle,2);
        DTItracts.pennation1 = DTItracts.penangle(:,1);
        DTItracts.pennation2 = DTItracts.penangle(:,2);
    end
    % Add absolute extension (sum of extensions at either end).
    if isfield(DTItracts,'ext')
        DTItracts.abs_ext =  nansum(DTItracts.ext,2);
    end
    
    % Select all fibres between the thresholds
    is_within_tresholds = DTItracts.pct_ext(selection) < pct_threshold & ...
        DTItracts.abs_ext(selection)       <= mm_threshold & ...
        DTItracts.fibrelength(selection)   >= min_length & ...
        DTItracts.fibrelength(selection)   <= max_length & ...
        DTItracts.curvature(selection)     >= min_curv & ...
        DTItracts.curvature(selection)     <= max_curv & ...
        DTItracts.ang(selection)           >= min_ang & ...
        DTItracts.ang(selection)           <= max_ang & ...
        DTItracts.pennation(selection)     >= min_pennation & ...
        DTItracts.pennation(selection)     <= max_pennation & ...
        all(DTItracts.endpoints(selection,:,3)     <= max_z,2) & ...
        all(DTItracts.endpoints(selection,:,3)     >= min_z,2);
    
    % Only include fibres that run from the aponeurosis to the muscle surface
    is_from_apo_to_mus = true(size(is_within_tresholds));
    if apo_to_mus == true && any(DTItracts.attach_type(:) == 2)
        is_from_apo_to_mus = DTItracts.attach_type(selection,1) ~= DTItracts.attach_type(selection,2);
    end
    % Remove the excluded fibres from the list of selected fibres.
    selection = selection(is_within_tresholds & is_from_apo_to_mus);
end

if isempty(vars)
    % If no custom selection of variables is provided, use the following
    % default values.
    vars = {'N_incl','pct_incl','abs_ext','pct_ext','length_mm',...
        'fibrelength','pennation','pennation1','pennation2','curvature',...
        'fa','md','lambda1','lambda2','lambda3'};
else
    if ischar(vars);vars = {vars};end
end


C = struct2cell(DTItracts);
F = fieldnames(DTItracts);

mean_value   = NaN(1,length(vars));
median_value = NaN(1,length(vars));
SD_value     = NaN(1,length(vars));
for i = 1 : length(vars)
    
    switch vars{i}
        case 'N_incl'
            data = length(selection);
        case 'pct_incl'
            data = length(selection) / N * 100;
            
        otherwise
            % Check if variable exists
            idx = find(strcmp(vars{i},F));
            if isscalar(idx)
                data = C{idx}(selection);
            else
                % Data is not available for this variable
                data = NaN;
            end
            
    end
    mean_value(i)   = nanmean(data);
    median_value(i) = nanmedian(data);
    SD_value(i)     = nanstd(data);
end


end % of function





