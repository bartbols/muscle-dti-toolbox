function [mean_value,median_value,SD_value,vars] = getAverage( DTItracts,varargin )
%GETAVERAGE Calculates the average value of the variables in DTItracts.

p = inputParser;
addRequired(p,'DTItracts')
addParameter(p,'selection',[],@(x) isnumeric(x)||iscell(x))
addParameter(p,'pct_threshold',[],@(x) isscalar(x) || isempty(x))
addParameter(p,'mm_threshold',[],@(x) isscalar(x) || isempty(x))
addParameter(p,'min_length',[],@(x) isscalar(x) || isempty(x))
addParameter(p,'max_length',[],@(x) isscalar(x) || isempty(x))
addParameter(p,'apo_to_mus',true,@(x) islogical(x) || x==0 || x==1)
addParameter(p,'vars',[],@(x) ischar(x) || iscell(x))
parse(p,DTItracts,varargin{:})
pct_threshold = p.Results.pct_threshold;
mm_threshold  = p.Results.mm_threshold;
min_length    = p.Results.min_length;
max_length    = p.Results.max_length;
vars          = p.Results.vars;
selection     = p.Results.selection;
apo_to_mus    = p.Results.apo_to_mus;

% Load data if a filename is provided instead of a structure array with the
% DTItract data
if ~isstruct(DTItracts)
    if exist(DTItracts,'file') == 2
        DTItracts = load(DTItracts);
    else
        error('Tract file %s does not exist. Check if extension was included.',DTItracts)
    end
end

N = size(DTItracts.fibindex,1);
if isempty(selection) && isfield(DTItracts,'pct_ext')
    % No custom selection was provided. Check if inclusion thresholds were
    % provided.
    if isempty(pct_threshold);pct_threshold = +Inf;end
    if isempty(mm_threshold); mm_threshold  = +Inf;end
    if isempty(min_length);   min_length    = -Inf;end
    if isempty(max_length);   max_length    = +Inf;end
    
    % First, remove fibres that are extrapolated by more than the thresholds
    % and don't have endpoints at both the muscle surface and the aponeurosis.
    is_within_tresholds = DTItracts.pct_ext < pct_threshold & ...
        nansum(DTItracts.ext,2) < mm_threshold & ...
        DTItracts.fibrelength >= min_length & ...
        DTItracts.fibrelength <= max_length;
    
    % Only include fibres that run from the aponeurosis to the muscle surface
    if apo_to_mus == true
        if any(DTItracts.attach_type(:) == 2)
            is_from_apo_to_mus = DTItracts.attach_type(:,1) ~= DTItracts.attach_type(:,2);
        end
    else
        is_from_apo_to_mus = true(size(is_within_tresholds));
    end
    selection = find(is_within_tresholds & is_from_apo_to_mus);
    
else
    % Use all fibres
    selection = 1 : N;
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

% Add pennation angle as single-column vectors (mean pennation and
% pennation at endpoint 1 and 2).
if isfield(DTItracts,'penangle')
    DTItracts.pennation = nansum(DTItracts.penangle,2);
    DTItracts.pennation1 = DTItracts.penangle(:,1);
    DTItracts.pennation2 = DTItracts.penangle(:,2);
end
% Add absolute extension (sum of extensions at either end).
if isfield(DTItracts,'ext')
    DTItracts.abs_ext =  nansum(DTItracts.ext,2);
end
C = struct2cell(DTItracts);
F = fieldnames(DTItracts);

mean_value   = NaN(1,length(vars));
median_value = NaN(1,length(vars));
SD_value     = NaN(1,length(vars));
for i = 1 : length(vars)
    
    switch vars{i}
        case 'N_incl'
            data = N;
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





