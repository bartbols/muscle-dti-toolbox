function [incl,settings,pre_selected] = include_fibres( DTItracts,varargin )
%EXCLUDE_FIBRES This function selects fibres for inclusion based on
%quantitative comparison of architecture parameters with neighbouring
%tracts. For each fibre, fibres within a selected radius (default = 2 mm)
%around the midpoint of the fibre are selected. If the length of the fibre
%is in the lowest or highest 10th percentile (default value) of the lengths
%of all neighbouring fibres, the fibres is excluded. The fibre is also
%marked for exclusion if the fibre length deviates by more than 5 mm
%(default) of the median length of its neighbours.
%
% Bart Bolsterlee, Neuroscience Research Australia
% March 2018
%
% USAGE:
% idx = exclude_fibres(DTItracts,varargin)
%
%
% INPUT:
% Required inputs
% - DTItracts: MATLAB structure (or filename) containing the DTI tracts and
%              architectural parameters (fibrelength)

% Optional inputs, provided as 'parameter',<value> pairs:
% - selection : list of fibre indices used for the quantitative analysis.
% - r         : radius (in mm) around the fibre midpoint to mark fibres as
%               neighbours. Fibres with at least one tract point within
%               distance r of the midpoint of the selected fibre are
%               considered neighbours. Default = 2.5 mm
% pct_threshold: percentage extension above which fibres are excluded.
%                Default: 30
% mm_threshold : absolute extension (in mm) above which fibres are
%                excluded. Default: 20
% min_length   : minimum fibre length. Default is minimum length used for
%                fibre tracking.
% max_length   : maximum fibre length. Default is maximum length used for
%                fibre tracking,
% apo_to_mus   : logical (true/false) deciding whether only fibres that run
%                from the aponeurosis to the surface model are included.
%                Default is true if different attachment types are present
%                in the variable 'attach_type'.
% prctile      : scalar indicating the lower and upper percentile of lengths
%                to exclude fibres. For example, if prctile = 10, a fibre
%                is excluded if it falls in the bottom or top 10 percent of
%                the lengths of all its neighbours. Default = 10.
% max_diff     : scalar indicating the maximum difference between the
%                median fibre length of all neighbouring fibres and the
%                selected fibre to be included.
%
% OUTPUT:
% excl        : list of fibre indices excluded based on the selected exclusion
%               criteria.

p = inputParser;
addRequired(p,'DTItracts')
addParameter(p,'selection',[],@isnumeric)
addParameter(p,'r',2.5,@(x) isscalar(x) || isempty(x))
addParameter(p,'pct_threshold',30,@(x) isscalar(x) || isempty(x))
addParameter(p,'mm_threshold',20,@(x) isscalar(x) || isempty(x))
addParameter(p,'min_length',[],@(x) isscalar(x) || isempty(x))
addParameter(p,'max_length',[],@(x) isscalar(x) || isempty(x))
addParameter(p,'max_diff',7.5,@(x) isscalar(x) || x > 0)
addParameter(p,'apo_to_mus',true,@(x) islogical(x) || x==0 || x==1)
addParameter(p,'diagnostics',false,@(x) islogical(x) || x==0 || x==1)
addParameter(p,'prctile',5,@(x) isscalar(x) || isempty(x))
addParameter(p,'constraint',true,@(x) islogical(x) || x==0 || x==1)

% Read inputs
parse(p,DTItracts,varargin{:});
selection     = p.Results.selection;
settings.pct_threshold = p.Results.pct_threshold;
settings.mm_threshold  = p.Results.mm_threshold;
settings.min_length    = p.Results.min_length;
settings.max_length    = p.Results.max_length;
settings.r             = p.Results.r;
settings.prct          = p.Results.prctile;
settings.max_diff      = p.Results.max_diff;
settings.apo_to_mus    = p.Results.apo_to_mus;
diagnostics            = p.Results.diagnostics;
constraint             = p.Results.constraint;

tic

% Check if filename or structure with DTItracts is provided.
if ~isstruct(DTItracts)
    DTItracts = load(DTItracts);
end
% By default, use the minimum and maximum length from the fibre track
% settings.
if isempty(settings.min_length)
    settings.min_length = DTItracts.TrackSettings.MinLength;
end
if isempty(settings.max_length)
    settings.max_length = DTItracts.TrackSettings.MaxLength;
end

N = length(DTItracts.length_mm);
% Report selection criteria
fprintf('Inclusion settings:\n')
formatstr = '%-35s: %.2f %s\n';
fprintf(formatstr,'Percentage ext. threshold',settings.pct_threshold,'%')
fprintf(formatstr,'Absolute ext. threshold',settings.mm_threshold,'mm')
fprintf(formatstr,'Minimum length',settings.min_length,'mm')
fprintf(formatstr,'Maximum length',settings.max_length,'mm')
fprintf(formatstr,'radius',settings.r,'mm')
fprintf(formatstr,'Exclusion percentile',settings.prct,'')
fprintf(formatstr,'Maximum difference from median',settings.max_diff,'mm')

if constraint == true
    fibre_length = DTItracts.fibrelength;
    
    % First, remove fibres that are extrapolated by more than the thresholds
    % and don't have endpoints at both the muscle surface and the aponeurosis.
    is_within_tresholds = DTItracts.pct_ext < settings.pct_threshold & ...
        nansum(DTItracts.ext,2) < settings.mm_threshold & ...
        fibre_length >= settings.min_length & ...
        fibre_length <= settings.max_length;
    
    if isempty(selection)
        selection = 1 : length(is_within_tresholds);
    end
    
    is_selected = false(size(is_within_tresholds));
    is_selected(selection) = true;
    
    % Only include fibres that run from the aponeurosis to the muscle surface
    if settings.apo_to_mus == true
        if any(DTItracts.attach_type(:) == 2)
            is_from_apo_to_mus = DTItracts.attach_type(:,1) ~= DTItracts.attach_type(:,2);
        end
    else
        is_from_apo_to_mus = true(size(is_within_tresholds));
    end
    
    % Include only those fibres that are within the thresholds, are within the
    % user-provided selection (if provided at all) and run from the aponeurosis
    % to the muscle surface.
    pre_selected = find(is_within_tresholds & is_selected & is_from_apo_to_mus);    
    fibindex = DTItracts.fibindex_trunc;
    
    % Calculate midpoints of fibre tracts. Exclude NaNs.
    MP = zeros(3,N);
    tmp = round(mean(fibindex,2));
    MP(:,~isnan(tmp)) = DTItracts.tracts_xyz(:,~isnan(tmp));
else
    % Use tract length (instead of fascicle length) and don't use
    % pre-selection based on extension thresholds.
    fibre_length = DTItracts.length_mm;
    pre_selected = 1 : size(DTItracts.fibindex,1);
    fibindex = DTItracts.fibindex;
    % Calculate midpoints of fibre tracts
    MP = DTItracts.tracts_xyz(:,round(mean(fibindex,2)));
end

incl = false(1,N);
incl(pre_selected) = true;

nFib = length(pre_selected);
hwait = waitbar(0,sprintf('Selecting %d fibres',nFib),...
    'Name','Progress bar include_fibres');

for c = 1 : nFib
    
    if any(round(linspace(1,nFib,50))==c)
        waitbar(c/nFib,hwait)
    end
    
    fibnr = pre_selected(c);
%     % Find midpoint of polynomial curve
%     P = [polyval(DTItracts.PolyCoeff(fibnr).x,(DTItracts.PolyCoeff(fibnr).t0 + DTItracts.PolyCoeff(fibnr).t1)/2),...
%         polyval(DTItracts.PolyCoeff(fibnr).y,(DTItracts.PolyCoeff(fibnr).t0 + DTItracts.PolyCoeff(fibnr).t1)/2),...
%         polyval(DTItracts.PolyCoeff(fibnr).z,(DTItracts.PolyCoeff(fibnr).t0 + DTItracts.PolyCoeff(fibnr).t1)/2)]';
    
    % Find fibres with at least one point within a radius r around P
    idx = find_fibres(DTItracts.tracts_xyz,fibindex,MP(:,fibnr),settings.r);
    
    % Remove the fibres that were excluded based on length and extension
    % criteria.
    idx = idx(ismember(idx,pre_selected));
    
    if isempty(idx)
        % No neighbouring fibres found. Remove fibre and continue.
        incl(fibnr) = false;
        continue
    end
    
    % Calculate fibre length thresholds as percentiles
    threshold = prctile(fibre_length(idx),[settings.prct 100-settings.prct]);
    
    % Calculate the difference in length between the current fibre tract
    % and the median of all neighbouring tracts
    diff_from_median = abs(fibre_length(fibnr) - median(fibre_length(idx)));
    if fibre_length(fibnr) < threshold(1) || ...
            fibre_length(fibnr) > threshold(2) || ...
            diff_from_median > settings.max_diff
        
        % Fibre is excluded
        incl(fibnr) = false;
        txt = 'no';
        
    else
        % Fibre is included
        incl(fibnr) = true;
        txt = 'yes';
    end
    if incl(fibnr) == false && diagnostics == true
        fig = figure;
        subplot(1,2,1)
        hold on
        plot3(P(1),P(2),P(3),'o','MarkerSize',20,'MarkerFaceColor','k',...
            'MarkerEdgeColor','none')
        h1 = InspectTracts('Tracts',DTItracts,'PlotStats',false,...
            'selection',fibnr,'ToPlot','poly','color','b','linewidth',5);
        h2 = InspectTracts('Tracts',DTItracts,'PlotStats',false,...
            'selection',idx,'ToPlot','poly','color','r','linewidth',0.5);
        %     plot_surface(filename.apo_model,'FaceColor','none','FaceAlpha',1,'EdgeColor','k','EdgeAlpha',0.3)
        subplot(1,2,2)
        MyHist(fibre_length(idx),0.2,'Fascicle length','mm',...
            'Text','all')
        text(0.95,0.50,{sprintf('Fascicle length = %.2f\nIncluded: %s\nThreshold = %.2f to %.2f\nDiff. from median = %.2f',...
            fibre_length(fibnr),txt,threshold(1),threshold(2),diff_from_median)},...
            'Units','Normalized',...
            'FontWeight','normal',...
            'FontSize',10,...
            'HorizontalAlignment','right',...
            'VerticalAlignment','bottom')
        pause
        close(fig)
    end
    
end

incl = find(incl);
t_elapsed = toc;
fprintf('%d out of %d (%.2f%%) pre-selected fibres were included.\n',...
    length(incl),length(pre_selected),length(incl) / length(pre_selected) * 100)
fprintf('It took %.2f seconds to include fibres.\n',t_elapsed)
close(hwait)


end % of function

