function varargout = TruncateTracts( DTItracts,surf_model)
%%TRUNCATETRACTS Truncates the fiber tracts so that they terminate inside
% the muscle surface model. The truncated length and truncated fibre
% indices are returned.
%
% Bart Bolsterlee, Neuroscience Research Australia (NeuRA)
% February 2017
%
% ----------------- USAGE -----------------
% [ fibindex_trunc,length_trunc ] = TruncateTracts( DTItracts,surface)
% 
% ----------------- INPUT -----------------
% - DTItracts : a structure array containing at least the fields tracts_xyz
%               and fibindex.
% - surf_model   : a structure array 'surface' containing the fields
%   'vertices','faces' with the vertices and faces of the surface used 
%   for truncating the tracts.
%         
% ----------------- OUTPUT -----------------
% - fibindex_trunc: n x 2 array with the first and last index of the
%    truncated fibre tracts that are within the surface
% - length_trunc: n x 1 array with the length of the truncated fibre tracts
%
% If one output argument is provided, the fields 'fibindex_trunc' and
% 'length_trunc' are added to DTItracts. If two outputs are provided, the
% first output is 'fibindex_trunc' and the second output 'length_trunc'
%
% Uses the function inside_surface.m to determine whether a point is inside
% or outside the surface.

%% Check inputs
if nargin ~= 2
    error('Wrong number of input arguments. Usage is TruncateTracts(DTItracts,surface)')
end

tic
p = inputParser;
addRequired(p,'DTItracts',@isstruct)
addRequired(p,'surf_model',@isstruct)
parse(p,DTItracts,surf_model)

% Check if all required fields are available in DTItracts and surface
if ~isfield(DTItracts,'tracts_xyz')
    error('Required field tracts_xyz not found in DTItracts.')
end
if ~isfield(DTItracts,'fibindex')
    error('Required field fibindex not found in DTItracts.')
end

if ~isfield(surf_model,'vertices')
    error('Required field vertices not found in surf_model.')
end

if ~isfield(surf_model,'faces')
    error('Required field faces not found in surf_model.')
end

%% Fibre truncation
nFib = size(DTItracts.fibindex,1);
fibindex_trunc = NaN(nFib,2);
length_trunc   = NaN(nFib,1);
hwait = waitbar(0,'Calculating which tract points are inside and outside surface model',...
    'Name','Progress bar TruncateTracts');

% Calculate whether point is inside or outside of surface model. It is much
% faster to do this for all points at once then inside the loop for each
% fibre individually.

% TSIG =  inpolyhedron(surface,DTItracts.tracts_xyz');
INSIDE =  inside_surface(surf_model,DTItracts.tracts_xyz');

for fibnr = 1:nFib
    waitbar(fibnr/nFib,hwait,sprintf('Truncating fibre %d of %d',fibnr,nFib))
    
    % Remove tract points outside the muscle volume
    if DTItracts.fibindex(fibnr,1) < DTItracts.fibindex(fibnr,2)
        sgn = 1;
    else
        sgn = -1;
    end
    p = DTItracts.fibindex(fibnr,1) : sgn : DTItracts.fibindex(fibnr,2);
    in = INSIDE(p);
    
    dsig = diff([1 ~in' 1]);
    startIndex = find(dsig < 0);
    endIndex   = find(dsig > 0)-1;
    len = endIndex - startIndex+1;
    [nSteps,maxIndex] = max(len);
    
    if isempty(len) || nSteps < 5
        % No or too few points are inside the surface
        continue
    else
        % Select the longest continuous section inside the muscle surface
        fibindex_trunc(fibnr,1) = p(startIndex(maxIndex));
        fibindex_trunc(fibnr,2) = p(endIndex(maxIndex));
        length_trunc(fibnr)     = (nSteps-1) * DTItracts.stepsize;
    end
    % Make sure the startpoint (column 1 in fibindex_trunc) has a lower z-value than
    % the endpoint (column2 in fibindex_trunc).
    if DTItracts.tracts_xyz(3,fibindex_trunc(fibnr,1)) > DTItracts.tracts_xyz(3,fibindex_trunc(fibnr,2))
        fibindex_trunc(fibnr,:) = fliplr(fibindex_trunc(fibnr,:));
        DTItracts.fibindex(fibnr,:) = fliplr(DTItracts.fibindex(fibnr,:));
        
    end
end

% Return as separate output arguments if 2 outputs are requested or add to
% structure 'DTItracts' when one output is requested.
if nargout == 1
    DTItracts.fibindex_trunc = fibindex_trunc;
    DTItracts.length_trunc   = length_trunc;
    varargout{1} = DTItracts;
elseif nargout == 2
    varargout{1} = fibindex_trunc;
    varargout{2} = length_trunc;
end
t_elapsed = toc;
fprintf('Time used for truncating fibres: %.2f\n',t_elapsed)
% fprintf('completed\n')
close(hwait)
end % of function

