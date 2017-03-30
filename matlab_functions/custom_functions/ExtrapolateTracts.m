function varargout = ExtrapolateTracts( DTItracts,surf_model,varargin )
%EXTRAPOLATETRACTS This function fits a polynomial on the truncated fibre
% tracts and then extrapolates this polynomial linearly at its endpoints
% until the surface is intersected.
%
% Bart Bolsterlee, Neuroscience Research Australia (NeuRA)
% February 2017
%
% ----------------- USAGE ----------------- 
% DTItracts = ExtrapolateTracts(DTItracts,surf_model,order)
% or 
% [PolyCoeff,fibrelength,ext,pct_ext,endpoints,residual] = ExtrapolateTracts(DTItracts,surf_model,order)
%
% ----------------- INPUT -----------------
% Required:
% - DTItracts : a structure array containing at least the fields tracts_xyz 
%               and fibindex_trunc (created with TruncateTracts.m).
% - surf_model: a structure array containing the fields
%               'vertices','faces' with the vertices and faces of the 
%               surface used for extrapolating the tracts.
%
% Optional: 
% - order     : order of the polynomial fit. Default = 3
%
% ----------------- OUTPUT -----------------
% If one output argument is provided DTItracts is returned with the 
% following fields appended (n = number of fibres)
%
% - PolyCoeff   : n x 1 structure array with fields x, y and z containing the
%                 polynomial coefficients of the polynomial curve for the x,
%                 y and z coordinates, respectively. Also has fields t0
%                 and t1 with the first and last t-value of the curve.
% - fibrelength : n x 1 array with fibre length (curved length + extrapolated parts)
% - ext         : n x 2 array with extension (in mm) at each of the endpoints
% - pct_ext     : n x 1 array with total extension as percentage of the total fibre length
% - endpoints   : n x 2 x 3 array with xyz-location of extrapolated endpoints on the surface.
% - residual    : residual of polynomial fit (absolute mean distance of tracts
%                points to nearest point on the polynomial curve).
%
% If multiple output arguments are provided, the outputs are (in this
% order):
%
% [PolyCoeff,fibrelength,ext,pct_ext,endpoints,residual]
%
% Uses the function inside_surface.m to determine whether a point is inside
% or outside the surface.

%% Check inputs
if nargin == 1 || nargin > 3 
    error('Wrong number of input arguments. Usage is ExtrapolateTracts(DTItracts,surf_model) or ExtrapolateTracts(DTItracts,surf_model,order)')
end

tic
p = inputParser;
addRequired(p,'DTItracts',@isstruct)
addRequired(p,'surf_model',@isstruct)
addOptional(p,'order',3,@(x)(isnumeric(x) && mod(x,1)==0 && x>0))
parse(p,DTItracts,surf_model,varargin{:})
order = p.Results.order;

% Check if all required fields are available in DTItracts and surface
if ~isfield(DTItracts,'tracts_xyz')
    error('Required field tracts_xyz not found in DTItracts.')
end
if ~isfield(DTItracts,'fibindex_trunc')
    error('Required field fibindex_trunc not found in DTItracts.')
end

if ~isfield(surf_model,'vertices')
    error('Required field vertices not found in surf_model.')
end

if ~isfield(surf_model,'faces')
    error('Required field faces not found in surf_model.')
end

%%
% Number of fibres
nFib = size(DTItracts.fibindex_trunc,1);

% Use the truncated tracts for fitting a polynomial and extrapolation.
if ~isfield(DTItracts,'fibindex_trunc')
    error('''fibindex_trunc'' is a required field in DTItracts')
end

endpoints   = NaN(nFib,2,3);
ext         = NaN(nFib,2);
fibrelength = NaN(nFib,1);
residual    = NaN(nFib,1);
hwait = waitbar(0,'','Name','Progress bar ExtrapolateTracts');
for fibnr = 1:1:nFib
    waitbar(fibnr/nFib,hwait,sprintf('Extrapolating fibre %d of %d',fibnr,nFib))
    
    first = DTItracts.fibindex_trunc(fibnr,1);
    last  = DTItracts.fibindex_trunc(fibnr,2);
    if any(isnan([first,last]))
        continue
    end
    if first < last
        sgn = 1;
    else
        sgn = -1;
    end
    tractpoints = [DTItracts.tracts_xyz(1,first:sgn:last)' ...
                   DTItracts.tracts_xyz(2,first:sgn:last)' ...
                   DTItracts.tracts_xyz(3,first:sgn:last)'];

    nPoints = size(tractpoints,1);
    T = zeros(nPoints,order+1);
    for o = order:-1:0
        T(:,order-o+1) = (0:nPoints-1).^o;
    end

    % Fit the polynomial
    coeff.x = (T \ tractpoints(:,1))';
    coeff.y = (T \ tractpoints(:,2))';
    coeff.z = (T \ tractpoints(:,3))';
    coeff.t0 = 0;
    coeff.t1 = nPoints-1;   
    PolyCoeff(fibnr,1) = coeff;            
    
%     if nargout > 4
    % Calculate residual distance of data points to polynomial fit
    [~,dist] = FindNearestT(coeff,tractpoints);
    residual(fibnr) = abs(mean(dist));

%     end
   
% Find projection of linearly extrapolated line from the endpoints to the
% surface model
        % ----- Endpoint 1 -----
        p1 =  [polyval(coeff.x,0),...
               polyval(coeff.y,0),...
               polyval(coeff.z,0)];
        d1 = [polyval(polyder(coeff.x),coeff.t0),...
              polyval(polyder(coeff.y),coeff.t0),...
              polyval(polyder(coeff.z),coeff.t0)];
%         d1 = d1 / norm(d1); 

        if inside_surface(surf_model,p1) == 0
            disp('point outside surface')
        end
        
        tract_dir1 =  [polyval(coeff.x,1),polyval(coeff.y,1),polyval(coeff.z,1)] - p1;
        if sign(dot(d1,tract_dir1)) == 1
            % Change sign of direction
            d1 = -d1;
        end

        % ----- Endpoint 2 -----
        p2 =  [polyval(coeff.x,nPoints-1),...
               polyval(coeff.y,nPoints-1),...
               polyval(coeff.z,nPoints-1)];
        d2 = [polyval(polyder(coeff.x),coeff.t1),...
              polyval(polyder(coeff.y),coeff.t1),...
              polyval(polyder(coeff.z),coeff.t1)];
%         d1 = d1 / norm(d1); 
        
        tract_dir2 =  [polyval(coeff.x,nPoints-2),polyval(coeff.y,nPoints-2),polyval(coeff.z,nPoints-2)] - p2;
        if sign(dot(d2,tract_dir2)) == 1
            % Change sign of direction
            d2 = -d2;
        end
        
        % Calculate intersection of projection of endpoint1 with surface model
        endpoints(fibnr,:,:) = intersectRaySurface(surf_model,[p1;p2],[d1;d2]);

% %%       Some code that could be uncommented for diagnostic purposes        
%         figure;
%         hold on
%         axis equal
%         patch('vertices',surf_model.vertices,'faces',surf_model.faces,...
%             'FaceAlpha',0.2,'FaceColor','g')
%         t_plot = linspace(coeff.t0,coeff.t1,100);
% 
%         plot3(polyval(coeff.x,t_plot),...
%             polyval(coeff.y,t_plot),...
%             polyval(coeff.z,t_plot),...
%             'LineWidth',3,'Color','b')
%         view(-34,-9)
% 
%         plot3([p1(1) endpoints(fibnr,1,1)],[p1(2) endpoints(fibnr,1,2)],[p1(3) endpoints(fibnr,1,3)],...
%             'LineWidth',2,'Color','r','MarkerSize',8,'Marker','o')
%         
%         plot3([p2(1) endpoints(fibnr,2,1)],[p2(2) endpoints(fibnr,2,2)],[p2(3) endpoints(fibnr,2,3)],...
%             'LineWidth',2,'Color','r','MarkerSize',8,'Marker','o')
% 
%%    
    % Calculate the total length of the fibre, which is the sum of the
    % extrapolated parts and the length along the polynomial-fitted curve.
    % 1) Length of polynomial fitted section (this should be very close to the
    % length of the tract on which the polynomial was fitted).
    
    [ CurveLength] = FuncCurveLength( coeff,0,nPoints-1 );
    ext(fibnr,:) = sqrt(sum(([p1;p2] - squeeze(endpoints(fibnr,:,:))).^2,2))';
    fibrelength(fibnr) = CurveLength + sum(ext(fibnr,:));

end

% Calculate the extension as percentage of the total fibre length.
pct_ext = sum(ext,2) ./ fibrelength * 100;
%% Create output arguments
% If only one output argument is provided, add outputs to the structure array
% DTItracts as fields.
if nargout == 1
    DTItracts.PolyCoeff   = PolyCoeff;
    DTItracts.fibrelength = fibrelength;
    DTItracts.ext         = ext;
    DTItracts.pct_ext     = pct_ext;
    DTItracts.endpoints   = endpoints;
    DTItracts.residual    = residual ;
    varargout{1}          = DTItracts;
end
if nargout > 1
    varargout{1} = PolyCoeff;
    varargout{2} = fibrelength;
    if nargout > 2
        varargout{3} = ext;
        if nargout > 3
            varargout{4} = pct_ext;
            if nargout > 4
                varargout{5} = endpoints;
                if nargout > 5
                    varargout{6} = residual; 
                end
            end
        end
    end
end
    

t_elapsed = toc;
fprintf('Time used for polynomial fitting and extrapolating fibres: %.2f\n',t_elapsed)
close(hwait)
end % of function
