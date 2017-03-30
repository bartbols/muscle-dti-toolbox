function penangle = CalcPenAngle( DTItracts,surf_model,varargin )
%CALCPENANGLE3D calculates the pennation angle of fibres described by the
%polynomial coefficient in PolyCoeff and endpoints, which are both fields in
%DTItracts, with the surface model surf_model, described by its fields
%'vertices' and 'faces'
%
% Bart Bolsterlee, Neuroscience Research Australia (NeuRA)
% February 2017
%
% ----------------- USAGE -----------------
% DTItracts = TruncateTracts( DTItracts,surf_model)
% or 
% DTItracts = TruncateTracts( DTItracts,surf_model,radius)
%
% ----------------- INPUT -----------------
% Required
% - DTItracts  : a structure array containing at least the fields PolyCoeff
%               and endpoints, as generated with ExtrapolateTracts
% - surf_model : a structure array 'surface' containing the fields
%                'vertices' and 'faces' and optionally 'normals'
%
% Optional
% - radius    : the radius around the fibre endpoint in which triangles are
%               included for pennation angle calculation. Default = 2.5
%
% ----------------- OUTPUT -----------------
% penangle    : n x 2 array of pennation angle of endpoint 1 and 2. n is
%               the number of fibres

%% Check inputs
p = inputParser;
addRequired(p,'DTItracts',@isstruct)
addRequired(p,'surf_model',@isstruct)
addOptional(p,'radius',2.5,@(x) validateattributes(x,{'numeric'},{'scalar'}) )
parse(p,DTItracts,surf_model,varargin{:})
radius = p.Results.radius;

%% Calculate the centre of all triangles
tic
centre = squeeze(mean(reshape(surf_model.vertices(surf_model.faces,:)',3,size(surf_model.faces,1),3),3))';

% Calculate face normals if not provided
if ~isfield(surf_model,'normals')
    surf_model.normals = facenormals(surf_model);
end

% Select triangles within 2.5 mm of the current point
radius = 2.5;
nFibres = size(DTItracts.endpoints,1);
penangle = NaN(nFibres,2);
for fibnr = 1:nFibres
    for ep = 1:2
        tmp = centre - repmat(squeeze(DTItracts.endpoints(fibnr,ep,:))',size(centre,1),1);
        dist2centre = sqrt(sum(tmp.^2,2));
        idx = dist2centre < radius;
        
        if sum(idx) == 0
            continue
        end
        
        if ep == 1
            t = DTItracts.PolyCoeff(fibnr).t0;
        else
            t = DTItracts.PolyCoeff(fibnr).t1;
        end
        slope = [polyval(polyder(DTItracts.PolyCoeff(fibnr).x),t) ,...
            polyval(polyder(DTItracts.PolyCoeff(fibnr).y),t) ,...
            polyval(polyder(DTItracts.PolyCoeff(fibnr).z),t)];
        slope = slope/norm(slope); % make unit vector
        
        % calculate angle with normal vectors of selected triangles
        penangle(fibnr,ep) = abs(mean(asind(sum(surf_model.normals(idx,:) .* repmat(slope,sum(idx),1),2))));
%         penangle(fibnr,ep) = abs(median(asind(sum(surface.normals(idx,:) .* repmat(slope,sum(idx),1),2))));
        
        %                 patch('vertices',surface.vertices,'faces',surface.faces(idx,:),'facecolor','g')
        %                 quiver3(surface.centre(idx,1),surface.centre(idx,2),surface.centre(idx,3),...
        %                         surface.normals(idx,1),surface.normals(idx,2),surface.normals(idx,3))
    end
    
end

t_elapsed = toc;
fprintf('Time used for calculating pennation angle: %.2f\n',t_elapsed)
end % of function
