function varargout = plotTensor(tensor,varargin)
%PLOTTENSOR plots the tensor field 'tensor' as ellipsoids. Support 1D, 2D
%and 3D vector fields.
%
%-USE--------------------------------------------------------------------------------
% example 1: plotTensor(tensor)
% where tensor is of size 3x3 or ix3x3xN, ixjx3x3 or ixjxk3x3 tensor field
%
% example 2: plotTensor(tensor,'spacing',<spacing>)
% where spacing is a scalar (for isotropic voxels) or nx1 vector (for
% anisotropic voxel sizes) that controls the size of a voxel in the field.
% Default: spacing=1
%
%
% Bart Bolsterlee
% Neuroscience Research Australia
% May 2021
% (with inspiration from PlotDTI in the fanDTasia Toolbox)

% Read input arguments
p = inputParser;
addRequired(p,'tensor')
addParameter(p,'spacing',[])
addParameter(p,'scale',[])
addParameter(p,'X',[])
addParameter(p,'Y',[])
addParameter(p,'Z',[])
addParameter(p,'np',10) % detail of ellipsoid patch
addParameter(p,'mask',[])
addParameter(p,'type','ev1') % plot eigenvector / diffusion ellipsoid
addParameter(p,'linewidth',2,@isscalar)
addParameter(p,'color','ev1')
addParameter(p,'cdata',[])
parse(p,tensor,varargin{:});

spacing  = p.Results.spacing;
scale    = p.Results.scale;
X        = p.Results.X;
Y        = p.Results.Y;
Z        = p.Results.Z;
np       = p.Results.np;
mask     = p.Results.mask;
plot_type = p.Results.type;
color    = p.Results.color;
custom_cdata = p.Results.cdata;

if isempty(spacing);spacing = 1;end
if isempty(scale);scale=1;end

% Reshape tensor field into one long array of tensors for easy processing.
sz = size(tensor);
nd = ndims(tensor);

field_dim = sz(1:nd-2);
n      = prod(field_dim);
tensor = reshape(tensor,[n,3,3]);

% Build up grid of center locations of tensor ellipsoids
if isempty(X) && isempty(Y) && isempty(Z) % full grid not provided
    if isscalar(spacing)
        spacing = ones(1,3)*spacing;
    end
    if length(field_dim) == 1
        nx=field_dim(1);
        ny=1;
        nz=1;
    elseif length(field_dim) == 2
        nx=field_dim(1);
        ny=field_dim(2);
        nz=1;
    elseif length(field_dim) == 3
        nx=field_dim(1);
        ny=field_dim(2);
        nz=field_dim(3);
    end
    [X,Y,Z] = meshgrid(...
        (0:nx-1)*spacing(1),...
        (0:ny-1)*spacing(2),...
        (0:nz-1)*spacing(3));
end

% If any of the tensor components is a NaN, remove the tensor from the list
% to be plotted.
notNaN_mask = ~any(any(isnan(tensor),3),2);
if ~isempty(mask)
    mask = mask(:) & notNaN_mask;
else
    mask = notNaN_mask;
end

% Keep only values inside the mask.
Xv = X(mask==1);
Yv = Y(mask==1);
Zv = Z(mask==1);
tensor = tensor(mask(:)==1,:,:);
n = size(tensor,1);

if isnumeric(color)
    color2 = color;
    color = 'uniform';
end

if strcmp(color,'fa')
    fa = tensor2fa(tensor);
    % elseif strcmp(color,'l1')
    %     l1 = tensor2lambda(tensor);
    % elseif strcmp(color,'l2')
    %     {~,l2] = tensor2lambda(tensor);
    % elseif strcmp(color,'l3')
    %     [~,~,l3] = tensor2lambda(tensor);
    % elseif strcmp(color,'md')
    %     [~,~,~,md] = tensor2lambda(tensor);
end

first = 0;
cdata=[];
for i=1:n
    
    % Get current tensor
    S = squeeze(tensor(i,:,:));
    if any(any(isnan(S)));continue;end
    
    % extract eigenvalues and -vectors of the current tensor.
    [v,l]=eig(S);
    
    % Sort in descending order.
    [l,idx] = sort(diag(l),'descend');
    v = v(:,idx);
    
    switch plot_type
        case {'ev1','ev2','ev3'} % plot only primary eigenvector
            if first == 0
                P = NaN(3*n,3);
                switch color
                    case  {'ev1','ev2','ev3'}
                        cdata       = NaN(3*n,3);
                    case {'fa','az','l1','l2','l3','md','custom'}
                        cdata       = NaN(3*n,1);
                end
                
                first = 1;
            end
            nr = str2double(plot_type(3));
            P((3*i)-2,:)           = [Xv(i) Yv(i) Zv(i)]-0.5*v(:,nr)'*scale;
            P((3*i)-1,:)           = [Xv(i) Yv(i) Zv(i)]+0.5*v(:,nr)'*scale;
            switch color
                case 'ev1'
                    cdata(3*(i-1)+(1:3),1:3) = repmat(abs(v(:,1)'),3,1);
                case 'ev2'
                    cdata(3*(i-1)+(1:3),1:3) = repmat(abs(v(:,2)'),3,1);
                case 'ev3'
                    cdata(3*(i-1)+(1:3),1:3) = repmat(abs(v(:,3)'),3,1);
                case 'fa'
                    cdata(3*(i-1)+(1:3),1) = fa(i);
                case 'az'
                    sgn = sign(v(3,1));
                    cdata(3*(i-1)+(1:3),1) = atan2(v(2,1).*sgn,v(1,1).*sgn);
                case 'custom'
                    cdata(3*(i-1)+(1:3),1) = custom_cdata(i);
            end
            
        case 'ell' % plot ellipsoid for current tensor
            if first == 0
                % Create unit ellipsoid and empty matrices.
                [Xe,Ye,Ze]  = ellipsoid(0,0,0,1/2,1/2,1/2,np);
                face_offset = 0;
                nv          = numel(Xe); % number of vertices per ellipsoid
                nf          = np^2;      % number of faces per ellipsoid
                vertices    = NaN(nv*n,3);
                faces       = NaN(nf*n,4);
                switch color
                    case {'ev1','ev2','ev3'}
                        cdata       = NaN(nv*n,3);
                    case {'fa','az','l1','l2','l3','md','custom'}
                        cdata       = NaN(nv*n,1);
                end
                first       = 1;
            end
            % Scale unit ellipsoid by the eigenvalues.
            Xc = Xe * l(1) * scale;
            Yc = Ye * l(2) * scale;
            Zc = Ze * l(3) * scale;
            
            % Rotate by rotation matrix formed by eigenvectors.
            Xr = v(1,1)*Xc + v(1,2)*Yc + v(1,3)*Zc;
            Yr = v(2,1)*Xc + v(2,2)*Yc + v(2,3)*Zc;
            Zr = v(3,1)*Xc + v(3,2)*Yc + v(3,3)*Zc;
            
            % Offset by location in grid.
            Xp = Xr + Xv(i);
            Yp = Yr + Yv(i);
            Zp = Zr + Zv(i);
            
            % Create a patch object.
            fvc = surf2patch(Xp,Yp,Zp);
            
            % Matlab's graphics are much faster when all ellipsoids are placed into
            % one patch object rather than one object per ellipsoid, so all
            % ellipsoids are stacked together here.
            
            % Add current ellipsoid to patch structure.
            faces((1:nf) + (i-1)*nf,:)    = fvc.faces+face_offset;
            vertices((1:nv) + (i-1)*nv,:) = fvc.vertices;
            
            switch color
                case 'ev1'
                    %       Color the ellipsoid with the direction of the primary eigenvector.
                    cdata((1:nv) + (i-1)*nv,:)    = repmat(abs(v(:,1)'),nv,1);
                case 'ev2'
                    %       Color the ellipsoid with the direction of the secondary eigenvector.
                    cdata((1:nv) + (i-1)*nv,:)    = repmat(abs(v(:,2)'),nv,1);
                case 'ev3'
                    %       Color the ellipsoid with the direction of the secondary eigenvector.
                    cdata((1:nv) + (i-1)*nv,:)    = repmat(abs(v(:,3)'),nv,1);
                case 'fa'
                    % Color according to fractional anisotropy
                    cdata((1:nv) + (i-1)*nv)      = fa(i);
                case 'az'
                    % Azimuthal angle
                    sgn = sign(v(3,1));
                    cdata((1:nv) + (i-1)*nv)      = atan2(v(2,1).*sgn,v(1,1).*sgn);
                case 'custom'
                    cdata((1:nv) + (i-1)*nv)      = custom_cdata(i);
            end
            face_offset                   = face_offset + nv;
    end
end

switch plot_type
    case {'ev1','ev2','ev3'} % plot only eigenvector lines
        
        if isempty(cdata)
            cdata = repmat(color2,size(P,1),1);
        end
        handle = patch(P(:,1),P(:,2),P(:,3),'k',...
            'EdgeColor','flat',...
            'FaceVertexCData',cdata,...
            'linewidth',p.Results.linewidth);
    case 'ell' % plot ellipsoid
        if isempty(cdata)
            cdata = repmat(color2,nv*n,1);
        end
        handle = patch('Vertices',vertices,...
            'Faces',faces,...
            'FaceColor','flat',...
            'EdgeColor','none',...
            'FaceVertexCData',cdata);
end
if any(strcmp(color,{'fa','az'}))
    colorbar
end
axis equal
view([30 45]);
lighting gouraud

% Create light that moves with camera.
c = camlight('headlight','infinite');
hr = rotate3d;                 % Create rotate3d-handle
hr.ActionPostCallback = @RotationCallback; % assign callback-function
hr.Enable = 'on';              % no need to click the UI-button


if nargout > 0
    varargout{1} = handle;
end


% Sub function for callback
    function RotationCallback(~,~)
        c = camlight(c,'headlight');
    end

end % of function