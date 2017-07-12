function IN = inside_surface(FV,pts)
%%INSIDE_SURFACE calculates whether the points given in pts (n x 3 matrix
%%of x y z coordinates) are inside (1) or outside (0) the triangulated
%%surface given by FV (containing fields faces and vertices)
%
% Bart Bolsterlee, Neuroscience Research Australia (NeuRA)
% February 2017
%
% This algorithm projects a 'ray' from the query point in a given direction
% and counts how many triangles of the surface model are intersected by the
% ray. An odd number of intersections means the point is inside the surface
% and an even number means the point is outside the surface.
%
% The algorithm to calculate the intersection of a ray and a triangle was
% inspired by some linear algebra wizardry which I stole from the c-code by
% Dan Sunday available on:
% http://geomalgorithms.com/a06-_intersect-2.html#intersect3D_RayTriangle()

VERT = FV.vertices;
FAC  = FV.faces;

% Get all edge vectors
V0 = VERT(FAC(:,1),:);
V1 = VERT(FAC(:,2),:);
V2 = VERT(FAC(:,3),:);

IN = false(size(pts,1),1);
for i = 1 : size(pts,1)
    point = pts(i,:);
    
    % The ray is projected along the x-axis from point [px py pz].
    % Triangles that are completely below or above the plane y = py and z =
    % pz can never intersect with this ray, so these faces can be excluded
    % now. This eliminates most faces from the surface and will make the
    % computation much faster.
    
    d  = [1 0 0]; % ray direction vector in positive x-direction
    excl = (V0(:,2) < point(2) & V1(:,2) < point(2) & V2(:,2) < point(2)) | ... % all vertices below y = yp (i.e. triangle is completely below plane y = yp)
        (V0(:,2) > point(2) & V1(:,2) > point(2) & V2(:,2) > point(2)) | ... % all vertices above y = yp
        (V0(:,3) < point(3) & V1(:,3) < point(3) & V2(:,3) < point(3)) | ... % all vertices below z = zp
        (V0(:,3) > point(3) & V1(:,3) > point(3) & V2(:,3) > point(3));      % all vertices above z = zp
    incl = ~excl;
    %
    %     figure;hold on axis equal patch('vertices',VERT,'faces',FAC,...
    %         'FaceAlpha',0.2,'FaceColor','g')
    %     patch('vertices',VERT,'faces',FAC(incl,:),...
    %         'FaceAlpha',1,'FaceColor','r')
    %     %     plot3(point(1),point(2),point(3),'ro',... %
    %     'MarkerFaceColor','r','MarkerEdgeColor','none',... %
    %     'MarkerSize',8) quiver3(point(1),point(2),point(3),...
    %         d(1),d(2),d(3),10,'-','filled','LineWidth',3,... 'Color','r')
    %
    %     view(-34,-9)
    
    % For the leftover faces, calculate whether the ray intersects the
    % face. The following code is the vector implementation of the c-code
    % given at.:
    % http://geomalgorithms.com/a06-_intersect-2.html#intersect3D_RayTriangle()
    nF = sum(incl);
    
    U = V1(incl,:) - V0(incl,:); % edge vectors
    V = V2(incl,:) - V0(incl,:);
    
    UU = sum(U.*U,2);
    UV = sum(U.*V,2);
    VV = sum(V.*V,2);
    DD = UV.*UV - UU.*VV;
    
    % Calculate normal of faces
    N = cross(U,V);
    
    W0 = ones(nF,1) * point - V0(incl,:);
    A = -sum(N .* W0,2);
    B =  sum(N .* (ones(nF,1) * d),2);
    
    % if R < 0 the ray points away from the triangle and no intersection
    % will be found.
    R = A ./ B;
    
    % calculate intersection of ray and plane
    II = ones(nF,1) * point + (R*ones(1,3)) .* (ones(nF,1) * d);
    
    W = II - V0(incl,:);
    WV = sum(W.*V,2);
    WU = sum(W.*U,2);
    S = (UV .* WV - VV .* WU) ./ DD;
    T = (UV .* WU - UU .* WV) ./ DD;
    
    % Count how many times the surface was intersected
    nint = sum((S >= 0 & S <= 1 & T >= 0 & (S+T) <= 1) & R >= 0);
    % return false if even, and true if odd
    IN(i) = logical(mod(nint,2));
    
end

end % of function



