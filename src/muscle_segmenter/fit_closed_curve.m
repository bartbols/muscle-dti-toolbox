function [Vs,Ns,ts,t,N] = fit_closed_curve( V,n )
%FIT_CLOSED_CURVE uses csape to fit a closed curve around the points given
%by xi and yi. The curve is then interpolated at n points. At these
%locations, the normal vectors are calculated as well.
% Vs : location of points sampled along the fitted curve
% Ns : normal vectors of points sampled along the fitted curve
% t  : t-values of control points
% ts : t-values of points sampled along the curve
% N  : normal vectors at control point location


% t-values of control points
t = cumsum(sqrt([0,diff(V(:,1)')].^2 + [0,diff(V(:,2)')].^2));

if nargin == 1
    % define number of samples based on the length of the curve so that
    % each segment is approximately 0.5 voxel long.
    n = 2* ceil(max(t));
end

px1 = csape(t,V(:,1),'periodic');
py1 = csape(t,V(:,2),'periodic');

% t-values of sampled points along the curve
ts = linspace(0,max(t),n);
Vs(:,1) = fnval(px1,ts);
Vs(:,2) = fnval(py1,ts);

if nargout > 1
    Ns(:,1) = -fnval(fnder(py1),linspace(0,max(t),n));
    Ns(:,2) =  fnval(fnder(px1),linspace(0,max(t),n));
    
        % normalize to unit vector
        Ns = Ns ./ (sqrt(sum(Ns.^2,2))*[1 1]);
    if nargout > 4
        N(:,1) = -fnval(fnder(py1),t);
        N(:,2) =  fnval(fnder(px1),t);
        
        % normalize to unit vector
        N = N ./ (sqrt(sum(N.^2,2))*[1 1]);
    end
end

end

