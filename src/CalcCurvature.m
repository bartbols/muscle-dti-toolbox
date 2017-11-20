function curvature = CalcCurvature( PolyCoeff)
%CALCCURVATURE calculates the curvature of all fibres in input 
% structure PolyCoeff. The curvature is the mean curvature of 100 points
% sampled along the polynomial curve at equal intervals of t.
%
% Bart Bolsterlee, Neuroscience Research Australia (NeuRA)
% February 2017
% ----------------- USAGE ----------------- 
% kappa = CalcCurvature(PolyCoeff)
%
% ----------------- INPUT -----------------
% - PolyCoeff : n x 1 structure array with fields x,y, z, t0 and t1 (output
%               of ExtrapolateTracts)
%         
% ----------------- OUTPUT -----------------
% - curvature : n x 1 array of mean curvatures per fibre.
%
nFib = length(PolyCoeff);

hwait = waitbar(0,'',...
    'Name','Progress bar CalcCurvature');

tic
curvature = NaN(nFib,1);
for fibnr = 1:nFib
    if isempty(PolyCoeff(fibnr).x)
        continue
    end
    waitbar(fibnr/nFib,hwait,sprintf('Calculating curvature of fibre %d of %d',fibnr,nFib))

    % sample at 100 equal intervals of t
    t = linspace(PolyCoeff(fibnr).t0,...
        PolyCoeff(fibnr).t1,100);
    curv = CalcKappa(PolyCoeff(fibnr),t);
    curvature(fibnr) = mean(curv) * 1000; % convert from 1/mm to 1/m
end
t_elapsed = toc;
fprintf('It took %.2f seconds to calculate curvatures.\n',t_elapsed)

close(hwait)
end

