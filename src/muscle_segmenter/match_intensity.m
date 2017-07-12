function [V2_opt, dist_along_n,best_fit_metric] = ...
    match_intensity(V1,V2,img1,img2, r1,r2,ds,opt_metric )
%MATCH_INTENSITY

% Make vector with locations along the normal to sample from
rng1 = -r1 : ds : r1;
rng2 = -r2 : ds : r2;

% Calculate normal vectors

[~,~,~,~,N1] = fit_closed_curve( V1);
[~,~,~,~,N2] = fit_closed_curve( V2);

% Calculate intensity profile in reference image
I_prof1 = calc_int( V1,N1,img1,rng1 );

% Calculate intensity profile in new image
I_prof2 = calc_int( V2,N2,img2,rng2 );

% Find best match
n1 = length(rng1);
n2 = length(rng2);

metric = zeros(size(I_prof1,1),n2-n1+1);
for ii = 1 : 1 : (n2-n1+1)
    %     A = I_prof2(:,(0:n1-1)+ii);
    switch ndims(img1)
        case 2
            A = I_prof2(:,(0:n1-1)+ii);
        case 3
            A = I_prof2(:,(0:n1-1)+ii,:);
    end
    switch opt_metric
        
        case 'corr'
            % Calculate the Pearson's correlation coefficient.
            % dev_meanA = bsxfun(@minus, A, nanmean(A,2)); % zero-mean all profiles
            err  = nansum(A .* I_prof1,2) ./ ...
                (sqrt(nansum(A.^2,2)) .* sqrt(nansum(I_prof1.^2,2)));
            if ndims(err) == 3
                err = sum(err,3);
            end
            metric(:,ii) = err;
            
        case 'lsq'
            % mean square difference as error metric.
            err = nanmean((A - I_prof1).^2,2);
            if ndims(err) == 3
                err = mean(err,3);
            end
            metric(:,ii) = err;
    end
end
% Calculate where the sum squared difference is minimal.
switch opt_metric
    case 'corr'
        % Calculate where the correlation coefficient is maximal
        [best_fit_metric, idx] = max(metric,[],2);
    case 'lsq'
        % Calculate where the sum squared difference is minimal.
        [best_fit_metric, idx] = min(metric,[],2);
end
% [best_fit_retric, idx] = min(metric,[],2);
dist_along_n = (rng2(idx) - rng1(1))';

% calculate new locations
V2_opt = V2 + N2 .* (dist_along_n*[1 1]);



end

