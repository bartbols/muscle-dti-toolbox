function [ idx ] = find_fibres( DTItracts,p,radius )
%FIND_FIBRES finds the indices of all fibres in 'DTItracts' with at least 
% one tract point with the radius 'r' from point 'p' (3x1 vector). The
% structure DTItracts should at least contain the fields tracts_xyz and
% fibindex_trunc (i.e. the tracts should have been truncated by
% TruncateTracts.m)
%
% Bart Bolsterlee, Neuroscience Research Australia
% November 2017
%
% USAGE:
% idx = find_tracts(DTItracts,p,radius)

% Find indices of tractpoints within the defined distance from p.
tract_idx = find(sqrt(sum((DTItracts.tracts_xyz - p * ones(1,size(DTItracts.tracts_xyz,2))).^2,1)) < radius);

% Now loop through all the selected tract points and find the fibre they 
% belong to.
fibindex_sorted = sort(DTItracts.fibindex_trunc,2);
idx             = zeros(size(tract_idx));

for i = tract_idx
    tmp = find(i >= fibindex_sorted(:,1)  & i <= fibindex_sorted(:,2) );
    if ~isempty(tmp)
        idx(i) = tmp;
    end
end
idx(idx==0) = [];
idx = unique(idx);

end % of function

