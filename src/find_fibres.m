function [ idx ] = find_fibres( tracts_xyz,fibindex,p,r )
%FIND_FIBRES finds the indices of all fibres provided by tracts_xyz and
%fibindex with at least one tract point within the radius 'r' from point 'p' 
%(3x1 vector).
%
% Bart Bolsterlee, Neuroscience Research Australia
% November 2017
%
% USAGE:
% idx = find_tracts(tracts_xyz,fibindex,p,radius)

% Find indices of tractpoints within the defined distance from p.
A = find(sqrt(sum((tracts_xyz - p * ones(1,size(tracts_xyz,2))).^2,1)) < r);

% Now loop through all the selected tract points and find the fibre they 
% belong to.
B   = sort(fibindex,2);
idx = zeros(size(A));

%%
% % This vectorized version gives the same answer but is slower than the loop
% BB1 = B(:,1) * ones(1,length(A));
% BB2 = B(:,2) * ones(1,length(A));
% AA = ones(size(B,1),1) * A;
% [idx,~] = find( sign(BB1 - AA) ~= ...
%                 sign(BB2 - AA));
           
for ii = 1 : length(A)
    tmp = find(A(ii) >= B(:,1)  & A(ii) <= B(:,2) );
    if ~isempty(tmp)
        idx(ii) = tmp;
    end
end
% % Remove zeros
idx(idx==0) = [];
% % Remove duplicates fibres (because they have multiple tract points within
% % radius r around P and were therefore selected multiple times.)
idx = unique(idx);
end % of function

