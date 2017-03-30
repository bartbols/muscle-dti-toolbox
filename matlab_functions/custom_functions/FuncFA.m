function FA = FuncFA( L1,L2,L3 )
%FUNCFA Calculates the fractional anisotropy FA given the eigenvalues L1, L2
%and L3.
%
% Bart Bolsterlee, Neuroscience Research Australia (NeuRA)
% February 2017
% ----------------- USAGE ----------------- 
% FA = FuncFA( L1,L2,L3 )
%
% ----------------- INPUT ----------------- 
% L1: array (1-D or multi-dimensional) with first eigenvalue
% L2: array (1-D or multi-dimensional) with second eigenvalue
% L3: array (1-D or multi-dimensional) with third eigenvalue
% 
% ----------------- OUTPUT ----------------- 
% FA: array (same dimensions as input) with fractional anisotropy

p = inputParser;
addRequired(p,'L1',@(x) validateattributes(x,{'numeric'},{'nonempty'}))
addRequired(p,'L2',@(x) validateattributes(x,{'numeric'},{'nonempty'}))
addRequired(p,'L3',@(x) validateattributes(x,{'numeric'},{'nonempty'}))
parse(p, L1,L2,L3);


FA = sqrt(1/2) * sqrt((L1 - L2).^2 + (L2 - L3).^2 + (L3 - L1).^2 ) ./ ...
                         sqrt(L1.^2 + L2.^2 + L3.^2);



end

