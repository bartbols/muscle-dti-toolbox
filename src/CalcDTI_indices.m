function varargout = CalcDTI_indices( DTItracts, fib_fname )
%CALCDTI_INDICES Loads the eigenvalue data from the the .fib.gz-file and
% interpolates the tracts to obtain FA, MD, lambda1, lambda2 and lambda3
% values along each of the tracts in 'DTItracts'. The mean value for a
% tracts is returned.
%
% Bart Bolsterlee, Neuroscience Research Australia (NeuRA)
% February 2017
%
% ----------------- USAGE -----------------
% DTItracts = CalcDTI_indices( DTItracts,FIB )
% or
% [fa,md,l1,l2,l3] = CalcDTI_indices( DTItracts,FIB )
%
% ----------------- INPUT -----------------
% - DTItracts : a structure array containing at least the fields tracts
%               and fibindex_trunc.
% - fib_fname : the filename of the fibre file created with DSI studio.
%               Must have extension .fib.gz
%
% ----------------- OUTPUT -----------------
% If one output argument is provided DTItracts is returned with the 
% following fields appended (n = number of fibres)
% fa          : n x 1 array of mean FA values for each fibre 
% md         : n x 1 array of mean Mean Diffusivity values for each fibre
% l1          : n x 1 array of mean primary eigenvalue for each fibre
% l2          : n x 1 array of mean secondary eigenvalue for each fibre
% l3          : n x 1 array of mean tertiary eigenvalue for each fibre
%
% If multiple output arguments are provided, the outputs are (in this
% order):
%
% [ fa,md,lambda1,lambda2,lambda3 ]
%

%% Check inputs
if nargin ~= 2
    error('Wrong number of input arguments. Usage is CalcTractIndices( DTItracts, fib_fname )')
end

p = inputParser;
addRequired(p,'DTItracts',@isstruct)
addRequired(p,'fib_fname',@(x) ~isempty(strfind(x,'.fib.gz')))
parse(p,DTItracts,fib_fname)

% Check if all required fields are available in DTItracts and surface
if ~isfield(DTItracts,'tracts')
    error('Required field tracts not found in DTItracts.')
end
if ~isfield(DTItracts,'fibindex_trunc')
    error('Required field fibindex_trunc not found in DTItracts.')
end

if exist(fib_fname,'file') ~= 2
    error('Fibre-file %s not found.',fib_fname)
end
%%
% Load the fibre file by first unzipping, then loading as a MAT-file.
% Delete the unzipped fibre-file again.
tic
hwait = waitbar(0,'Loading fibre file...','Name','Progress bar CalcDTI_indices');
gunzip(fib_fname);
data = load(fib_fname(1:end-3),'-mat');
delete(fib_fname(1:end-3))

% The FA, MD and eigenvalue maps are stored as a vector in the fibre-file 
% and need to be reshaped to form the 3D map, which can then be interpolated at
% the points along the fibre tracts.
%
% Note: DSI studio calls the primary eigenvalue 'axial_dif', the secondary
% eigenvalue 'radial_dif1' and the tertiary eigenvalue 'radial_dif2'.

FA_map  = reshape(data.fa0,data.dimension);
MD_map  = reshape(data.md,data.dimension);
L1_map  = reshape(data.ad,data.dimension);
L2_map  = reshape(data.rd1,data.dimension);
L3_map  = reshape(data.rd2,data.dimension);

nFib    = size(DTItracts.fibindex_trunc,1);
fa      = NaN(nFib,1);
md      = NaN(nFib,1);
lambda1 = NaN(nFib,1);
lambda2 = NaN(nFib,1);
lambda3 = NaN(nFib,1);

clear data
for fibnr = 1 : 1 : nFib
    waitbar(fibnr/nFib,hwait,sprintf('Calculating DTI indices for fibre %d of %d',fibnr,nFib))
    first = DTItracts.fibindex_trunc(fibnr,1);
    last  = DTItracts.fibindex_trunc(fibnr,2);
    if any(isnan([first,last]))
        continue
    end
    if first < last
        d = 1;
    else
        d = -1;
    end
    
    % Interpolate maps of the DTI index at the tract point. Add 1 because
    % the tract points are indexed from 0 (voxel 1  = index 0), while the
    % maps used for interpolation are indexed from 1.
    x = DTItracts.tracts(1,first:d:last)+1;
    y = DTItracts.tracts(2,first:d:last)+1;
    z = DTItracts.tracts(3,first:d:last)+1;
    fa(fibnr)      = nanmean(interp3(FA_map,y,x,z));
    md(fibnr)     = nanmean(interp3(MD_map,y,x,z));
    lambda1(fibnr) = nanmean(interp3(L1_map,y,x,z));
    lambda2(fibnr) = nanmean(interp3(L2_map,y,x,z));
    lambda3(fibnr) = nanmean(interp3(L3_map,y,x,z));
    
end
t_elapsed = toc;
fprintf('Time used for calculating DTI indices: %.2f\n',t_elapsed)

close(hwait)

%% Create output arguments
if nargout == 1
    % Add outputs to DTI tracts
    DTItracts.fa      = fa;
    DTItracts.md      = md;
    DTItracts.lambda1 = lambda1;
    DTItracts.lambda2 = lambda2;
    DTItracts.lambda3 = lambda3;
    varargout{1} = DTItracts;
end
if nargout > 1
    % Provide DTI indices as separate outputs.
    varargout{1} = fa;
    varargout{2} = md;
    if nargout > 2
        varargout{3} = lambda1;
        if nargout > 3
            varargout{4} = lambda2;
            if nargout > 4
                varargout{5} = lambda3;
                if nargout > 5
                    error('Too many output arguments. Usage is: [fa,md,l1,l2,l3] = CalcDTI_indices( DTItracts,FIB )')
                end
            end
        end
    end
end
end % of function



