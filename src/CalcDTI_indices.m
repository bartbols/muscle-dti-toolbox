function varargout = CalcDTI_indices( DTItracts, fib_fname,varargin)
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
% Optional, as 'parameter',<value> pair:
% - flip      : scalar value 1, 2 or 3 indicating the dimension to flip the
%               image over prior to interpolation.
% 
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
p = inputParser;
addRequired(p,'DTItracts',@isstruct)
addRequired(p,'fib_fname',@(x) contains(x,'.fib'))
addParameter(p,'flip',[],@isscalar)
parse(p,DTItracts,fib_fname,varargin{:})

flip2 = p.Results.flip;

if isempty(flip2)
    if isfield(DTItracts.TrackSettings,'algorithm')
        if strcmp(DTItracts.TrackSettings.algorithm,'matlab')
            % When fibre tracking is performed in MATLAB from DSI-studio
            % derived FIB-files, the parameter maps needs to be flipped
            % around the second dimension prior to interpolation.
            flip2 = 2;
        end
    end
end
    

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

if ~isempty(flip2) && flip2 ~= 0
    % Flip the dimensions
    FA_map = flip(FA_map,flip2);
    MD_map = flip(MD_map,flip2);
    L1_map = flip(L1_map,flip2);
    L2_map = flip(L2_map,flip2);
    L3_map = flip(L3_map,flip2);
end


nFib    = size(DTItracts.fibindex_trunc,1);
fa      = NaN(nFib,1);
md      = NaN(nFib,1);
lambda1 = NaN(nFib,1);
lambda2 = NaN(nFib,1);
lambda3 = NaN(nFib,1);

clear data
waitbar(0,hwait,sprintf('Calculating DTI indices for %d fibres',nFib))

% Interpolate maps of the DTI index at the tract points. Add 1 because
% the tract points are indexed from 0 (voxel 1  = index 0), while the
% maps used for interpolation are indexed from 1.
FA = interp3(FA_map,DTItracts.tracts(2,:)+1,DTItracts.tracts(1,:)+1,DTItracts.tracts(3,:)+1);
MD = interp3(MD_map,DTItracts.tracts(2,:)+1,DTItracts.tracts(1,:)+1,DTItracts.tracts(3,:)+1);
L1 = interp3(L1_map,DTItracts.tracts(2,:)+1,DTItracts.tracts(1,:)+1,DTItracts.tracts(3,:)+1);
L2 = interp3(L2_map,DTItracts.tracts(2,:)+1,DTItracts.tracts(1,:)+1,DTItracts.tracts(3,:)+1);
L3 = interp3(L3_map,DTItracts.tracts(2,:)+1,DTItracts.tracts(1,:)+1,DTItracts.tracts(3,:)+1);

for fibnr = 1 : 1 : nFib
    
    if any(round(linspace(1,nFib,11))==fibnr)
        waitbar(fibnr/nFib,hwait)
    end
    
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
    
    % Calculate mean value along tract
    fa(fibnr)      = nanmean(FA(first:d:last));
    md(fibnr)      = nanmean(MD(first:d:last));
    lambda1(fibnr) = nanmean(L1(first:d:last));
    lambda2(fibnr) = nanmean(L2(first:d:last));
    lambda3(fibnr) = nanmean(L3(first:d:last));
end
t_elapsed = toc;
fprintf('It took %.2f seconds to calculate DTI indices.\n',t_elapsed)

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



