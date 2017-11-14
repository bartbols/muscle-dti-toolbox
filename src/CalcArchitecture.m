function [ DTItracts ] = CalcArchitecture( DTItracts,SurfModel,varargin)
%CALCARCHITECTURE truncates and extrapolates raw tracts (output of
%TrackFibres) and calculates fibre length, pennation angle and curvature.
%
% ----------------- USAGE ----------------- 
% DTItracts = CalcArchitecture(DTItracts,surf_model)
% or
% DTItracts = CalcArchitecture(DTItracts,surf_model,order)
%
% ----------------- INPUT -----------------
% Required:
% - DTItracts : a structure array containing at least the fields tracts_xyz 
%               and fibindex_trunc (created with TruncateTracts.m)
% - SurfModel : a structure array containing the fields
%               'vertices','faces' with the vertices and faces of the 
%               surface used for extrapolating the tracts.
%
% Optional inputs, provided as 'parameter',<value> pairs:
% - order     : order of the polynomial fit. Default = 3
% - aponeurosis: a structure array containing the fields
%               'vertices','faces' with the vertices and faces of the 
%               aponeurosis used for extrapolating the tracts.
%
% ----------------- OUTPUT -----------------
% - DTItracts : a structure array with reconstruction parameters and the
%               architectural measures added as fields.
%
% This function calls TruncateTracts, ExtrapolateTracts, CalcPenAngle and
% CalcCurvature. Type help <function_name> for more information on the
% outputs of each of these functions.

% Optional: 
% - order     : order of the polynomial fit. Default = 3

%% Check inputs
p = inputParser;
addRequired(p,'DTItracts',@(x) isstruct(x) || exist(x,'file')==2)
addRequired(p,'SurfModel',@(x) isstruct(x) || endsWith(x,'.stl','IgnoreCase',true))
addParameter(p,'order',3,@(x) isscalar(x) && x>0)
addParameter(p,'aponeurosis',[],@(x) isstruct(x) || endsWith(x,'.stl','IgnoreCase',true))
parse(p,DTItracts,SurfModel,varargin{:})

order       = p.Results.order;
aponeurosis = p.Results.aponeurosis;

%% Read inputs
% if tract filename is provided, read the file.
if ~isstruct(DTItracts)
    DTItracts = load(DTItracts);
end

% if surface model filename is provided, read the stl file
if ~isstruct(SurfModel)
    if exist(SurfModel,'file') == 2
        SurfModel = stlread(SurfModel);
    else
        error('%s does not exist.',SurfModel)
    end
end

% Read the aponeurosis surface, if provided
if ~isempty(aponeurosis) && ~isstruct(aponeurosis)
    if exist(aponeurosis,'file') == 2
        aponeurosis = stlread(aponeurosis);
    else
        error('%s does not exist.',aponeurosis)
    end
end
%%
% Truncate tracts
DTItracts   = TruncateTracts(DTItracts,SurfModel,...
    'aponeurosis',aponeurosis);

% Extrapolate tracts
DTItracts   = ExtrapolateTracts(DTItracts,SurfModel,...
    'order',order,...
    'aponeurosis',aponeurosis);

% Calculate pennation angle
DTItracts.penangle = CalcPenAngle( DTItracts,SurfModel,...
    'aponeurosis',aponeurosis);
    
% Calculate curvature
DTItracts.curvature = CalcCurvature( DTItracts.PolyCoeff);



end

