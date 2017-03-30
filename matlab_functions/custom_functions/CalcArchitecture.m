function [ DTItracts ] = CalcArchitecture( DTItracts,surf_model,varargin)
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
% - surf_model: a structure array containing the fields
%               'vertices','faces' with the vertices and faces of the 
%               surface used for extrapolating the tracts.
%
% Optional: 
% - order     : order of the polynomial fit. Default = 3
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
addRequired(p,'DTItracts',@isstruct)
addRequired(p,'surf_model',@isstruct)
addOptional(p,'order',3,@(x)(isnumeric(x) && mod(x,1)==0 && x>0))
parse(p,DTItracts,surf_model,varargin{:})
order = p.Results.order;

%%
% Truncate tracts
DTItracts   = TruncateTracts(DTItracts,surf_model);

% Extrapolate tracts
DTItracts   = ExtrapolateTracts(DTItracts,surf_model,order);

% Calculate pennation angle
DTItracts.penangle = CalcPenAngle( DTItracts,surf_model);
    
% Calculate curvature
DTItracts.curvature = CalcCurvature( DTItracts.PolyCoeff);



end

