function [ controlPoints3,controlPoints2 ] = segment_slice( controlPoints1,img1,img2,varargin )
%SEGMENT_SLICE predicts the location of control points defined on image 1
%to image 2.

% Read inputs and define default if parameters are not provided
p = inputParser;

% -- Required inputs
addRequired(p,'controlPoints1',@(x) isnumeric(x) && size(x,2) == 2)

addRequired(p,'img1',@(x) validateattributes(x,{'numeric'},...
    {'nonempty','nonsparse'},...
    'segment_slice','img1'))

addRequired(p,'img2',@(x) validateattributes(x,{'numeric'},...
    {'nonempty','nonsparse'},...
    'segment_slice','img2'))

% -- Optional inputs
addParameter(p,'regionSize',[15 15],@(x) validateattributes(x,{'numeric'},...
    {'positive','nonempty','nonsparse','numel',2,'vector'},...
    'segment_slice','regionSize'))

addParameter(p,'blockSize',[11 11],@(x) validateattributes(x,{'numeric'},...
    {'positive', 'nonempty', 'nonsparse', 'numel', 2, 'integer', ...
    'vector',  'odd', '>=', 5},...
    'segment_slice','blockSize'))

addParameter(p,'r1',2,   @(x) validateattributes(x,{'numeric'},...
    {'scalar','positive'},'segment_slice','r1'))

addParameter(p,'r2',4,    @(x) validateattributes(x,{'numeric'},...
    {'scalar','positive'},'segment_slice','r2'))

addParameter(p,'dr',0.05, @(x) validateattributes(x,{'numeric'},...
    {'scalar','positive'},'segment_slice','dr'))

addParameter(p,'opt_metric','lsq',@(x) any(strcmp(x,{'lsq','corr'})))

addParameter(p,'channels',[1 3 4],@(x) validateattributes(x,{'numeric'},...
    {'vector','positive','numel',3},'segment_slice','channels'))

parse(p,controlPoints1,img1,img2,varargin{:});

% Check if 2D or 3D data is provided.
img1 = double(img1);
img2 = double(img2);
switch ndims(img1)
    case 2
        % Image has only 1 channel: convert to RGB map with 3 channels.
        img1_tracking = ind2rgb(img1,gray(max(round(img1(:)))));
        img2_tracking = ind2rgb(img2,gray(max(round(img2(:)))));
    case 3
        % Multiple channels in data.
        % Only 3 channels can be used for KLT tracking. These channels are
        % provided in the optional input argument 'channels'. If not provided,
        % channel 1, 3 and 4 are used, because these are the channels in mdixon
        % data with the highest contrast between muscle and surrounding
        % tissues.
        img1_tracking = zeros(size(img1,1),size(img1,2),3);
        for i = 1 : 3
            C = img1(:,:,p.Results.channels(i));
            img1_tracking(:,:,i) = C  / max(C(:));
        end
        
        img2_tracking = zeros(size(img2,1),size(img2,2),3);
        for i = 1 : 3
            C = img2(:,:,p.Results.channels(i));
            img2_tracking(:,:,i) = C  / max(C(:));
        end
end

% Step 1: track control points from image 1 to image 2 using
%         Kanade-Lucas-Tomasi (KLT) feature-tracking algorithm.
%         Tracking is done for each control point independently in a region
%         around each control point defined by regionSize.

controlPoints2 = zeros(size(controlPoints1));

rng(0); % Reset random number generator for reproducible results.
for ii = 1 : (size(controlPoints1,1)-1)
    
    % Make regular grid of control points
    [Xg,Yg] = ndgrid( linspace(-p.Results.regionSize(1)/2,p.Results.regionSize(1)/2,21) + controlPoints1(ii,1),...
        linspace(-p.Results.regionSize(2)/2,p.Results.regionSize(2)/2,21) + controlPoints1(ii,2));
    pt1 = [Xg(:) Yg(:)];
    
    % Create a tracker object.
    tracker = vision.PointTracker('MaxBidirectionalError',1,...
        'BlockSize',p.Results.blockSize);
    
    % Initialize the tracker.
    initialize(tracker,pt1,img1_tracking);
    [pt2, isFound] = step(tracker,img2_tracking);
    
    % Estimate the geometric transformation between the old points
    % and the new points and eliminate outliers.
    xform = estimateGeometricTransform(...
        pt1(isFound,:), pt2(isFound,:), 'similarity', 'MaxDistance', 1);
    
    % Estimate location of centre point.
    C_est = [controlPoints1(ii,:) 1] * xform.T;
    controlPoints2(ii,:) = C_est(1:2);
    
end

controlPoints2(end,:) = controlPoints2(1,:); % last point is the same as the first point
% Step 2: Optimize location using intensity profile matching.

[controlPoints3, dist_along_n,best_fit_metric] = ...
    match_intensity(controlPoints1,controlPoints2,...
    img1, img2, p.Results.r1,p.Results.r2,p.Results.dr,p.Results.opt_metric );

end