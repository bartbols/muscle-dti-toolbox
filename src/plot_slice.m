function varargout = plot_slice( img,sliceOrientation,sliceNumber,varargin )
%PLOT_SLICE operates on a 3D image. It display an axial, coronal or
%sagittal slice.
%
% INPUTS:
% -img- is a structure containing a 3D image file in nii format (loaded with
%             load_untouch_nii)
% -sliceOrientation- indicates the slice plane. May be
%   1. axial ('A'),
%   2. sagittal ('S'), or
%   3. coronal ('C').
% -sliceNumber-, a positive integer <= the size of the slice dimension.
%
% Optional inputs as 'parameter',<value> pairs or as fields in a structure.
%   'Stack'             - index of the stack to be plotted if img is 4D 
%                         (for instance,to plot the B0 map of DTI data,
%                          use 'stack',1)
%
%   'Overlay'           - nifti structure to overlay on img using imshowpair
%   'OverlayStack'      - Stack of image to overlay (if overlay image is
%                         4D)
%
% To overlay a vectorfield, add : 
%    'VectorField'      - nifti structure with vectorfield to overlay
%    'sF'               - sampling fraction of vectors to plot, i.e. plots 
%                         every sF-th voxel in the rectangular grid.
%                         Default = 5
%    'AutoScaleVector'  - 'on' / 'off': if on, vectors in field are
%                         automatically scaled. Default = 'off'
%    'VectorColor'      - color of the vectors in the vector field.
%                         Default = 'r'
%    'Mask'             - binary mask of values in vector field to display.
%
% To overlay a second vectorfield, add: 
%    'VectorField2'     - nifti structure with 2nd vectorfield to overlay
%    'sF2'              - sampling fraction of vectors to plot, i.e. plots 
%                         every sF-th voxel in the rectangular grid.
%                         Default = 5
%    'AutoScaleVector2' - 'on' / 'off': if on, vectors in 2nd field are
%                         automatically scaled. Default = 'off'
%    'VectorColor2'     - color of the vectors in the 2nd vector field.
%                         Default = 'b'
%    'Mask2'            - binary mask of values in 2nd vector field to display.
%
% OUTPUT:
% -h- handle to the plot object(s)
%
% The image orientation is:
%  x points from right --> left
%  y points from anterior --> posterior
%  z points from inferior --> superior

%% read inputs
p = inputParser;
p.CaseSensitive = false;
addRequired(p,'img')
addRequired(p,'sliceOrientation',@ischar)
addRequired(p,'sliceNumber',@isnumeric)
addParameter(p,'Stack',[])
addParameter(p,'Overlay',[])
addParameter(p,'OverlayStack',[])

addParameter(p,'VectorField',[])
addParameter(p,'sF',5)
addParameter(p,'AutoScaleVector','off')
addParameter(p,'VectorColor','r')
addParameter(p,'Mask',[])

addParameter(p,'VectorField2',[])
addParameter(p,'sF2',5)
addParameter(p,'AutoScaleVector2','off')
addParameter(p,'VectorColor2','b')
addParameter(p,'Mask2',[])
parse(p,img,sliceOrientation,sliceNumber,varargin{:});

%% -------------- IMAGE -------------
% Get the slice data and spatial reference object for plotting
R = get_spatial_ref(img,sliceOrientation);

[slice_data,y_axis_orientation,label_x,label_y] = ...
    get_slice_data(img,sliceOrientation,sliceNumber,p.Results.Stack);

%% -------------- IMAGE OVERLAY -------------
% If Overlay is provided as input, show both images using imshowpair.
% Otherwise, just show image 'img'
if ~isempty(p.Results.Overlay)
    
    R2 = get_spatial_ref(p.Results.Overlay,sliceOrientation);
    
    % Get slice closest to the selected slice in img.
    idx = strfind('SCA',sliceOrientation)+1;
    sliceNumber2 = round( (sliceNumber-0.5) * img.hdr.dime.pixdim(idx) ./...
        p.Results.Overlay.hdr.dime.pixdim(idx)+0.5);
    
    [slice_data2] = ...
      get_slice_data(p.Results.Overlay,sliceOrientation,sliceNumber2,p.Results.OverlayStack);
  h(1) = imshowpair(slice_data,R,slice_data2,R2);
else
    h(1) = imshow(slice_data,R);
end

%% -------------- DEFORMATION FIELD -------------
% Add vector field if provided
if ~isempty(p.Results.VectorField)
    % Get slice closest to the selected slice in img.
    idx = strfind('SCA',sliceOrientation)+1;
    sliceNumber2 = round( (sliceNumber-0.5) * img.hdr.dime.pixdim(idx) ./...
        p.Results.VectorField.hdr.dime.pixdim(idx)+0.5);
    
    hold on
    h(end+1) = add_vector_field(p.Results.VectorField,...
        sliceOrientation,...
        sliceNumber2,...
        p.Results.sF,...
        p.Results.AutoScaleVector,...
        p.Results.VectorColor,...
        p.Results.Mask);
    
end

% Add second vector field if provided
if ~isempty(p.Results.VectorField2)
    % Get slice closest to the selected slice in img.
    idx = strfind('SCA',sliceOrientation)+1;
    sliceNumber2 = round( (sliceNumber-0.5) * img.hdr.dime.pixdim(idx) ./...
        p.Results.VectorField2.hdr.dime.pixdim(idx)+0.5);
    
    hold on
    h(end+1) = add_vector_field(p.Results.VectorField2,...
        sliceOrientation,...
        sliceNumber2,...
        p.Results.sF2,...
        p.Results.AutoScaleVector2,...
        p.Results.VectorColor2,...
        p.Results.Mask2);
end

%% ---------------- AXES PROPERTIES ----------------
% Add labels to axes
xlabel(label_x)
ylabel(label_y)

% Auto-scale the intensity of the image
set(gca,'CLimMode','auto','YDir',y_axis_orientation,...
    'XLim',R.XWorldLimits,...
    'YLim',R.YWorldLimits)
axis on % show the coordinates on the axes

%% SUBFUNCTIONS

    function [slice_data,y_axis_orientation,label_x,label_y] =...
               get_slice_data(imdata,sliceOrientation,sliceNumber,stack)
        
        if ndims(imdata.img) == 4 && isempty(stack)
            % if data is 4D and no stack is selected, use the first stack.
            stack = 1; 
        end
        if ~isempty(stack) && ~strcmpi(stack,'rgb')
                % select the requested stack from the 4D image
                imdata.img = imdata.img(:,:,:,stack);
        end
        switch sliceOrientation
            case 'A'
                % axial slice (x-axis in plot = RL, y-axis in plot = AP)
                if strcmpi(stack,'rgb')
                    slice_data = permute(squeeze(imdata.img(:,:,sliceNumber,1:3)),[2 1 3]);
                else
                    slice_data = imdata.img(:,:,sliceNumber)';
                end
                    
                y_axis_orientation = 'reverse';
                label_x = '<-- right  X  left -->';
                label_y = '<-- posterior  Y  anterior -->';
            case 'S'
                % sagittal slice (x-axis in plot = AP, y-axis in plot = IS)
                y_axis_orientation = 'normal';
                if strcmpi(stack,'rgb')
                    slice_data = permute(squeeze(imdata.img(sliceNumber,:,:,1:3)),[2 1 3]);
                else
                    slice_data = squeeze(imdata.img(sliceNumber,:,:))';
                end
                label_x = '<-- anterior  Y  posterior -->';
                label_y = '<-- inferior  Z  superior -->';
            case 'C'
                % coronal slice
                if strcmpi(stack,'rgb')
                    slice_data = permute(squeeze(imdata.img(:,sliceNumber,:,1:3)),[2 1 3]);
                else
                    slice_data = squeeze(imdata.img(:,sliceNumber,:))';
                end
                y_axis_orientation = 'normal';
                label_x = '<-- right  X   left -->';
                label_y = '<-- inferior  Z  superior -->';
            otherwise
                error('Unknown slice orientation ''%s''. Permitted values are ''A'',''S'' and ''C''',sliceOrientation)
                
        end
    end % of subfunction get_slice_data

    function handle = add_vector_field(D,sliceOrientation,sliceNumber,...
        sF,AutoScaleVector,VectorColor,mask)
        % Get components of the deformation field in the requested slice
        % The vectors point FROM the fixed image TO the moving image.
        % The 1st, 2nd and 3rd value in the 5th dimension of the 
        % deformation image represent the x, y and z component.
        D.img = squeeze(D.img); % remove singleton dimension, if there is one.
        
        % Set vectors outside the mask to 0 (if a mask is provided).
        if ~isempty(mask)
            if isstruct(mask)
                % If the mask is provided as a nifti structure, make binary
                % mask of all non-zero values
                mask = mask.img ~= 0;
            end
            % Set values outside the mask to 0
            D.img = D.img .* cast(repmat(mask~=0,1,1,1,size(D.img,4)),'like',D.img);
        end
        
        switch sliceOrientation
            case 'A'
                % axial slice
                V1 = squeeze(D.img(:,:,sliceNumber,1));
                V2 = squeeze(D.img(:,:,sliceNumber,2));
                V3 = squeeze(D.img(:,:,sliceNumber,3));
                
                % Make grid of x and y coordinates of voxel centres
                [X,Y] = ndgrid( ((1:D.hdr.dime.dim(2))-0.5) .* D.hdr.dime.pixdim(2),...
                      ((1:D.hdr.dime.dim(3))-0.5) .* D.hdr.dime.pixdim(3) );
            case 'S'
                % sagittal slice
                V1 = squeeze(D.img(sliceNumber,:,:,2))';
                V2 = squeeze(D.img(sliceNumber,:,:,3))'; 
                V3 = squeeze(D.img(sliceNumber,:,:,1))';
                
                % Make grid of x and y coordinates of voxel centres
                [Y,X] = ndgrid( ((1:D.hdr.dime.dim(4))-0.5) .* D.hdr.dime.pixdim(4),...
                      ((1:D.hdr.dime.dim(3))-0.5) .* D.hdr.dime.pixdim(3));
            case 'C'
                % coronal slice
                V1 = squeeze(D.img(:,sliceNumber,:,1))';
                V2 = squeeze(D.img(:,sliceNumber,:,3))'; 
                V3 = squeeze(D.img(:,sliceNumber,:,2))'; 
                
                % Make grid of x and y coordinates of voxel centres
                [Y,X] = ndgrid( ((1:D.hdr.dime.dim(4))-0.5) .* D.hdr.dime.pixdim(4),...
                      ((1:D.hdr.dime.dim(2))-0.5) .* D.hdr.dime.pixdim(2));
                
        end
        Z = zeros(size(X));
        
        % subsample a fraction of the vectors
        Xs =   X(1:sF:end,1:sF:end);
        Ys =   Y(1:sF:end,1:sF:end);
        Zs =   Z(1:sF:end,1:sF:end);
        V1s = V1(1:sF:end,1:sF:end);
        V2s = V2(1:sF:end,1:sF:end);
        V3s = V3(1:sF:end,1:sF:end);
        
        % Plot as quivers
        handle = quiver3(Xs(:),Ys(:),Zs(:),V1s(:),V2s(:),V3s(:),...
                'AutoScale',AutoScaleVector,...
                'Color',VectorColor,'LineWidth',1);
    end

    function R = get_spatial_ref(imdata,sliceOrientation)
        % Create spatial reference object and get the data for the selected slice.
        % See 'doc imref2d' for more information about spatial referencing  
       switch sliceOrientation
           case 'A'
                R = imref2d(imdata.hdr.dime.dim([3 2]),... % image dimensions
                    [0 imdata.hdr.dime.dim(2) * imdata.hdr.dime.pixdim(2)],... % bounding box in x-direction (RL)
                    [0 imdata.hdr.dime.dim(3) * imdata.hdr.dime.pixdim(3)]);   % bounding box in y-direction (AP)
            case 'S'
                % sagittal slice (x-axis in plot = AP, y-axis in plot = IS)
                R = imref2d(imdata.hdr.dime.dim([4 2]),... % image dimensions
                    [0 imdata.hdr.dime.dim(3) * imdata.hdr.dime.pixdim(3)],... % bounding box in y-direction (AP)
                    [0 imdata.hdr.dime.dim(4) * imdata.hdr.dime.pixdim(4)]);   % bounding box in z-direction (IS)
            case 'C'
                % coronal slice
                R = imref2d(imdata.hdr.dime.dim([4 3]),... % image dimensions
                    [0 imdata.hdr.dime.dim(2) * imdata.hdr.dime.pixdim(2)],... % bounding box in x-direction (RL)
                    [0 imdata.hdr.dime.dim(4) * imdata.hdr.dime.pixdim(4)]);   % bounding box in z-direction (IS)
       end


    end
if nargout > 0
    varargout{1} = h;
end

end % of function

