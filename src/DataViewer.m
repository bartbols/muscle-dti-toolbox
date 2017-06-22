function DataViewer( img,options )
%DATAVIEWER Interactive viewer for browsing through 3D image data and
%vector fields.
%
% Bart Bolsterlee, NeuRA
% April 2017
%
% INPUT:
% -img : nifti structure with image data (loaded with load_untouch_nii)
% -options -  optional input structure with one or more of the following fields:
%
%    The default values are applied if the field does not exist in the
%    structure 'options', or if the field is empty.
%
%   'Overlay'           - nifti structure to overlay on img using imshowpair
%   'Stack'             - index of the stack to be plotted if img is 4D 
%                         (for instance,to plot the B0 map of DTI data,
%                          use 'stack',1)
%
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
%    'VectorField2'     - nifti structure with 2nd vectorfield to overlay
%    'sF2'              - sampling fraction of vectors to plot, i.e. plots 
%                         every sF-th voxel in the rectangular grid.
%                         Default = 5
%    'AutoScaleVector2' - 'on' / 'off': if on, vectors in 2nd field are
%                         automatically scaled. Default = 'off'
%    'VectorColor2'     - color of the vectors in the 2nd vector field.
%                         Default = 'b'
%    'Mask2'            - binary mask of values in 2nd vector field to display.


% % If no options are provided, create empty structure 'options'
if nargin == 1
    options = struct();
end

if ~isfield(options,'Stack')
    options.Stack = [];
end

% Open a figure window
fig_dim = [50 50 500 600];
figure('Name',' Data viewer',...
    'NumberTitle','Off',...
    'Units','Pixels',...
    'Position',fig_dim);

% set some spacings and define the plot objects
spacing_small = 10;
spacing_left  = 60;% px
spacing_right = 10;
spacing_top   = 10;
spacing_bottom= 60;

% slider bar for slice number selection
xPos = spacing_left;
yPos = spacing_small;
slider_bar_height = 25;
axes_width = fig_dim(3)-spacing_left-spacing_right;
handle.sliceNumber = uicontrol(  'Style','slider',...
            'Min',1','Max',10,...
            'Value',5,...
            'Units','Pixels',...
            'Position',[xPos yPos axes_width slider_bar_height],...
            'CallBack', {@setSliceNumber},...
            'Tag','changeSlice');
       
% drop-down menu for slice orientation selection
popup_width = 100;
xPos = fig_dim(3) - spacing_right - popup_width;
yPos = yPos + spacing_small + slider_bar_height;
handle.sliceOrientation = uicontrol('Style','popup',...
    'String',{'Axial','Sagittal','Coronal'},...
    'Units','Pixels',...
    'Position',[xPos yPos popup_width slider_bar_height],...
    'CallBack',{@setOrientation},...
    'Tag','sliceOrientation');

text_width = 50;
xPos = xPos - spacing_small - text_width;

% text for slice orientation selection
handle.textOrientation = uicontrol('Style','text',...
    'String','Select orientation',...
    'Units','Pixels',...
    'Position',[xPos yPos text_width slider_bar_height],...
    'HorizontalAlignment','Right');
xPos_t = xPos;

if ~isfield(options,'Stack') && ndims(img.img) > 3
    options.Stack = [];
end

if ndims(img.img) > 3 && ~strcmpi(options.Stack,'rgb') && isempty(options.Stack)
    % text for stack selection
    options.Stack = 1;
    xPos = spacing_left;
    handle.textStack = uicontrol('Style','text',...
        'String','Select stack',...
        'Units','Pixels',...
        'Position',[xPos yPos text_width slider_bar_height],...
        'HorizontalAlignment','Left');

    % slider bar for stack selection
    xPos = xPos + text_width - spacing_small;
    width = xPos_t - xPos - spacing_small;

    handle.sliceNumber = uicontrol(  'Style','slider',...
                'Min',1','Max',size(img.img,4),...
                'Value',1,...
                'Units','Pixels',...
                'Position',[xPos yPos width slider_bar_height],...
                'CallBack', {@setStackNumber},...
                'Tag','changeStack',...
                'SliderStep',[1/(size(img.img,4)-1) 1/(size(img.img,4)-1)]);
end
% axes
xPos = spacing_left;
yPos = yPos + spacing_bottom + slider_bar_height;
axes_height = fig_dim(4) - yPos - spacing_top;

handle.axes = axes('Units','Pixels',...
    'Position',[xPos yPos axes_width axes_height],...
    'NextPlot','replace');

% Set all units to normalized so that resizing the figure window scales
% everything proportionally.
set(findall(gcf,'Units','Pixels'),'Units','Normalized')

setOrientation(handle.sliceOrientation,0)
%% Subfunctions
    function setOrientation(src,event)

        % Get the slice orientation from the popup menu
        str = get(src,'String');
        idx = get(src,'Value');
        sliceOrientation = str{idx}(1);
        
        % Set slider maximum to the number of slices in the current orientation
        maxSlice = round(img.hdr.dime.dim(strfind('SCA',sliceOrientation)+1));
        sliceNumber = round(maxSlice/2);
        set(findobj(gcf,'Tag','changeSlice'),...
            'Max',maxSlice,...
            'SliderStep',[1/(maxSlice-1) 10/(maxSlice-1)],...
            'Value',sliceNumber)
        
        % update the slice orientation in userdata
        userdata = get(gcf,'Userdata');
        userdata{1} = sliceOrientation;
        userdata{3} = maxSlice;           % maximum number of slices    
        userdata{2} = round(maxSlice/2);  % current slice number
        set(gcf,'UserData',userdata);
        
        update_slice %(img,options)
        
    end

    function setSliceNumber(src,event)
        % Get current slice number from the slider bar
        sliceNumber = round(get(src,'Value'));
        
        % update the slice number in userdata
        userdata = get(gcf,'Userdata');
        userdata{2} = sliceNumber;        
        set(gcf,'UserData',userdata);
        
        update_slice %(img,options)
    end

    function setStackNumber(src,event)
        % Get current stack number from the slider bar
        stackNumber = round(get(src,'Value'));
        % Update the slice number in 'options'
        options.Stack = stackNumber;
        update_slice %(img,options)
    end
    function update_slice %(img,options)
        cla
        userdata = get(gcf,'Userdata');
        plot_slice(img,userdata{1},userdata{2},options)
        
        % Add slice number 
        txt{1} = sprintf('Slice %d / %d',userdata{2},userdata{3});
        if ~isempty(options.Stack) && ~strcmpi(options.Stack,'rgb')
            txt{2} = sprintf('Stack %d',options.Stack);
        end
        text(0.02,0.98,txt,...
            'Units','normalized',...
            'FontSize',14,...
            'FontWeight','bold',...
            'Color','y',...
            'Tag','SliceNumberTxt',...
            'HorizontalAlignment','left',...
            'VerticalAlignment','top')
    end
        
end

