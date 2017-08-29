function muscle_segmenter(filename)
%% MUSCLE_SEGMENTER is a graphical user interface for semi-automated 
% segmentation of 3D/4D images.
%
% Bart Bolsterlee
% Neuroscience Research Australia
% 21/06/2017
% Version 1.0

% Define variables that will be shared across the main and subfunctions
handle = struct();
C=[]; % cell of accepted control points for each slice
Cs=[]; % interpolated points along the curve (for plotting)
Ns=[]; % normal vectors to the curve (at interpolated values)
t=[]; % t-values of control points
ts=[]; % t-values of interpolated points

idx=[];
channel = 1;
mask_alpha = 0.2;
slice_nr=[];slice_nr2=[];slice_prev1=[];slice_prev2=[];
slice_data=[];slice_data2=[];
pars = struct();
binaryMask = [];
data=struct();
voxelSize = [];
imdim = [];
mask_color = [0 255 0];
%%
% ---------------- MAKE AXES AND GUI OBJECTS --------------------------
fig_size = [900 600];
figure('Position',[50 150 fig_size],...
    'NumberTitle','off',...
    'Name','Muscle segmenter v1.0')

% ----------- Axis window and plot objects
spacing = 25;
ax_dim = [fig_size(2)-2*spacing fig_size(2)-2*spacing];
handle.axes = axes(...
    'Units','pixels',...
    'Position',[spacing spacing ax_dim]);
handle.img = imshow(NaN(2,2));
% handle.img = image(NaN(2,2));
hold on

width = 40;
height = 20;
uicontrol('Style','popup',...
    'Units','Pixels',...
    'Position',[spacing+ax_dim(1)-width-5 spacing+ax_dim(2)-height-5 width height],...
    'Tag','channel',...
    'Callback',@set_channel,...
    'String','',...
    'Visible','off')

handle.mask = image(NaN(2,2),'AlphaData',0);

handle.curvefit  = plot(NaN,NaN,'r','LineWidth',2,...
    'Tag','curvefit',...
    'PickableParts','none');
handle.cp  = plot(NaN,NaN,'gs','LineWidth',2,...
    'ButtonDownFcn',@move_controlPoint,...
    'Tag','controlPoints');
handle.slice_nr = text(0.05,0.95,'',...
    'FontWeight','bold',...
    'FontSize',12,...
    'Color','y',...
    'HorizontalAlignment','left',...
    'VerticalAlignment','top');
caxis auto

handle.slider = uicontrol('Style','slider',...
    'Min',1','Max',10,...
    'Value',5,...
    'Units','Pixels',...
    'Position',[spacing 0 ax_dim(1) spacing],...
    'Callback',@set_slicenumber);

% ----------- Side panels
back_color = 0.95*ones(1,3);

% -- File
spacing2 = 10;
panel_width  = fig_size(1) - ax_dim(1) -  3*spacing;
height       = 40;
panel_height = 3*height+spacing+3*spacing2;
xPos1        = 2*spacing + ax_dim(1);
yPos1        = fig_size(2) - panel_height - spacing2;

hp1 = uipanel('Title','File',...
    'FontSize',12,...
    'Backgroundcolor',back_color,...
    'Units','Pixels',...
    'Position',[xPos1 yPos1 panel_width panel_height]);

xPos   = spacing2;
width  = (panel_width - 3*spacing2)/2;
yPos   = panel_height-spacing-height;
uicontrol('Parent',hp1,...
    'Style','push',...
    'String','Load file',...
    'Units','pixels',...
    'BackGroundColor',back_color,...
    'Position',[xPos yPos width height],...
    'Callback',@load_data)

xPos = xPos + width + spacing2;
% Export mask button
handle.export = uicontrol('Parent',hp1,...
    'Style','push',...
    'String','Export NIFTI',...
    'Units','pixels',...
    'BackGroundColor',back_color,...
    'Position',[xPos yPos width height],...
    'Callback',@export_mask);

% Load segmentation button
yPos   = yPos-spacing2-height;
handle.load_segm = uicontrol('Parent',hp1,...
    'Style','push',...
    'String','Load segmentation',...
    'Units','pixels',...
    'BackGroundColor',back_color,...
    'Position',[xPos yPos width height],...
    'Callback',@load_segmentation);

% Save segmentation button
yPos   = yPos-spacing2-height;
handle.save_segm = uicontrol('Parent',hp1,...
    'Style','push',...
    'String','Save segmentation',...
    'Units','pixels',...
    'BackGroundColor',back_color,...
    'Position',[xPos yPos width height],...
    'Callback',@save_segmentation);

xPos   = spacing2;
height = 60;
yPos   = panel_height-spacing-2*height-spacing2;
% width = panel_width - 2*spacing2;
uicontrol('Parent',hp1,...
    'Style','text',...
    'String','<no file selected>',...
    'Units','pixels',...
    'BackGroundColor',back_color,...
    'Position',[xPos yPos width height],...
    'HorizontalAlignment','left',...
    'Tag','info_text')

% -- Segmentation
height = 40;
panel_height = 3*height + 3*spacing2 + spacing;
yPos2 = yPos1 - panel_height - spacing2;
hp2 = uipanel('Title','Segmentation',...
    'FontSize',12,...
    'Backgroundcolor',back_color,...
    'Units','Pixels',...
    'Position',[xPos1 yPos2 panel_width panel_height]);

xPos   = spacing2;
width  = (panel_width - 3*spacing2)/2;
yPos = panel_height-spacing-height;
uicontrol('Parent',hp2,...
    'Style','push',...
    'String','Set control points',...
    'Units','pixels',...
    'BackGroundColor',back_color,...
    'Position',[xPos yPos width height],...
    'Callback',@set_controlPoints)

xPos = xPos + spacing2 + width;
uicontrol('Parent',hp2,...
    'Style','push',...
    'String','Add control point',...
    'Units','pixels',...
    'BackGroundColor',back_color,...
    'Position',[xPos yPos width height],...
    'Callback',@add_controlPoints)

xPos   = spacing2;
yPos = yPos - height - spacing2;
handle.segment = uicontrol('Parent',hp2,...
    'Style','push',...
    'String','Segment',...
    'Units','pixels',...
    'BackGroundColor',back_color,...
    'Position',[xPos yPos width height],...
    'Callback',@segment,...
    'Tag','segment');

xPos = xPos + spacing2 + width;
width1 = width/2;
% width = 20;
handle.up_or_down = uicontrol('Parent',hp2,...
    'Style','popup',...
    'String',{'down','up'},...
    'Value',1,...
    'Units','pixels',...
    'BackGroundColor',back_color,...
    'Position',[xPos yPos width1 height]);

xPos = xPos + width1+spacing2/2;
width1 = 15;
uicontrol('Parent',hp2,...
    'Style','text',...
    'String','by',...
    'Units','pixels',...
    'BackGroundColor',back_color,...
    'Position',[xPos yPos width1 height],...
    'HorizontalAlignment','left')

xPos = xPos + width1+spacing2/2;
width1 = panel_width - xPos - spacing2;
handle.stepsize = uicontrol('Parent',hp2,...
    'Style','popup',...
    'String',{'1','2','3','4','5'},...
    'Units','pixels',...
    'BackGroundColor',back_color,...
    'Position',[xPos yPos width1 height],...
    'HorizontalAlignment','left');


xPos   = spacing2;
yPos = yPos - height - spacing2;
handle.binarize_slice = uicontrol('Parent',hp2,...
    'Style','push',...
    'String','Binarize slice',...
    'Units','pixels',...
    'BackGroundColor',back_color,...
    'Position',[xPos yPos width height],...
    'Callback',@binarize,...
    'Tag','binarize_slice');

xPos = xPos + spacing2 + width;
uicontrol('Parent',hp2,...
    'Style','push',...
    'String','Binarize all',...
    'Units','pixels',...
    'BackGroundColor',back_color,...
    'Position',[xPos yPos width height],...
    'Callback',@binarize,...
    'Tag','binarize_all')

% -- Parameters
height = 20;
nPars = 3;
panel_height = nPars*height + nPars*spacing2 + spacing;
yPos3 = yPos2 - panel_height - spacing2;
hp3 = uipanel('Title','Parameters',...
    'FontSize',12,...
    'Backgroundcolor',back_color,...
    'Units','Pixels',...
    'Position',[xPos1 yPos3 panel_width panel_height]);

% text - edit - text
width1 = 80;
width2 = 30;
width3 = 40;
parInfo = {'Block size','blockSize','voxels';...
    'Region size','regionSize','voxels';...
    'Profile length','r1','voxels'};
for i = 1 : nPars
    if i == 1
        yPos = panel_height-spacing-height;
    else
        yPos = yPos - height-spacing2;
    end
    xPos = spacing2;
    % text box - parameter name
    uicontrol('Parent',hp3,...
        'Style','text',...
        'String',[parInfo{i,1} ':'],...
        'HorizontalAlignment','Left',...
        'Units','Pixels',...
        'Position',[xPos yPos width1 height])
    
    % edit box - parameter value
    xPos = xPos + width1 + spacing2;
    uicontrol('Parent',hp3,...
        'Style','edit',...
        'String',NaN,...
        'Tag',parInfo{i,2},...
        'Units','Pixels',...
        'Position',[xPos yPos width2 height],...
        'Callback',@set_parameters)
    
    % text box - parameter unit
    xPos = xPos + width2 + spacing2;
    uicontrol('Parent',hp3,...
        'Style','text',...
        'String',parInfo{i,3},...
        'HorizontalAlignment','Left',...
        'Units','Pixels',...
        'Position',[xPos yPos width3 height])
end

xPos = xPos + spacing2 + width3;
width = panel_width - xPos - spacing2;
yPos = spacing2;
height = panel_height - spacing - spacing2;
uicontrol('Parent',hp3,...
    'Style','push',...
    'String','Update',...
    'Units','pixels',...
    'BackGroundColor',back_color,...
    'Position',[xPos yPos width height],...
    'Callback',@segment,...
    'Tag','segment_update')

% -- Contrast
panel_height = yPos3-spacing;
yPos4 = yPos3 - panel_height - spacing2;
hp4 = uipanel('Title','Contrast',...
    'FontSize',12,...
    'Backgroundcolor',back_color,...
    'Units','Pixels',...
    'Position',[xPos1 yPos4 panel_width panel_height]);

xPos = spacing2;
height = 20;
width1  = (panel_width-6*spacing2)/5;
% width1  = 20;
% text box - min contrast
w = 0.6*width1;
uicontrol('Parent',hp4,...
    'Style','text',...
    'String','Min',...
    'HorizontalAlignment','Left',...
    'Units','Pixels',...
    'Position',[xPos yPos w height])

% edit box - min contrast value
xPos = xPos + w + spacing2;
w = 1.4*width1;
handle.MinContrast = uicontrol('Parent',hp4,...
    'Style','edit',...
    'String',NaN,...
    'Tag','MinContrast',...
    'Units','Pixels',...
    'Position',[xPos yPos w height],...
    'Callback',@set_contrast);

% text box - max contrast
xPos = xPos + w + spacing2;
w = 0.6*width1;
uicontrol('Parent',hp4,...
    'Style','text',...
    'String','Max',...
    'HorizontalAlignment','Left',...
    'Units','Pixels',...
    'Position',[xPos yPos w height])

% edit box - min contrast value
xPos = xPos + w + spacing2;
w = 1.4*width1;
handle.MaxContrast = uicontrol('Parent',hp4,...
    'Style','edit',...
    'String',NaN,...
    'Tag','MaxContrast',...
    'Units','Pixels',...
    'Position',[xPos yPos w height],...
    'Callback',@set_contrast);

% push button - auto contrast
xPos = xPos + w + spacing2;
uicontrol('Parent',hp4,...
    'Style','push',...
    'String','Auto',...
    'Tag','AutoContrast',...
    'Units','Pixels',...
    'Position',[xPos yPos width1 height],...
    'Callback',@set_contrast)


set(findall(gcf,'Units','pixels'),'Units','normalized')


% Read data, if provided
if nargin > 0
    load_data(filename)
end

%%
% --------------------------------------------------------------------
% ----------------------- SUBFUNCTIONS -------------------------------
% --------------------------------------------------------------------

    function load_data(input1,events)
        % Loads a new image and updates the axes and segmentation variables
        
        if isobject(input1)
            % Load-button was clicked. Open a dialog box to select a file.
            [fname, pathname] = uigetfile( ...
                {'*.nii;*.nii.gz;',...
                'Image files (*.nii,*.nii.gz)'}, ...
                'Select an image');
            if isnumeric(fname)
                % cancel button was clicked
                return
            end
            filename = fullfile(pathname,fname);
            data = load_untouch_nii(filename);
        else
            % input is filename of nifti structure. Load the data
            filename = input1;
            data = load_untouch_nii(input1);
        end
        
        if ndims(data.img) == 5
            % 5D data - remove singleton dimension
            data.img = squeeze(data.img);
            data.hdr.dime.dim(1) = 4;
            data.hdr.dime.dim(5) = 4;
            data.hdr.dime.dim(6) = 1;
        end
        
        [pathstr,fname,ext] = fileparts(filename);
        voxelSize = data.hdr.dime.pixdim(2:4);
        imdim     = size(data.img);
        
        if length(imdim) == 4
            % Activate drop-down menu for channel selection
            channel = 1;
            set(findobj(gcf,'Tag','channel'),...
                'Visible','on',...
                'String',num2cell(1:imdim(4)),...
                'Value',1)
        elseif length(imdim) == 3
            % Deactiviate drop-down menu.
            channel  = [];
            set(findobj(gcf,'Tag','channel'),...
                'visible','off',...
                'String','')
        end
        
        % Automatically define parameters
        pars.regionSize = max([11,round(10 ./ voxelSize(1))]);
        pars.blockSize  = pars.regionSize;
        pars.r1 = 5 / voxelSize(1);
        pars.r2 = 2 * pars.r1;
        pars.dr = pars.r1 / 51;
        pars.stepsize = max([1,floor(3 / voxelSize(3))]);
        set_parameters
        
        % Update info text
        dimstr = sprintf('%dx%dx%d',imdim(1:3));
        if length(imdim)>=4
            dimstr = [dimstr 'x' int2str(imdim(4))];
        end
        set(findobj(gcf,'Tag','info_text'),...
            'String',{...
            [fname ext],...
            sprintf('%.2fx%.2fx%.2fmm',voxelSize),...
            dimstr})
        
        % Show middle slice
        slice_nr = round(imdim(3)/2);
        slice_nr2 = slice_nr;
        
        % Set slider dimensions
        set(handle.slider,...
            'Max',imdim(3),...
            'SliderStep',[1/(imdim(3)-1) 10/(imdim(3)-1)],...
            'Value',slice_nr)
        
%         if ndims(data.img) == 4
%             slice_data = data.img(:,:,slice_nr,channel);
%         else
%             slice_data = data.img(:,:,slice_nr);
%         end
%         set(handle.img,'CData',slice_data)
%         set(handle.slice_nr,'String',sprintf('%d / %d',slice_nr,imdim(3)))
        
        % Make binary mask
        binaryMask = false(imdim(1:3));
        
        % Make empty mask for plotting
        M = zeros([imdim(1) imdim(2) 3]);
        for c = 1 : 3
            M(:,:,c) = mask_color(c);
        end
        set(handle.mask,'AlphaData',binaryMask(:,:,slice_nr),...
            'CData',M)
        
        % Make empty matrix with control points C
        C = cell(imdim(3),1);
        set([handle.curvefit handle.cp],...
            'XData',NaN,'YData',NaN);
        axis equal auto
        
        % update plot
        set_slicenumber
        
    end

    function set_channel(src,events)
        channel = get(src,'Value');
        slice_data = data.img(:,:,slice_nr,channel);
        set(handle.img,'CData',slice_data)
        
    end

    function set_slicenumber(varargin)
        % Get current slice number from the slider bar
        slice_nr = round(get(handle.slider,'Value'));
        
        if ndims(data.img) == 4
            slice_data = data.img(:,:,slice_nr,channel);
        else
            slice_data = data.img(:,:,slice_nr);
        end
        
        set(handle.mask,'alphaData',binaryMask(:,:,slice_nr)*mask_alpha)
        
        set(handle.img,'CData',slice_data)
        
        update_curve(C{slice_nr});
        set(handle.slice_nr,'String',sprintf('%d / %d',slice_nr,imdim(3)))
        
        if isempty(C{slice_nr})
            set([handle.binarize_slice,handle.segment],...
                'Enable','off')
        else            
            set([handle.binarize_slice,handle.segment],...
                'Enable','on')
        end
        
        clim = get(handle.axes,'CLim');
        set(handle.MinContrast,...
            'String',sprintf('%d',clim(1)));
        set(handle.MaxContrast,...
            'String',sprintf('%d',clim(2)));
        
    end

    function segment(src,events)
        % slice_nr      : slice number of reference slice
        % slice_nr2     : slice number of slice to be predicted
        pars.stepsize = get(handle.stepsize,'Value');
        up_down = get(handle.up_or_down,'Value');
        if strcmp(get(src,'Tag'),'segment') && up_down == 2
            % segment up
            slice_nr2 = slice_nr + pars.stepsize;
            if slice_nr2 > imdim(3)
                return
            end
        elseif strcmp(get(src,'Tag'),'segment') && up_down == 1
            slice_nr2 = slice_nr - pars.stepsize;
            if slice_nr2 <= 0
                return
            end
        elseif strcmp(get(src,'Tag'),'segment_update')
            % Use slices used for previous segmentation
            slice_nr  = slice_prev1;
            slice_nr2 = slice_prev2;
        end
        
        if ndims(data.img) == 4
            slice_data  = data.img(:,:,slice_nr,channel);
            slice_data2 = data.img(:,:,slice_nr2,channel);
        else
            slice_data  = data.img(:,:,slice_nr);
            slice_data2 = data.img(:,:,slice_nr2);
        end
        
        % Predict control point location in the selected slice.
        C{slice_nr2} = segment_slice(C{slice_nr},slice_data,slice_data2,...
            'regionSize',[1 1]*pars.regionSize,...
            'blockSize',[1 1]*pars.blockSize,...
            'r1',pars.r1,...
            'r2',pars.r2,...
            'dr',pars.dr);
        
        % Remove mask from plot
        binaryMask(:,:,slice_nr2) = false;
        set(handle.mask,'alphaData',binaryMask(:,:,slice_nr2)*mask_alpha);
        
        % Update slice numbers
        slice_prev1 = slice_nr;
        slice_prev2 = slice_nr2;
        slice_nr   = slice_nr2; % currently displayed
        
        % Update plot
        set(handle.img,'CData',slice_data2);
        update_curve(C{slice_nr});
        set(handle.slice_nr,'String',sprintf('%d / %d',slice_nr,imdim(3)))
        
        % Update slider bar
        set(handle.slider,'Value',slice_nr)
        
    end

%     function accept(src,events)
%     end

    function binarize(src,events)
        
        switch get(src,'Tag')
            case 'binarize_slice'
                % make binary mask of current slice only by setting voxels
                % within the fitted curve to 1
                binaryMask(:,:,slice_nr) = poly2mask(Cs(:,1), Cs(:,2), imdim(1), imdim(2));
                set(handle.mask,'alphaData',binaryMask(:,:,slice_nr)*mask_alpha);
            case 'binarize_all'
                % make binary mask of all slices with control points
                for k = 1 : imdim(3)
                    if isempty(C{k});continue;end
                    Ck = fit_closed_curve(C{k});
                    binaryMask(:,:,k) = poly2mask(Ck(:,1), Ck(:,2), imdim(1), imdim(2));
                end
                
        end
        set(handle.mask,'alphaData',binaryMask(:,:,slice_nr)*mask_alpha);
    end

    function export_mask(varargin)
        % Write current mask to a nifti file. Use same header information
        % as currently loaded image.
        mask_nii = data;
        
        mask_nii.img = uint16(binaryMask);
        mask_nii.hdr.dime.dim(1) = 3;
        mask_nii.hdr.dime.dim(5) = 1;
        mask_nii.hdr.dime.bitpix = 16;
        mask_nii.hdr.dime.scl_slope = 1;
        mask_nii.hdr.dime.scl_slope = 1;
        mask_nii.hdr.dime.glmax = max(binaryMask(:));
        mask_nii.hdr.dime.datatype = 16;
        
        [pathstr,fname,ext] = fileparts(filename);
        [~,fname,~] = fileparts(fname);
        mask_save_path = uigetdir(pathstr,'Select folder to save mask to');
        if mask_save_path == 0
            return
        end
        
        mask_filename = inputdlg('Enter filename (including extension .nii.gz)',...
            'Filename',1,{[fname '_mask.nii.gz']});
        mask_full_filename = fullfile(char(mask_save_path),char(mask_filename));
        save_untouch_nii(mask_nii,mask_full_filename)
        fprintf('Mask exported to %s\n',mask_full_filename)
    end

    function save_segmentation(varargin)
        % Write current mask and control points to a MAT file.
        
        [pathstr,fname] = fileparts(filename);
        [~,fname,~] = fileparts(fname);
        segm_save_path = uigetdir(pathstr,'Select folder to save mask to');
        if segm_save_path == 0
            return
        end
        
        segm_filename = inputdlg('Enter filename (including extension .mat)',...
            'Filename',1,{[fname '_segm.mat']});
        segm_full_filename = fullfile(char(segm_save_path),char(segm_filename));
        save(segm_full_filename,'binaryMask','C')
        fprintf('Segmentation saved as %s\n',segm_full_filename)
    end

    function load_segmentation(varargin)
        % Load current mask and control points from a MAT file.        
        pathstr = fileparts(filename);
        [segm_filename,segm_path] = uigetfile(pathstr,'Select segmentation file');
        if segm_path == 0
            return
        end
        load(fullfile(segm_path,segm_filename))
        set_slicenumber % update plot
    end


    function set_controlPoints(varargin)
        
        % reset existing control points on current slice
        update_curve([])
        h = impoly(handle.axes);
        wait(h);
        C{slice_nr} = getPosition(h);
        C{slice_nr} = [C{slice_nr};C{slice_nr}(1,:)];
        delete(h)
        
        % Update curve and plot
        update_curve        
        set_slicenumber
    end

    function add_controlPoints(src,events)
        handle.point = impoint(handle.axes);
        setColor(handle.point,'b');
        wait(handle.point);
        
        % find closest location
        d = sqrt(sum((Cs - ones(size(Cs,1),1) * getPosition(handle.point)).^2,2));
        [mindist,id] = min(d);
        
        % Find control points before and after added point
        idx1 = find(diff(sign( ts(id) - t))~=0);
        idx2 = idx1+1;
        
        % Add point between idx1 and idx2
        C{slice_nr} = [C{slice_nr}(1:idx1,:);...
            getPosition(handle.point);...
            C{slice_nr}(idx2:end,:)];
        
        delete(handle.point);
        update_curve
        set_slicenumber
    end

    function move_controlPoint(~,events)
        %MOVE_POINT allows the user to interactively drag a point to a new location
        
        % find point nearest to click
        [~,idx] = min(sqrt(sum((C{slice_nr} - ones(size(C{slice_nr},1),1) * events.IntersectionPoint(1:2)).^2,2)));
        
        % drag to new position
        handle.point = impoint(gca,C{slice_nr}(idx,:));
        setColor(handle.point,'b');
        addNewPositionCallback(handle.point,@update_controlPoints);
        wait(handle.point);
        
        try
            handle.point.Deletable;
            delete(handle.point)
        catch
            % Point was deleted. Update curve.
            if size(C{slice_nr},1) <= 4
                % Minimum number of points is 3.
                warning('Minimum number of control points is 3.')
                return
            end
            if idx == 1 || idx == size(C{slice_nr},1)
                % remove first and last point
                C{slice_nr}([1 end],:) = [];
                C{slice_nr} = [C{slice_nr};C{slice_nr}(1,:)];
            else
                C{slice_nr}(idx,:) = [];
            end
            update_curve
            set_slicenumber
        end
    end

    function update_controlPoints(pos)
        if idx == 1 || idx == size(C{slice_nr},1)
            C{slice_nr}([1 size(C{slice_nr},1)],:) = [1;1]*pos;
        else
            C{slice_nr}(idx,:) = pos;
        end
        update_curve
        set_slicenumber
    end

    function update_curve(varargin)
        if nargin == 1
            C{slice_nr} = varargin{1};
        end
        
        if isempty(C{slice_nr})
            set(handle.cp,...
                'XData',NaN,...
                'YData',NaN)
            set(handle.curvefit,...
                'XData',NaN,...
                'YData',NaN)
        else
            [Cs,Ns,ts,t] = fit_closed_curve(C{slice_nr});
            set(handle.curvefit,...
                'XData',Cs(:,1)',...
                'YData',Cs(:,2)')
            
            set(handle.cp,...
                'XData',C{slice_nr}(:,1)',...
                'YData',C{slice_nr}(:,2)')
        end
    end

    function set_parameters(src,events)
        if nargin > 0
            parname = get(src,'Tag');
            value = cell2mat(textscan(get(src,'String'),'%f',1));
            switch parname
                case 'regionSize'
                    if value < 11
                        warning('Minimum region size is 11')
                        pars.regionSize = 11;
                    else
                        pars.regionSize = value;
                    end
                case 'blockSize'
                    if value < 11
                        warning('Minimum block size is 11')
                        pars.blockSize = 11;
                    else
                        pars.blockSize = value;
                    end
                case 'r1'
                    if value <= 0
                        warning('profile length has to be positive')
                        pars.r1 = 5 / voxelSize(1);
                    else
                        pars.r1 = value;
                    end
                    % Also update r2 and dr
                    pars.r2 = 2 * pars.r1;
                    pars.dr = pars.r1 / 51;
                case 'stepsize'
                    if value <= 0
                        warning('step size has to be positive')
                    else
                        pars.stepsize = round(value);
                    end
            end
            set_parameters; % update edit boxes
        else
            P = struct2cell(pars);
            F = fieldnames(pars);
            
            % Update edit boxes with current parameters
            for p = 1 : size(parInfo,1)
                
                set(findobj(gcf,'Tag',parInfo{p,2}),...
                    'String',P{strcmp(F,parInfo{p,2})}(1))
            end
        end
    end


    function set_contrast(src,events)
        % Get current limits
        clim = get(handle.axes,'CLim');        
        oldvalue = get(src,'Value');
        newvalue    = str2double(get(src,'String'));
        switch get(src,'Tag')
            case 'MinContrast'
                
                if newvalue >= clim(2)
                    warning('Minimum contrast cannot be larger than the maximum contrast')
                    % Restore previous value
                    newvalue = oldvalue;
                end
                
                set(handle.axes,'CLim',[newvalue clim(2)])
            case 'MaxContrast'
                
                if newvalue <= clim(1)
                    warning('Maximum contrast cannot be smaller than the minimum contrast')
                    % Restore previous value
                    newvalue = oldvalue;
                end
                
                set(handle.axes,'CLim',[clim(1) newvalue])
            case 'AutoContrast'
                set(handle.axes,'CLimMode','auto')
                
        end
        clim = get(handle.axes,'CLim');
        set(handle.MinContrast,'String',sprintf('%d',clim(1)),...
            'Value',clim(1));
        set(handle.MaxContrast,'String',sprintf('%d',clim(2)),...
            'Value',clim(2));
    end
end % of main function
