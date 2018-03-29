function varargout = FibreInspecter(varargin)
% FIBREINSPECTER MATLAB code for FibreInspecter.fig
%      FIBREINSPECTER, by itself, creates a new FIBREINSPECTER or raises the existing
%      singleton*.
%
%      H = FIBREINSPECTER returns the handle to a new FIBREINSPECTER or the handle to
%      the existing singleton*.
%
%      FIBREINSPECTER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FIBREINSPECTER.M with the given input arguments.
%
%      FIBREINSPECTER('Property','Value',...) creates a new FIBREINSPECTER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before FibreInspecter_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to FibreInspecter_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help FibreInspecter

% Last Modified by GUIDE v2.5 29-Mar-2018 11:28:55

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @FibreInspecter_OpeningFcn, ...
    'gui_OutputFcn',  @FibreInspecter_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

% --- Executes just before FibreInspecter is made visible.
function FibreInspecter_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to FibreInspecter (see VARARGIN)

% Choose default command line output for FibreInspecter
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

set(gcf,'MenuBar','figure')
axes(handles.axes3D)
axis equal vis3d
view(-40,5)
% UIWAIT makes FibreInspecter wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = FibreInspecter_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in raw.
function raw_Callback(hObject, eventdata, handles)
% hObject    handle to raw (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of raw


% --- Executes on button press in ext.
function ext_Callback(hObject, eventdata, handles)
% hObject    handle to ext (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ext


% --- Executes on button press in checkbox7.
function checkbox7_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox7


% --- Executes on button press in checkbox8.
function checkbox8_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox8



function tract_edit_Callback(hObject, eventdata, handles)
% hObject    handle to tract_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tract_edit as text
%        str2double(get(hObject,'String')) returns contents of tract_edit as a double


% --- Executes during object creation, after setting all properties.
function tract_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tract_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in select_tract.
function select_tract_Callback(hObject, eventdata, handles)
% hObject    handle to select_tract (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~isfield(handles,'previous_tract_path')
    handles.previous_tract_path = pwd;
end
[filename,pathstr]= uigetfile('*.mat','Select a DTI tract file',handles.previous_tract_path);

if filename == 0;return;end
handles.previous_tract_path = pathstr;
set(findobj(gcf,'Tag','tract_edit'),...
    'String',fullfile(pathstr,filename),...
    'TooltipString',fullfile(pathstr,filename))

% Load the tract data
tractfilename = get(findobj(gcf,'Tag','tract_edit'),'String');
if ~isempty(tractfilename)
    if exist(tractfilename,'file') == 2
        handles.D = load(tractfilename);
    else
        warning('Tract file %s not found.\n',tractfilename)
    end
end
if isfield(handles.D,'fibrelength')    
    % Switch off poly and extensions
    set(handles.poly,'Value',1,'Enable','on')
    set(handles.ext,'Value',0,'Enable','on')
    set(handles.raw,'Value',0,'Enable','on')
    set(handles.from_apo_to_mus,'Enable','off','Value',0)
    if isfield(handles.D,'attach_type')
        if any(handles.D.attach_type(:) == 2)            
            set(handles.from_apo_to_mus,'Enable','on','Value',1)
        end
    end
else    
    % Only raw tracts are available. Switch off poly and extension options.
    set(handles.poly,'Value',0,'Enable','off')
    set(handles.ext,'Value',0,'Enable','off')
    set(handles.raw,'Value',1,'Enable','on')
end
guidata(hObject,handles);
update_tracts_Callback(hObject, eventdata, handles)

% --- Executes on button press in select_muscle.
function select_muscle_Callback(hObject, eventdata, handles)
% hObject    handle to select_muscle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~isfield(handles,'previous_surf_path')
    handles.previous_surf_path = pwd;
end

[filename,pathstr]= uigetfile('*.stl','Select the muscle model',handles.previous_surf_path);
if filename == 0;return;end
handles.previous_surf_path = pathstr;
set(findobj(gcf,'Tag','mus_edit'),...
    'String',fullfile(pathstr,filename),...
    'TooltipString',fullfile(pathstr,filename))

% Load the muscle model
musfilename = get(findobj(gcf,'Tag','mus_edit'),'String');
muscle = [];
if ~isempty(musfilename)
    if exist(musfilename,'file') == 2
        muscle = stlread(musfilename);
    else
        warning('Muscle model %s not found.\n',musfilename)
    end
end
% Plot the muscle
if ~isempty(muscle)
    if ~isfield(handles,'muscle')
        axes(handles.axes3D)
        % Create a new muscle handle.
        handles.muscle = patch('Vertices',muscle.vertices,'Faces',muscle.faces,...
            'FaceColor',get(handles.set_color_mus,'BackgroundColor'),...
            'FaceAlpha',str2double(get(handles.set_alpha_mus,'String')),...
            'EdgeColor','k',...
            'EdgeAlpha',0.1);
        guidata(hObject,handles)
    else
        % Handle already exists. Just update.
        set(handles.muscle,'Vertices',muscle.vertices,'Faces',muscle.faces,...
            'FaceColor',get(handles.set_color_mus,'BackgroundColor'),...
            'FaceAlpha',str2double(get(handles.set_alpha_mus,'String')),...
            'EdgeColor','k',...
            'EdgeAlpha',0.1);
    end
end
guidata(hObject,handles);

% --- Executes on button press in select_apo.
function select_apo_Callback(hObject, eventdata, handles)
% hObject    handle to select_apo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~isfield(handles,'previous_surf_path')
    handles.previous_surf_path = pwd;
end

[filename,pathstr]= uigetfile('*.stl','Select the aponeurosis model',handles.previous_surf_path);
if filename == 0;return;end
handles.previous_surf_path = pathstr;
set(findobj(gcf,'Tag','apo_edit'),...
    'String',fullfile(pathstr,filename),...
    'TooltipString',fullfile(pathstr,filename))

% Load the aponeurosis model
apo = [];
apofilename = get(findobj(gcf,'Tag','apo_edit'),'String');
if ~isempty(apofilename)
    if exist(apofilename,'file') == 2
        apo = stlread(apofilename);
    else
        warning('Aponeurosis model %s not found.\n',apofilename)
    end
end
% Plot the aponeurosis
if ~isempty(apo)
    if ~isfield(handles,'apo')
        % Create a new aponeurosis handle.
        axes(handles.axes3D)
        handles.apo = patch('Vertices',apo.vertices,'Faces',apo.faces,...
            'FaceColor',get(handles.set_color_apo,'BackgroundColor'),...
            'FaceAlpha',str2double(get(handles.set_alpha_apo,'String')),...
            'EdgeColor','k',...
            'EdgeAlpha',0.1);
        guidata(hObject,handles)
    else
        % Handle already exists. Just update.
        set(handles.apo,'Vertices',apo.vertices,'Faces',apo.faces,...
            'FaceColor',get(handles.set_color_apo,'BackgroundColor'),...
            'FaceAlpha',str2double(get(handles.set_alpha_apo,'String')),...
            'EdgeColor','k',...
            'EdgeAlpha',0.1);
    end
end

function mus_edit_Callback(hObject, eventdata, handles)
% hObject    handle to mus_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mus_edit as text
%        str2double(get(hObject,'String')) returns contents of mus_edit as a double


% --- Executes during object creation, after setting all properties.
function mus_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mus_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function apo_edit_Callback(hObject, eventdata, handles)
% hObject    handle to apo_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of apo_edit as text
%        str2double(get(hObject,'String')) returns contents of apo_edit as a double


% --- Executes during object creation, after setting all properties.
function apo_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to apo_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in update_tracts.
function update_tracts_Callback(hObject, eventdata, handles)
% hObject    handle to update_tracts (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

axes(handles.axes3D)
% Update the tracts
if ~isempty(handles.D)
    hold on
    % Get selection from table
    if strcmp(hObject.Tag,'select_tract')
        % New tracts were loaded. Display all values and automatically set
        % the thresholds.
        if ~isfield(handles.D,'ang') && isfield(handles.D,'endpoints_dir')
            handles.D.ang = acosd(sum(squeeze(handles.D.endpoints_dir(:,1,:)) .* squeeze(handles.D.endpoints_dir(:,2,:)),2));
        end
        if isfield(handles.D,'ext') && ~isfield(handles.D,'abs_ext')
            handles.D.abs_ext = nansum(handles.D.ext,2);
        end
        if isfield(handles.D,'penangle') && ~isfield(handles.D,'pennation')
            handles.D.pennation = nanmean(handles.D.penangle,2);
        end
        if isfield(handles.D,'fibrelength')
            handles.threshold_table.Data = {...
                floor(min(handles.D.pct_ext))    , ceil(max(handles.D.pct_ext));...
                floor(min(handles.D.abs_ext))    , ceil(max(handles.D.abs_ext));...
                floor(min(handles.D.fibrelength)), ceil(max(handles.D.fibrelength));...
                floor(min(handles.D.pennation))  , ceil(max(handles.D.pennation));...
                floor(min(handles.D.curvature))  , ceil(max(handles.D.curvature));...
                floor(min(handles.D.ang))        , ceil(max(handles.D.ang));...
                floor(min(handles.D.length_mm))  , ceil(max(handles.D.length_mm))};
        else
            % Only raw tracts are available
            handles.threshold_table.Data = {NaN,NaN;NaN NaN;NaN NaN;...
                NaN,NaN;NaN NaN;NaN NaN;...
                floor(min(handles.D.length_mm)),ceil(max(handles.D.length_mm))};
        end
    end
    
    % Get thresholds from the table
    thresholds = handles.threshold_table.Data;
    
    if isfield(handles.D,'fibrelength')
        if get(handles.from_apo_to_mus,'Value')
            is_from_apo_to_mus = (handles.D.attach_type(:,1) ~= handles.D.attach_type(:,2));
        else
            is_from_apo_to_mus = true(size(handles.D.pct_ext,1),1);
        end
        selection = find(handles.D.pct_ext     > thresholds{1,1} & handles.D.pct_ext     < thresholds{1,2} & ...
            handles.D.abs_ext     > thresholds{2,1} & handles.D.abs_ext     < thresholds{2,2} & ...
            handles.D.fibrelength > thresholds{3,1} & handles.D.fibrelength < thresholds{3,2} & ...
            handles.D.pennation   > thresholds{4,1} & handles.D.pennation   < thresholds{4,2} & ...
            handles.D.curvature   > thresholds{5,1} & handles.D.curvature   < thresholds{5,2} & ...
            handles.D.ang         > thresholds{6,1} & handles.D.ang         < thresholds{6,2}  & ...
            handles.D.length_mm   > thresholds{7,1} & handles.D.length_mm   < thresholds{7,2} & ...
            is_from_apo_to_mus)';
    else
        selection = find(handles.D.length_mm   > thresholds{7,1} & handles.D.length_mm   < thresholds{7,2});
    end
    nSel = length(selection);
    nFib = size(handles.D.fibindex,1);
    set(handles.Nfibres,'String',sprintf('Displaying %d of %d fibres (%.1f%%)',...
        nSel,nFib,nSel/nFib*100))
    
    if get(handles.poly,'Value')
        % Plot the polynomial fitted tracts, including extrapolations
        PlotX = [];PlotY = [];PlotZ = [];
        if isfield(handles.D,'PolyCoeff')
            P = handles.D.PolyCoeff;
        else
            error('''PolyCoeff'' not found as a field in handles.D. Polynomials fits cannot be plotted.')
        end
        for fibnr = selection
            if isempty(P(fibnr).x);continue;end
            if any(isnan(handles.D.fibindex_trunc(fibnr,1:2)))
                continue
            end
            
            tmpX = [handles.D.endpoints(fibnr,1,1) polyval(P(fibnr).x,linspace(P(fibnr).t0,P(fibnr).t1,100)) handles.D.endpoints(fibnr,2,1)];
            tmpY = [handles.D.endpoints(fibnr,1,2) polyval(P(fibnr).y,linspace(P(fibnr).t0,P(fibnr).t1,100)) handles.D.endpoints(fibnr,2,2)];
            tmpZ = [handles.D.endpoints(fibnr,1,3) polyval(P(fibnr).z,linspace(P(fibnr).t0,P(fibnr).t1,100)) handles.D.endpoints(fibnr,2,3)];
            PlotX = [PlotX NaN tmpX];
            PlotY = [PlotY NaN tmpY];
            PlotZ = [PlotZ NaN tmpZ];
        end
        if ~isfield(handles,'tracts_poly')
            % Create new handle.
            handles.tracts_poly = plot3(PlotX(:),PlotY(:),PlotZ(:),...
                'LineWidth',1,...
                'Color',get(handles.set_color_poly,'BackgroundColor'));
        else

            % Update polynomials.
            set(handles.tracts_poly,'XData',PlotX(:),'YData',PlotY(:),'ZData',PlotZ(:))
        end
    else
        if isfield(handles,'tracts_poly')
            set(handles.tracts_poly,'XData',NaN,'YData',NaN,'ZData',NaN)
        end
        
    end
    if get(handles.raw,'Value')
        % Plot the raw tracts
        PlotX = [];PlotY =[];PlotZ = [];
        for fibnr = selection
            PlotX = [PlotX handles.D.tracts_xyz(1,min(handles.D.fibindex(fibnr,1:2)):max(handles.D.fibindex(fibnr,1:2))) NaN];
            PlotY = [PlotY handles.D.tracts_xyz(2,min(handles.D.fibindex(fibnr,1:2)):max(handles.D.fibindex(fibnr,1:2))) NaN];
            PlotZ = [PlotZ handles.D.tracts_xyz(3,min(handles.D.fibindex(fibnr,1:2)):max(handles.D.fibindex(fibnr,1:2))) NaN];
        end
        if ~isfield(handles,'tracts_raw')
            % Create new handle.
            handles.tracts_raw = plot3(PlotX(:),PlotY(:),PlotZ(:),...
                'LineWidth',1,...
                'Color',get(handles.set_color_raw,'BackgroundColor'));
        else
            % Update raw tracts.
            set(handles.tracts_raw,'XData',PlotX(:),'YData',PlotY(:),'ZData',PlotZ(:))
        end
    else
        if isfield(handles,'tracts_raw')
             set(handles.tracts_raw,'XData',NaN,'YData',NaN,'ZData',NaN)
        end
    end
    if get(handles.ext,'Value')
        % Plot the extended parts
        PlotX = [];PlotY = [];PlotZ = [];
        for fibnr = selection
            if isfield(handles.D,'PolyCoeff')
                P = handles.D.PolyCoeff;
            else
                error('''PolyCoeff'' not found as a field in handles.D. Polynomials fits cannot be plotted.')
            end
            
            tmpX = [handles.D.endpoints(fibnr,1,1) polyval(P(fibnr).x,P(fibnr).t0) NaN polyval(P(fibnr).x,P(fibnr).t1)  handles.D.endpoints(fibnr,2,1)];
            tmpY = [handles.D.endpoints(fibnr,1,2) polyval(P(fibnr).y,P(fibnr).t0) NaN polyval(P(fibnr).y,P(fibnr).t1)  handles.D.endpoints(fibnr,2,2)];
            tmpZ = [handles.D.endpoints(fibnr,1,3) polyval(P(fibnr).z,P(fibnr).t0) NaN polyval(P(fibnr).z,P(fibnr).t1)  handles.D.endpoints(fibnr,2,3)];
            PlotX = [PlotX NaN tmpX];
            PlotY = [PlotY NaN tmpY];
            PlotZ = [PlotZ NaN tmpZ];
            if ~isfield(handles,'tracts_ext')
                % Create new handle.
                handles.tracts_ext = plot3(PlotX(:),PlotY(:),PlotZ(:),...
                    'LineWidth',1,...
                    'Color',get(handles.set_color_ext,'BackgroundColor'));
            else
                % Update extensions.
                set(handles.tracts_ext,'XData',PlotX(:),'YData',PlotY(:),'ZData',PlotZ(:))
            end
        end
    else
        if isfield(handles,'tracts_ext')
            set(handles.tracts_ext,'XData',NaN,'YData',NaN,'ZData',NaN)
        end
    end
    
    % Update histograms
    VARS = {'Raw tract length','mm',2,'length_mm';...
        '% extension','%',2,'pct_ext';...
        'Abs. extension','mm',2,'abs_ext';...
        'Fascicle length','mm',2,'fibrelength';...
        'Pennation angle','degr',2,'pennation';...
        'Curvature','1/m',2,'curvature';...
        'Angle between endpoint slopes','degr',2,'ang';...
        'Fractional anisotropy','-',0.01,'fa';...
        'Mean diffusivity','1e-9 mm^2/s',0.05,'md';...
        '\lambda_1','1e-3 mm^2/s',0.05,'lambda1';...
        '\lambda_2','1e-3 mm^2/s',0.05,'lambda2';...
        '\lambda_3','1e-3 mm^2/s',0.05,'lambda3'};
    
    CELLDATA = struct2cell(handles.D);
    varNames = fieldnames(handles.D);
    for k = 1 : size(VARS,1)
        % Decide which subplot to display the data in
        eval(sprintf('axes(handles.axes%d);',k))
        set(gca,'Visible','on')
        % Check if variable exist. If not, continue with next.
        idx = strcmp(varNames,VARS{k,4});
        if ~any(idx)
            set(gca,'Visible','off')
            cla
            title([VARS{k,1} ' data not available'])
        else
            set(gca,'Visible','on')
            % Get the data
            data = mean(CELLDATA{idx},2);

            % Filter the data to only plot histograms for the selected fibres
            % (if not selection is provided, histograms for all fibres are
            % provided)
            data = data(selection);

            % Plot a histogram
            MyHist( data,VARS{k,3},VARS{k,1},VARS{k,2},...
                'TextPos',[0.95,0.65 ],...
                'Color','b','BarAlpha',0.5,...
                'FontWeight','bold',...
                'EdgeColor','none',...
                'TextColor','k');
            set(gca,'XLim',[min(data)-VARS{k,3} max(data+VARS{k,3})])
        end
        set(get(gca,'Title'),'Units','Normalized',...
            'Position',[0.95,0.95],...
            'HorizontalAlignment','right',...
            'VerticalAlignment','top',...
            'String',{VARS{k,1},sprintf('bin = %.2f',VARS{k,3})})
        ticklabelinside(gca,'y')
        if ~any(k==[1 5 9])
            set(get(gca,'YLabel'),'String','')
%             % Bring ylabel inwards
%             set(get(gca,'YLabel'),'Units','Normalized')
%             pos = get(get(gca,'YLabel'),'Position');
%             set(get(gca,'YLabel'),'Position',pos + [0.01 0 0])
        end
    end
    
end

guidata(hObject, handles);

% --- Executes on button press in set_color.
function set_color_Callback(hObject, eventdata, handles)
% hObject    handle to set_color_mus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
color = uisetcolor(get(hObject,'BackgroundColor'),'Select a color for the aponeurosis');
set(hObject,'BackgroundColor',color);
switch get(hObject,'tag')
    case 'set_color_mus'
        if isfield(handles,'muscle')
            set(handles.muscle,'FaceColor',color)
        end
    case 'set_color_apo'
        if isfield(handles,'apo')
            set(handles.apo,'FaceColor',color)
        end
    case 'set_color_poly'
        if isfield(handles,'tracts_poly')
            set(handles.tracts_poly,'color',color)
        end
    case 'set_color_raw'
        if isfield(handles,'tracts_raw')
            set(handles.tracts_raw,'color',color)
        end
    case 'set_color_ext'
        if isfield(handles,'tracts_ext')
            set(handles.tracts_ext,'color',color)
        end
end


% --- Executes on button press in poly.
function poly_Callback(hObject, eventdata, handles)
% hObject    handle to poly (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of poly
%

function set_alpha_Callback(hObject, eventdata, handles)
% hObject    handle to set_alpha (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mus_alpha_txt as text
%        str2double(get(hObject,'String')) returns contents of set_alpha as a double

alpha_value = str2double(get(hObject,'String'));
if isnan(alpha_value) || alpha_value < 0 || alpha_value > 1
    alpha_value = 0.2;
    set(hObject,'String',sprintf('%0.1f',alpha_value))
end

switch get(hObject,'tag')
    case 'set_alpha_mus'
        if isfield(handles,'muscle')
            set(handles.muscle,'FaceAlpha',alpha_value)
        end
    case 'set_alpha_apo'
        if isfield(handles,'apo')
            set(handles.apo,'FaceAlpha',alpha_value)
        end
end

% --- Executes during object creation, after setting all properties.
function mus_alpha_txt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mus_alpha_txt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function set_alpha_mus_CreateFcn(hObject, eventdata, handles)
% hObject    handle to set_alpha_apo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of set_alpha_apo as text
%        str2double(get(hObject,'String')) returns contents of set_alpha_apo as a double


% --- Executes during object creation, after setting all properties.
function set_alpha_apo_CreateFcn(hObject, eventdata, handles)
% hObject    handle to set_alpha_apo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in from_apo_to_mus.
function from_apo_to_mus_Callback(hObject, eventdata, handles)
% hObject    handle to from_apo_to_mus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of from_apo_to_mus


% --- Executes on button press in delete_tract.
function delete_tract_Callback(hObject, eventdata, handles)
% hObject    handle to delete_tract (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isfield(handles,'D')
    handles = rmfield(handles,'D');
end
set(handles.tract_edit,'String',[]);
% Delete tracts graphics objects
if isfield(handles,'tracts_poly')
    delete(handles.tracts_poly);
    handles = rmfield(handles,'tracts_poly');
end
if isfield(handles,'tracts_raw')
    delete(handles.tracts_raw);
    handles = rmfield(handles,'tracts_raw');
end
if isfield(handles,'tracts_ext')
    delete(handles.tracts_ext);
    handles = rmfield(handles,'tracts_ext');
end
for k = 1 : 12
    eval(sprintf('axes(handles.axes%d);',k))
    cla
    set(gca,'Visible','off')
end
set(handles.Nfibres,'String','')

guidata(hObject, handles);
    

% --- Executes on button press in delete_mus.
function delete_mus_Callback(hObject, eventdata, handles)
% hObject    handle to delete_mus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isfield(handles,'muscle')
    delete(handles.muscle)
    handles = rmfield(handles,'muscle');    
end
set(handles.mus_edit,'String',[]);
guidata(hObject, handles);


% --- Executes on button press in delete_apo.
function delete_apo_Callback(hObject, eventdata, handles)
% hObject    handle to delete_apo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isfield(handles,'apo')
    delete(handles.apo)
    handles = rmfield(handles,'apo');    
end
set(handles.apo_edit,'String',[]);
guidata(hObject, handles);
