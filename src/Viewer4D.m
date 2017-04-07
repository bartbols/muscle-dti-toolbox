function Viewer4D(varargin)
%% VIEWER4D_NII displays slices of a 4D nifti file. Using two slider bars 
% you can scroll through the 3rd and 4th dimension of the data.
%
% Bart Bolsterlee, NeuRA
% April 2017
%
% INPUT:
% first input : full filename of a nifti file (extension .nii or .nii.gz) OR
%               structure array with 4D data loaded from a nifti file by load_untouch_nii
%
% USAGE:
% Viewer4D_nii     (no input; dialog box opens to select a file)
% Viewer4D_nii('filename.nii.gz')
% Viewer4D_nii(nii_structure_array)

%% Check input arguments
if nargin == 0
    % No input arguments provided. Ask user to select a file.
    [fname,path] = uigetfile({'*.nii.gz;*.nii','NIFTI-files (*.nii, *.nii.gz' },'Select a 4D nifti file');
    nii = fullfile(path,fname);
elseif nargin == 1
    nii = varargin{1};
elseif nargin > 1
    error('Wrong number of input arguments (must be 0 or 1)')
end

%% Check if structure array or filename is provided as input
if ~isstruct(nii)
    % Input is not a structure, so must be a filename of a nifti file.
    % Load the file into the workspace, if the file exists.
    if exist(nii,'file') ~= 2
        error('NIFTI-file %s does not exist',nii)
    else
        filename = nii;
        fprintf('Loading %s... ',filename)
        nii = load_untouch_nii(filename);
        fprintf('completed.\n')
    end
else
    filename = '';
end

%% Get image dimensions
dim = nii.hdr.dime.dim(2:5);
origin = [0 0 0];
voxelsize = nii.hdr.dime.pixdim(2:4);

%% Create a figure and some plot objects
figure('Name',sprintf('4D image viewer: %s',filename))
hold on
colormap gray
im_handle = imagesc((0:1:dim(1)-1)*voxelsize(1) + origin(1),...
    (0:1:dim(2)-1)*voxelsize(2) + origin(2),...
    NaN(dim(2),dim(1)));
txt_handle = text(0.02,0.98,'','Units','Normalized',...
    'HorizontalAlignment','Left','VerticalAlignment','Top',...
    'FontWeight','bold','FontSize',14,'Color','w');
set(gca,'Position',[0.2 0.22 0.9 0.76],...
    'Units','Normalized',...
    'Tag','Axes',...
    'XLim',[0 (dim(1)-1)*voxelsize(1)] + origin(1),...
    'YLim',[0 (dim(2)-1)*voxelsize(2)] + origin(2))
axis equal tight
colorbar('eastoutside')

uicontrol('Style','Slider',...
    'Min',1,'Max',dim(3),...
    'SliderStep',[1/(dim(3)-1), 10/(dim(3)-1)],...
    'Units','Normalized',...
    'Position',[0.05 0.05 0.9 0.05],...
    'Value',round(dim(3)/2),...
    'CallBack',@UpdateSlice,...
    'Tag','dim3_slider');

uicontrol('Style','Slider',...
    'Min',1,'Max',dim(4),...
    'SliderStep',[1/(dim(4)-1), 1/(dim(4)-1)],...
    'Units','Normalized',...
    'Position',[0.05 0.12 0.9 0.05],...
    'Value',round(dim(4)/2),...
    'CallBack',@UpdateSlice,...
    'Tag','dim4_slider');

uicontrol('Style','pushbutton',...
    'Units','Normalized',...
    'Position',[0.05 0.23 0.15 0.05],...
    'CallBack',@ChangeCLim,...
    'String','Auto-update scaling on current slice')
% set(gca,'CLimMode','auto')
uicontrol('Style','checkbox',...
    'Units','Normalized',...
    'Position',[0.05 0.18 0.15 0.05],...
    'CallBack',@ChangeCLim,...
    'String','Hold current scaling across slices',...
    'Value',0)
% set(gca,'CLimMode','manual',...
%     'CLim',[0,0.95*max(nii.img(:))])

UpdateSlice
set(gca,'YDir','reverse')
% set(gca,'CLimMode','manual')

function UpdateSlice(varargin)
% Update the slice
k = round(get(findobj(gcf,'Tag','dim3_slider'),'Value'));
l = round(get(findobj(gcf,'Tag','dim4_slider'),'Value'));

set(im_handle,'CData',squeeze(nii.img(:,:,k,l))')
set(txt_handle,'String',sprintf('Slice %d/%d, gradient %d/%d',k,dim(3),l,dim(4)))
end

function ChangeCLim(src,event)
switch src.Style
    case 'pushbutton'
        % Push button was clicked. Auto-scale and reset CLimMode to
        % current value
        state_in = get(gca,'CLimMode');
        set(gca,'CLimMode','auto') % this auto-scales the axis
        drawnow
        set(gca,'CLimMode',state_in) % this sets the CLimMode to its incoming state
    case 'checkbox'
        % Checkbox was ticked
        if get(src,'Value') == 0
            set(gca,'CLimMode','auto')
        else
            set(gca,'CLimMode','manual')
        end
end

end

end% of function





