function [ output_args ] = reg_elastix( fixed, moving, parfile,transformfile,varargin )
%REG_ELASTIX Calls elastix to perform a registration on the fixed and moving
% image using parameter file 'parfile'. The resulting transformation file
% will be saved with filename 'transformfile'. 'parfile' can also be a cell
% with multiple parameter files, in which case a multi-step registration
% will be performed.
%
% ----------------- INPUT -----------------
% ----- REQUIRED -----
% fixed:    filename of the fixed image
% moving:   filename of the moving image
% parfile:  filename (char) or a cell string of filenames of the elastix
%           parameter file(s) with registration parameters. For example, to
%           perform a rigid and then a bspline registration: parfile =
%           {'rigid.txt','bspline.txt'}. Or to only perform a bspline, use
%           'bspline.txt'.
% transformfile:  filename of the final transformation file (will be
%                 created by elastix)
%
% ----- OPTIONAL -----
% Optional inputs are provided as pairs of 'ParameterName',<value> (e.g.
% 'mask','mask.nii.gz')
%
% - mask                 : filename of the mask file used for registration
% - foreground_threshold : threshold intensity for foreground. A foreground
%                          mask will be created from the fixed image using
%                          this threshold.
% - stack_f              : if fixed image is 4D, choose which stack is
%                          used for registration.
% - stack_m              : if moving is 4D, choose which stack is
%                          used for registration.
% - surface_in           : STL filename of surface model to be transformed
%                          with the resulting transformation
% - surface_out           : STL filename of surface model after
%                           transformation (will be created)

% Read the inputs
p = inputParser;
addRequired(p,'fixed',@(x) contains(x,'.nii.gz'))
addRequired(p,'moving',@(x) contains(x,'.nii.gz'))
addRequired(p,'parfile',@(x) ischar(x) || iscell(x))
addRequired(p,'transformfile',@(x) ischar(x))
addParameter(p,'mask',[],@(x) isempty(x) || contains(x,'.nii.gz'))
addParameter(p,'foreground_threshold',[],@(x) isscalar(x) || isempty(x))
addParameter(p,'stack_f',[],@(x) assert(isnumeric(x)))
addParameter(p,'stack_m',[],@(x) assert(isnumeric(x)))
addParameter(p,'surface_in',[],@(x) contains(x,'.stl'))
addParameter(p,'surface_out',[],@(x) contains(x,'.stl'))
parse(p,fixed, moving, parfile,transformfile,varargin{:});

mask        = p.Results.mask;
foreground_threshold = p.Results.foreground_threshold;
stack_f     = p.Results.stack_f;
stack_m     = p.Results.stack_m;
surface_out = p.Results.surface_out;
surface_in  = p.Results.surface_in;

% Create temporary working directory.
char_list = char(['a':'z' '0':'9']) ;
tmpdir = [];
while exist(tmpdir,'dir') == 7 || isempty(tmpdir)
    tmpdir = fullfile(pwd,char_list(ceil(length(char_list)*rand(1,8))));
end
mkdir(tmpdir)

try

    if ~isempty(stack_f)
        % Get the selected stack(s) from the 4D fixed image
        for i = 1 : length(stack_f)
            fprintf('In reg_elastix: Extracting stack %d from 4D fixed image... ',stack_f(i))
            fixed_3D{i} = fullfile(tmpdir,sprintf('fixed_image%02d.nii.gz',i));
            extract_3Dfrom4D(fixed,fixed_3D{i},stack_f(i));
            fprintf('completed.\n')
        end
    else
        fixed_3D{1} = fixed;
    end
    if ~isempty(stack_m)
        % Get the selected stack(s) from the 4D moving image
        for i = 1:length(stack_m)
            fprintf('In reg_elastix: Extracting stack %d from 4D moving image... ',stack_m(i))
            moving_3D{i} = fullfile(tmpdir,sprintf('moving_image%02d.nii.gz',i));
            extract_3Dfrom4D(moving,moving_3D{i},stack_m(i));
            fprintf('completed.\n')
        end
    else
        moving_3D{1} = moving;
    end
    
    
    if isempty(mask) && ~isempty(foreground_threshold)
        % Create mask based on foreground_threshold. Use the first fixed
        % image to make the mask.
        fprintf('In reg_elastix: Creating foreground mask with threshold %.1f...',foreground_threshold)
        mask_img = load_untouch_nii(fixed_3D{1});
        mask_img.img = cast(mask_img.img * mask_img.hdr.dime.scl_slope>foreground_threshold,'like',mask_img.img);
        mask_img.hdr.dime.scl_slope = 1;
        mask = fullfile(tmpdir,'mask.nii.gz');
        save_untouch_nii(mask_img,mask);
        fprintf('completed.\n')
        
    end
        
    % Check how many steps there are in the registration
    if ischar(parfile)
        % Only one parameter file is provided.
        parfile = {parfile};
    end
    
    nSteps = length(parfile);
    % Now, iterate through the steps, using the transformation from the
    % previous step as the initial transformation.
    for stepnr = 1 : nSteps
        mkdir(fullfile(tmpdir,sprintf('step%02d',stepnr)))
        
        % Build up the elastix command
        elastix_cmd = 'elastix ';
        if isempty(stack_f) || numel(stack_f) == 1
            % One fixed image is provided.
            elastix_cmd = [elastix_cmd ' -f ' fixed_3D{1}];
        else
            % Multiple fixed images are provided.
            for i = 1: length(stack_f)
                elastix_cmd = [elastix_cmd ' -f' int2str(i-1) ' ' fixed_3D{i}];
            end
        end
        if isempty(stack_m) || numel(stack_m) == 1
            % One fixed image is provided.
            elastix_cmd = [elastix_cmd ' -m ' moving_3D{1}];
        else
            % Multiple fixed images are provided.
            for i = 1: length(stack_m)
                elastix_cmd = [elastix_cmd ' -m' int2str(i-1) ' ' moving_3D{i}];
            end
        end
        
        % Add output directory and parameter file
        elastix_cmd = [elastix_cmd ' -out ' fullfile(tmpdir,sprintf('step%02d',stepnr)),...
                                   ' -p ' parfile{stepnr}];
        
        if ~isempty(mask)
            % Add fixed mask, if provided.
            elastix_cmd = [elastix_cmd sprintf(' -fMask %s',mask)];
            
        end
        if stepnr ~= 1
            % After first step, use previous parameter file as initial
            % transform
            elastix_cmd = [elastix_cmd sprintf(' -t0 %s',...
                fullfile(tmpdir,sprintf('step%02d',stepnr-1),'TransformParameters.0.txt'))];
        end
        % Run registration with Elastix
        system(elastix_cmd)
        pause(1)
    end
    
    if ~isempty(surface_in)
        FV = stlread(surface_in);        
        
        % Write vertices to transformix input points file. Because of 
        % differences in ITK and NIFTI coordinate system, x- and y-values 
        % need to be flipped before transformation.
        R =  [-1 0 0;0 -1 0;0 0 1];
        InputPointsWriter(fullfile(tmpdir,'inputpoints.txt'),...
            FV.vertices*R,'point')

        % Build up the transformix command.
        transformix_cmd = sprintf('transformix -def %s -out %s -tp %s',...
            fullfile(tmpdir,'inputpoints.txt'),...
            tmpdir,...
            fullfile(tmpdir,sprintf('step%02d',nSteps),...
            'TransformParameters.0.txt'));
        
        % Run transformix
        system(transformix_cmd);

        % Read the transformed vertices
        out = OutputPointsReader(fullfile(tmpdir,'outputpoints.txt'));
        
        % Flip vertices back to NIFTI coordinates and save as STL file.
        FV.vertices = out.OutputPoint*R;
        stlwrite(surface_out,FV)
        fprintf('In reg_elastix: Transformed surface saved as %s\n',surface_out)

    end
    % Copy the final transform file.
    movefile(fullfile(tmpdir,sprintf('step%02d',nSteps),...
        'TransformParameters.0.txt'),transformfile);
    
    % Delete the temporary working directory.
    rmdir(tmpdir,'s')
    
catch ME
    % remove temporary working directory, then throw error message
    rmdir(tmpdir,'s')
    error(ME.message)
end


end

