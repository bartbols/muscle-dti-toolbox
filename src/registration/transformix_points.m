function [ points_out ] = transformix_points( points,transform_file )
%TRANSFORMIX_POINTS transforms 'points' (n x 3 array with physical
% coordinates) and transforms these points with the Elastix transformation
% file 'transformfile' using Transformix. The transformed points are
% returned.
%
% Bart Bolsterlee, Neuroscience Research Australia
% November 2017
%
% Example:
% points_out = transformix_points(points_in,transform_file)

% Create temporary working directory.
char_list = char(['a':'z' '0':'9']) ;
tmpdir = [];
while exist(tmpdir,'dir') == 7 || isempty(tmpdir)
    tmpdir = fullfile(tempdir,char_list(ceil(length(char_list)*rand(1,8))));
end
mkdir(tmpdir)

try    
        
    % Check if n x 3 or 3 x n array is provided. An array of the same
    % dimensions will be returned.
    
    n = size(points,2);
    if n ~= 3
        points = points';
    end
    % Write vertices to transformix input points file. Because of 
    % differences in ITK and NIFTI coordinate system, x- and y-values 
    % need to be flipped before transformation.
    R =  [-1 0 0;0 -1 0;0 0 1];
    
    % Check for NaNs
    nanidx = any(isnan(points),2);   
    
    InputPointsWriter(fullfile(tmpdir,'inputpoints.txt'),...
        points(~nanidx,:)*R,'point')
    
    % Build up the transformix command.
    transformix_cmd = sprintf('transformix -def %s -out %s -tp %s',...
        fullfile(tmpdir,'inputpoints.txt'),...
        tmpdir,...
        transform_file);
    
    % Run transformix
    system(transformix_cmd);
    
    % Read the transformed vertices
    out = OutputPointsReader(fullfile(tmpdir,'outputpoints.txt'));
    
    % Flip vertices back to NIFTI coordinates and save as STL file.
    points_out = NaN(size(points));    
    points_out(~nanidx,:) = out.OutputPoint*R;
    if n ~= 3
        points_out = points_out';
    end
    rmdir(tmpdir,'s')
    fprintf('\n')
    
catch ME
    % remove temporary working directory, then throw error message
    rmdir(tmpdir,'s')
    error(ME.message)
end

end

