function InputPointsWriter( filename,points,type )
%INPUTPOINTSWRITER Write a file with points that can be read by
%transformix.


if nargin <3
    % type is not defined. Assume that the points are provided as voxel
    % indices.
    type = 'index';
end

if ~any(strcmp({'index','point'},type))
    error('Invalid input ''type''. ''type'' must either be ''index'' or ''point''.')
end

nPoints = size(points,1);

% Write a file with the points in the transformix file format
fid = fopen(filename,'w');
fprintf(fid,'%s\n',type);
fprintf(fid,'%d\n',nPoints);
for k = 1 : nPoints
    fprintf(fid,'%f %f %f\n',points(k,:));
end
fclose(fid);
end % of the function

