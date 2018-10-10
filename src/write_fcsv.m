function write_fcsv( points,filename )
%WRITE_FCSV writes an FCSV file with the coordinates in the n x 3 array
%(first input) to an FCSV file named <filename> (second input argument).

fid = fopen(filename,'w');

fprintf(fid,'# Markups fiducial file version = 4.7\n');
fprintf(fid,'# CoordinateSystem = 0\n');
fprintf(fid,'# columns = id,x,y,z,ow,ox,oy,oz,vis,sel,lock,label,desc,associatedNodeID\n');

% vtkMRMLMarkupsFiducialNode_0,23.3639,43.2236,77.2876,0,0,0,1,1,1,0,name1,description1,vtkMRMLModelNode8
formatstr = 'vtkMRMLMarkupsFiducialNode_%d,%.4f,%.4f,%.4f,0,0,0,1,1,1,0,,,vtkMRMLModelNode8\n';
for i = 1 : size(points,1)
    fprintf(fid,formatstr,i,points(i,1),points(i,2),points(i,3));
end
fclose(fid);
fprintf('%d points saved to %s\n',size(points,1),filename)


end

