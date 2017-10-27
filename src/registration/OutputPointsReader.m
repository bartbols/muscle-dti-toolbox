
function out =  OutputPointsReader(filename)
% OUTPUTPOINTSREADER  Reads in data from stored in 'filename', which should
% be a file with transformed points that is genereated with Transformix

fid = fopen(filename);
% formatstr = 'Point %d; %s = [ %d %d %d ] ; %s = [ %f %f %f ] ; %s = [ %d %d %d] ; %s = [ %f %f %f ] ; %s = [ %f %f %f ] ; %s = [ %d %d %d ]';
formatstr = 'Point %d; %s = [ %d %d %d ] ; %s = [ %f %f %f ] ; %s = [ %d %d %d] ; %s = [ %f %f %f ] ; %s = [ %f %f %f ]';
mydata = textscan(fid,formatstr);
fclose(fid);

% Index of points that were transformed
% PointNumber = mydata{1};
out.InputIndex(:,1) = mydata{3};
out.InputIndex(:,2) = mydata{4};
out.InputIndex(:,3) = mydata{5};

% Physical location of points that were transformed
out.InputPoint(:,1) = mydata{7};
out.InputPoint(:,2) = mydata{8};
out.InputPoint(:,3) = mydata{9};

% Index of output points in fixed image coordinates
out.OutputIndexFixed(:,1) = mydata{11};
out.OutputIndexFixed(:,2) = mydata{12};
out.OutputIndexFixed(:,3) = mydata{13};

% Physical location of output points
out.OutputPoint(:,1) = mydata{15};
out.OutputPoint(:,2) = mydata{16};
out.OutputPoint(:,3) = mydata{17};

% Deformation of points
out.Deformation(:,1) = mydata{19};
out.Deformation(:,2) = mydata{20};
out.Deformation(:,3) = mydata{21};

% % Index of output points in moving image coordinates
% out.OutputIndexMoving(:,1) = mydata{23};
% out.OutputIndexMoving(:,2) = mydata{24};
% out.OutputIndexMoving(:,3) = mydata{25};


