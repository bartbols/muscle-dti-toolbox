function [ LabelNumbers,LabelNames,FullLabelNames ] = read_lab( filename )
%READ_LAB reads the label file and returns the label numbers, short label names
% and full label names.
%
% Bart Bolsterlee
% Neuroscience Research Australia
% October 2018
%
% USAGE
% [LabelNumbers,LabelNames,FullLabelNames] = read_lab(filename)


fid = fopen(filename,'r');
LabelInfo = textscan(fid,'%d %s %s','Delimiter',',');
fclose(fid);
LabelNumbers   = LabelInfo{1};
LabelNames     = LabelInfo{2};
FullLabelNames = LabelInfo{3};


end

