function SetBsplineOrderToZero( ParFileName,ParFileNameCorrected )
%CORRECTBSPLINEORDER The parameter "FinalBSplineInterpolationOrder" in the
% resulting parameter file should be set to 0 to correctly transform a
% mask. This is explained in the Elastix-manual. When this is not done, the
% mask is not correctly transformed. This script opens the original
% parameter file and writes a corrected parameter file with identical
% parameters, but with FinalBSplineInterpolationOrder 0 (regardless of what
% the original value was).

% Read the file into cell A and make the change 
fid = fopen(ParFileName);
i = 1;
tline = fgetl(fid);
A{i} = tline;
while ischar(tline)
    i = i+1;
    tline = fgetl(fid);
    if ~isempty(findstr(tline,'FinalBSplineInterpolationOrder'))
        A{i} = '(FinalBSplineInterpolationOrder 0)';
    else
        A{i} = tline;
    end
end
fclose(fid);

% Write cell A into a new text file
fid = fopen(ParFileNameCorrected, 'w');
for i = 1:numel(A)
    if A{i+1} == -1
        fprintf(fid,'%s', A{i});
        break
    else
        fprintf(fid,'%s\n', A{i});
    end
end
fclose(fid);


end

