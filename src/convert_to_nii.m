function convert_to_nii(fname)
%%CONVERT_to_NII converts a .mha file to a NIfTI file with the same name.
if nargin == 0
    [filename,pathname] = uigetfile('*.mha','Select a distance map');
    fname = fullfile(pathname,filename);
end

system(sprintf('c3d %s -flip xy -o %s',fname,strrep(fname,'.mha','.nii.gz')))
fprintf('%s converted to %s\n',fname,strrep(fname,'.mha','.nii.gz'));
