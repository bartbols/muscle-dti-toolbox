function ModifyTransformFile( ParFileName, img_filename,ParFileNameCorrected )
%MODIFYTRANSFORMFILE modifies the parameters Size, Spacing and Origin
% in the elastix paramter file (first input) according to the values in the
% image file (second input) and saves the modified transformation file with
% a new name (third input)

% Load the image
img = load_untouch_nii(img_filename);

% Read the file into cell A and make the change 
fid = fopen(ParFileName);
i = 1;
tline = fgetl(fid);
A{i} = tline;
while ischar(tline) 
    i = i+1;
    tline = fgetl(fid);
    if ~ischar(tline)
        break
    end
    if startsWith(tline,'(Size')
        A{i} = sprintf('(Size %d %d %d)',img.hdr.dime.dim(2:4));
    elseif startsWith(tline,'(Spacing')
        A{i} = sprintf('(Spacing %.5f %.5f %.5f)',img.hdr.dime.pixdim(2:4));
    elseif startsWith(tline,'(Origin')
        % minus signs before the x- and y-component are necessary because
        % of the nifti coordinate system.
        A{i} = sprintf('(Origin %.5f %.5f %.5f)',...
            -img.hdr.hist.qoffset_x,...
            -img.hdr.hist.qoffset_y,...
             img.hdr.hist.qoffset_z);
         
    else
        A{i} = tline;
    end
end
fclose(fid);

% Write cell A into a new text file
fid = fopen(ParFileNameCorrected, 'w');
for i = 1:numel(A)
    fprintf(fid,'%s\n', A{i});
end
fclose(fid);


end

