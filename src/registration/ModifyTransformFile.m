function ModifyTransformFile( ParFileName, img_filename,ParFileNameCorrected )
%MODIFYTRANSFORMFILE modifies the parameters Size, Spacing and Origin
% in the elastix paramter file (first input) according to the values in the
% image file (second input) and saves the modified transformation file with
% a new name (third input)


if ~isstruct(img_filename)
    % Load the image.
    img = load_untouch_nii(img_filename);
else
    % NIfTI structure is provided.
    img = img_filename;
end

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
        % of differences in coordinate systems between NIFTI (our data) and
        % ITK (in which elastix works).
        A{i} = sprintf('(Origin %.5f %.5f %.5f)',...
            -img.hdr.hist.qoffset_x,...
            -img.hdr.hist.qoffset_y,...
             img.hdr.hist.qoffset_z);
    elseif startsWith(tline,'(Direction')
        % Change orientation of direction cosines because of differences in 
        % coordinate systems between NIFTI (our data) and
        % ITK (in which elastix works).
        [~,R] = makeT_from_quat(img);
        if img.hdr.dime.pixdim(1) == -1
            A{i} = sprintf('(Direction %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f)',...
                -R(1),-R(2),R(3),-R(4),-R(5),R(6),R(7),R(8),-R(9));
        elseif img.hdr.dime.pixdim(1) == 1
            A{i} = sprintf('(Direction %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f)',...
                -R(1),-R(2),R(3),-R(4),-R(5),R(6),-R(7),-R(8),--R(9));
        end
         
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

