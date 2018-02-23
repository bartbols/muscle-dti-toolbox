function varargout = remove_angulation( filename_in,filename_out )
%REMOVE_ANGULATION removes the angulation from the NIfTI header and saves
%as a new NIfTI file.
%
% Bart Bolsterlee, Neuroscience Research Australia
% December 2017

% Load the NIfTI file
if ~isstruct(filename_in)
    img_in = load_untouch_nii(filename_in);
else
    img_in = filename_in;
end
clear filename_in

sgn = zeros(1,3);

% Get sign of the x, y and z vectors.
sgn(1) = sign(img_in.hdr.hist.srow_x(1));
sgn(2) = sign(img_in.hdr.hist.srow_y(2));
sgn(3) = sign(img_in.hdr.hist.srow_z(3));

% Get the voxel dimensions.
pixdim = img_in.hdr.dime.pixdim(2:4);


% Overwrite the srow info after removal of the angulation.
img_out = img_in;
img_out.hdr.hist.srow_x = [sgn(1)*pixdim(1) 0 0 img_in.hdr.hist.srow_x(4)];
img_out.hdr.hist.srow_y = [0 sgn(2)*pixdim(2) 0 img_in.hdr.hist.srow_y(4)];
img_out.hdr.hist.srow_z = [0 0 sgn(3)*pixdim(3) img_in.hdr.hist.srow_z(4)];

% Change quaternion parameters accordingly. This has been cross-validated
% with the quaternion definitions in ITK-SNAP.
if sgn(1) == -1 && sgn(2) == -1
    img_out.hdr.hist.quatern_b = 0;
    img_out.hdr.hist.quatern_c = 0;
    img_out.hdr.hist.quatern_d = 1;
elseif sgn(1) == -1 && sgn(2) == 1
    img_out.hdr.hist.quatern_b = 0;
    img_out.hdr.hist.quatern_c = 1;
    img_out.hdr.hist.quatern_d = 0;
elseif sgn(1) == 1 && sgn(2) == -1
    img_out.hdr.hist.quatern_b = 1;
    img_out.hdr.hist.quatern_c = 0;
    img_out.hdr.hist.quatern_d = 0;
elseif sgn(1) == 1 && sgn(2) == 1
    img_out.hdr.hist.quatern_b = 0;
    img_out.hdr.hist.quatern_c = 0;
    img_out.hdr.hist.quatern_d = 0;
end

save_untouch_nii(img_out,filename_out)
fprintf('Angulation removed and result saved as %s.\n',filename_out)

if nargout > 0
    varargout{1} = img_out;
end
end

