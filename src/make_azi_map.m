function varargout = make_azi_map( fname_EV1,fname_azi )
%MAKE_AZI_MAP
EV1 = load_untouch_nii(fname_EV1);

azi_map = EV1;
azi_map.img = atan2d(EV1.img(:,:,:,2),EV1.img(:,:,:,1));
azi_map.hdr.dime.scl_slope = 1;
azi_map.hdr.dime.dim([1 5]) = [3 1];
save_untouch_nii(azi_map,fname_azi);
fprintf('Polar map saved as %s.\n',fname_azi)
if nargout > 0
    varargout{1} = azi_map;
end
end

