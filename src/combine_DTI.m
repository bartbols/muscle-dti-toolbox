function combine_DTI( B0,DW,newname )
%COMBINE_DTI This function combines the B0 image and diffusion weighted
%images into one scan.

B0 = load_untouch_nii(B0);
DW = load_untouch_nii(DW);

if B0.hdr.dime.scl_slope ~= DW.hdr.dime.scl_slope
    error('scl_slope is not equal')
end

new = B0;
new.img = zeros([B0.hdr.dime.dim(2:4) DW.hdr.dime.dim(5)+1]);
new.hdr.dime.dim(1) = 4;
new.hdr.dime.dim(5) = DW.hdr.dime.dim(5)+1;
new.img(:,:,:,1) = B0.img;
new.img(:,:,:,2:end) = DW.img;
new.hdr.dime.glmax = max(new.img(:));
save_untouch_nii(new,newname);
fprintf('Combined file saved as %s\n',newname)
end

