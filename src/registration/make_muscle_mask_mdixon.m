function make_muscle_mask_mdixon( mdixon,mask,threshold,stack,nErode )
%MAKE_MUSCLE_MASK_MDIXON creates a mask with only the muscle tissue. This
%can be used for targeted registration.


%% Load the data and perform the muscle segmentation
mdixon = load_untouch_nii(mdixon);
segm = (mdixon.img(:,:,:,stack)*mdixon.hdr.dime.scl_slope)>threshold;

% Erode a number of times
for i = 1 : nErode
    segm = imerode(segm,ones(3,3,3));
end
% Extract the largest connected component
CC = bwconncomp(segm,26);
[~,idx] = max(cellfun(@length,CC.PixelIdxList));

segm = zeros( mdixon.hdr.dime.dim(2:4));
segm(CC.PixelIdxList{idx}) = 1;

% Dilate again
for i = 1 : nErode
    segm = imdilate(segm,ones(3,3,3));
end

for k = 1 : size(segm,3)
    segm(:,:,k) = imfill(segm(:,:,k),'holes');
end
%% Save as mask
M = mdixon;
M.img = cast(segm,'like',mdixon.img);
M.hdr.dime.dim(1) = 3;
M.hdr.dime.dim(5) = 1;
M.hdr.dime.scl_slope = 1;
save_untouch_nii(M,mask)
fprintf('Mask saved as %s.\n',mask)

end

