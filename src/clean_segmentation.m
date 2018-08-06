function clean_segmentation(filename_in,filename_clean)
%% CLEAN_SEGMENTATION cleans a mask file by filling holes and selecting the
% largest connected component in each label.

if nargin == 0
    [fname,pathname] = uigetfile( {'*.nii;*.nii.gz',...
     'NIfTI Files (*.nii,*.nii.gz)'}, ...
       'Select a mask file');
    if fname == 0
        return
    end
    filename_in = fullfile(pathname,fname);
    [fname,pathname] = uiputfile({'*.nii;*.nii.gz'},'Save as',pathname);
    
    if fname == 0
        return
    end
    filename_clean = fullfile(pathname,fname);
end

fprintf('Reading %s... ',filename_in);
M = load_untouch_nii(filename_in);
fprintf('completed.\n')
label_nrs = unique(M.img(:));
label_nrs(label_nrs==0) = [];

if length(label_nrs) > 100
    error('A maximum of 100 label numbers per file can be processed. Make sure you selected a segmentation file and not an image file.')
end
new = M;
new.img = cast(zeros(size(M.img)),'like',M.img);
fprintf('The following labels were found: '),fprintf('%d ',label_nrs);fprintf('\n')
for label_nr = label_nrs'
    fprintf('Cleaning label %d...',label_nr)
    % Select the largest connected component and fill holes
    BW = bwconncomp(imfill((M.img == label_nr),'holes'));
    [~,idx] = max(cellfun(@numel,BW.PixelIdxList));
%     BW = imfill((M.img == label_nr),'holes');
    new.img(BW.PixelIdxList{idx}) = cast(label_nr,'like',M.img);
    fprintf(' completed.\n')
end
save_untouch_nii(new,filename_clean);
fprintf('Cleaned mask saved as %s\n',filename_clean)