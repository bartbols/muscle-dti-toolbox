clear
s = 1;
filename = fullfile('../data',sprintf('%03d',s),...
    sprintf('STDTI%03d_mdixon_left.nii.gz',s));

% filename = 'C:\Users\b.bolsterlee\Documents\Working\ShortLongStudy\data\MRI\MRI_s01_j01.nii.gz';


close all
muscle_segmenter(filename)
