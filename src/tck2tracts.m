function DTItracts = tck2tracts(filename)
%TCK2TRACTS Converts mrtrix-generated tracts in tck format to the 
% DTItracts format used for processing in Matlab.
%

if ~endsWith(filename,'.tck')
    error('Input to tck2tracts should be a .tck file')
end
tck = read_mrtrix_tracks(filename);

% Make variables tracts_xyz and fibindex.
DTItracts.tracts_xyz = cell2mat(tck.data')';

len = cellfun(@length,tck.data);
DTItracts.fibindex(:,2)     = cumsum(len)';
DTItracts.fibindex(1,1)     = 1;
DTItracts.fibindex(2:end,1) = DTItracts.fibindex(1:end-1,2)+1;

end % of function
