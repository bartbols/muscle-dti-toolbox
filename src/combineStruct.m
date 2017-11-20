function Sout = combineStruct(varargin)
%MERGESTRUCT Combine structures with dissimilar fields. Only those fields
%that all structures have in common are returned in a structure with
%multiple dimensions. Inputs can be:
% - a cell array containing a list of MAT-filenames
% - list of MAT-filenames or structure as separate inputs
%
% Bart Bolsterlee, Neuroscience Research Australia (NeuRA)
% November 2017
%
% Examples:
% Sout = combineStruct({'filename1.mat','filename2.mat'})
% Sout = combineStruct('filename.mat','filename2.mat','filename3.mat')
% Sout = combineStruct(S1,S2,S3,....,Sn)

F = cell(0);
C = cell(0);
if nargin == 1 && iscell(varargin{1})
    % Filenames have been provided. Load file and stack fields from all
    % files in one array
    fnames = varargin{1};
    N = length(fnames);
    for i = 1 : N
        S = load(fnames{i});
        F = [F;fieldnames(S)];
        C = [C;struct2cell(S)];
    end
    
else
   N = nargin;
   for i = 1 : N
       if ischar(varargin{i})
           % Filename is provided. Load the file in structure S.
           S = load(varargin{i});
       else
           S = varargin{i};
       end
       
       % Stack all fields together.
       F = [F; fieldnames(S)];
       C = [C; struct2cell(S)];
   end
end

[Fu,~,c] = unique(F);

COMB = cell(length(Fu),N);
names = cell(length(Fu),1);
k=0;
for f = 1 : length(Fu)
    if sum(c==f) == N
        % Include only if entry occurs the same number of times as the
        % number of files
        k=k+1;
        COMB(k,:) = C(c==f);
        names(k,1) = F(find(c==f,1));
    end
end
COMB = COMB(1:k,:);
names = names(1:k);
Sout = cell2struct(COMB,names);

end 