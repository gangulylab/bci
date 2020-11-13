function [file_names f_names] = findfiles(entry, directory, isRecursive)
% function [file_names f_names] = findfiles(entry, directory, isRecursive)
% findfiles - finds files with the specified entry in specified director and
% subdirectories (if recursive).
% Arguments: entry - the string which must be included in the file names
%            directory - the directory specified by user
%            isRecursive - search subdirectories or not; 1 - yes; 0 - no;
% Returns:  a cell array with the fully specified file names
% Examples:
%    1. file_names = findfiles(entry, directory)
%       searches in the specified directory, not in subdirectories.
%    2. file_names = findfiles(entry, directory, 1)
%       searches in the specified directory and all its subdirectories.
% Notes:
%    1. directory is complusory, use pwd for the current directory;
%    2. special entries:
%       if entry = '' then return all files
%       use entry = '.???' for extension matching, for example, '.txt'
%       use entry = '~???' for names not including '???'
%       entry is case sensitive
%    3. level limitation applies when recursive

% Copyright (C) 2005 Bioeng Institute, University of Auckland
% Author: Xiang Lin, x.lin@auckland.ac.nz
% Modified by J.C. Mizelle 03/02/2006

% check inputs
if nargin < 2
    error('findfiles.m: You must define entry and directory.');
elseif nargin == 2
    oldDir = pwd;
    if(~isdir(directory))
        error('findfiles.m: No such directory found.');
    end;
    isRecursive = (1==0);
elseif nargin == 3
    oldDir = pwd;
    if(~isdir(directory))
        error('findfiles.m: No such directory found.');
    end;
    if(~(isRecursive == 1))
        error('findfiles.m: isRecursive = 1 if searching subdirectories wanted.');
    end;
    recurse = (1==1);
else
    error('findfiles.m: No more than three inputs ');
end

d = dir(directory);


file_names = {};
f_names = {};
numMatches = 0;
for i=1:length(d)
  
    a_name = d(i).name;
    a_dir = d(i).isdir;
  
    % if the file is not a directory, and there is at least one
    % occurence in the file name or entry = ''
    if(~a_dir & isempty(findstr('.lnk', a_name)))
      
        if(~isempty(findstr(entry,a_name)) | isempty(entry) | (strcmp(entry(1),'~') & isempty(findstr(entry(2:end),a_name))))
            % add the file name to the list.
            numMatches = numMatches + 1;
            file_names{numMatches} = fullfile(directory, a_name);
          
            for ii = 1:length(file_names);
                [temp1 temp2] = fileparts(char(file_names(ii)));
                f_names{numMatches} = strcat(temp1,'\',temp2);
            end


        end;
          
    % if recursive is required and the file is a directory but not '.', '..' and links
    elseif(isRecursive & a_dir & ~strcmp(a_name,'.') & ~strcmp(a_name,'..') & isempty(findstr('.lnk',a_name)))
      
        % solved link problem in windows, not sure in Linux
        file_names = [file_names findfiles(entry, fullfile(directory,a_name), isRecursive)];

        for iii = 1:length(file_names);
        [temp3 temp4] = fileparts(char(file_names(iii)));
        f_names{iii} = strcat(temp3,'\',temp4);
        end
          

        numMatches = length(file_names);
      
    end


end

for i = 1:size(file_names,2);
    [idx idx2] = fileparts(file_names{i}) ;
    f_names{i} = idx2;
end

