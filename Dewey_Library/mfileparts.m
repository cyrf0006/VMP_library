function [pathStr,baseName,ext,versn] = mfileparts(fileList)
%
% mfileparts.m--Runs Mathworks' fileparts.m function on multiple filenames
% (fileparts.m only works on one file at a time). Result is returned in the
% form of cell arrays.
%
% If input argument is a single filename in the form of a string, rather
% than a cell array, mfileparts.m acts just like fileparts.m and
% returns strings rather than cell arrays.
%
% Syntax: [pathStr,baseName,ext,versn] = mfileparts(fileList)
%
% e.g.,   fileList = {'/home/bartlett/matlab/FileOps/where.m',...
%                     '/home/bartlett/matlab/bearings/count_rots.m',...
%                     '/home/bartlett/matlab/bearings/plot_rots.m',...
%                     '/home/bartlett/matlab/colours/getcolours.m',...
%                     '/home/bartlett/matlab/filters/decimate_me.m'};
%         [pathStr,baseName,ext,versn] = mfileparts(fileList)
%
% e.g.,   [pathStr,baseName,ext,versn] = mfileparts('/home/bartlett/matlab/FileOps/where.m')

% Developed in Matlab 7.3.0.298 (R2006b) on GLNX86.
% Kevin Bartlett (kpb@uvic.ca), 2007-04-27 16:08
%-------------------------------------------------------------------------

if iscell(fileList)
    isCellInput = 1;
elseif isstr(fileList)
    isCellInput = 0;
    fileList = cellstr(fileList);
else
    error([mfilename '.m--Input must be a cell array of strings or a single string.']);
end % if

nFiles = length(fileList);

[pathStr{1:nFiles}] = deal('');
baseName = pathStr;
ext = pathStr;
versn = pathStr;

for iFile = 1:nFiles    
    thisFile = fileList{iFile};
    [thisPathStr,thisBaseName,thisExt,thisVersn] = fileparts(thisFile);
    pathStr{iFile} = thisPathStr;
    baseName{iFile} = thisBaseName;
    ext{iFile} = thisExt;
    versn{iFile} = thisVersn;        
end % for

% If input was a string, make output just like fileparts.m.
if ~isCellInput
    pathStr = pathStr{1};
    baseName = baseName{1};
    ext = ext{1};
    versn = versn{1};
end % if
