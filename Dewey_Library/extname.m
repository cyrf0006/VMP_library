function [extendedName] = extname(varargin)
%
% extname.m--Given a basename and extension, extname.m returns the
% "extended name" of a file (basename plus extension). Alternatively,
% extname.m accepts a fully-qualified filename and returns the extended
% name.
%
% The basename and extension can be entered as strings (for single files)
% or as cell arrays of strings (for multiple files). A basename and
% extension of length N will result in N extended names being returned. If
% either the basename or extension is of length one, it will be expanded to
% length N to match the length of the longer of the two.
% 
% The extended name is returned as a cell array of strings, except when the
% basename or extension provided by the user is a simple string. In this
% case, the extended name is returned as a string, if possible (if multiple
% filenames are returned, they will be returned in a cell array of
% strings, regardless).
%
% Syntax: extendedName = extname(baseName,ext)
%   OR
%         extendedName = extname(fullname)
%
% e.g.,   extendedName = extname('readme','.txt')
% e.g.,   extendedName = extname({'readme' 'readmetoo'},'.txt')
% e.g.,   extendedName = extname({'readme' 'readmetoo'},{'.txt' 'dat'})
% e.g.,   extendedName = extname({'readme'},{'.txt' '.dat' '.asc'})
% e.g.,   extendedName = extname('/home/bartlett/readme.txt')
% e.g.,   extendedName = extname({'/home/bartlett/readme.txt' '/home/bartlett/readmetoo.txt'})

% Developed in Matlab 7.0.4.352 (R14) Service Pack 2 on GLNX86.
% Kevin Bartlett (kpb@uvic.ca), 2007-06-12 10:40
%-------------------------------------------------------------------------

% e.g.,   extendedName = extname('readme','.txt')
% e.g.,   extendedName = extname({'readme'},'.txt')
% e.g.,   extendedName = extname('readme',{'.txt'})
% e.g.,   extendedName = extname({'readme'},{'.txt'})
% e.g.,   extendedName = extname({'readme' 'readmetoo'},'.txt')
% e.g.,   extendedName = extname({'readme' 'readmetoo'},{'.txt'})
% e.g.,   extendedName = extname({'readme' 'readmetoo'},{'.txt' 'dat'})
% e.g.,   extendedName = extname({'readme'},{'.txt' '.dat' '.asc'})
% e.g.,   extendedName = extname('readme',{'.txt' '.dat' '.asc'})
% e.g.,   extendedName = extname('/home/bartlett/readme.txt')
% e.g.,   extendedName = extname({'/home/bartlett/readme.txt'})
% e.g.,   extendedName = extname({'/home/bartlett/readme.txt' '/home/bartlett/readmetoo.txt'})

% Run modes for fully-qualified filenames and basename/extension
% combinations:
FQ = 1;
BE = 2;

if nargin == 1
    fqName = varargin{1};
    runMode = FQ;
elseif nargin == 2
    baseName = varargin{1};
    ext = varargin{2};
    runMode = BE;
else
    error([mfilename '.m--Incorrect number of input arguments.']);
end % if

% If running on a fully-qualified filename, extract the basename and
% extension from it.
if runMode == FQ
    [dummy,baseName,ext] = mfileparts(fqName);
end % if

% Convert basename and extension to cells, if necessary.
returnAsStr = 0;
    
if ~iscell(baseName)
    baseName = {baseName};
    returnAsStr = 1;
end % if

if ~iscell(ext)
    ext = {ext};
end % if

% Test to see that baseName and ext are the same length. If they are not
% the same length, then one of them must be of length one.
if length(baseName) ~= length(ext)
    
    if length(baseName)~=1 & length(ext)~=1
        error([mfilename '.m--Incorrect lengths of baseName and/or ext.']);
    end % if
        
end % if

% The length of the variable returned is the length of baseName or ext,
% whichever is larger.
returnLength = max(length(baseName),length(ext));

% Pad the baseName or ext, if required.
if length(baseName) < returnLength
    [baseName{1:returnLength}] = deal(baseName{1});
elseif length(ext) < returnLength
    [ext{1:returnLength}] = deal(ext{1});    
end % if

% Combine the baseName and extension variables.
extendedName = baseName;

for iName = 1:returnLength

    %extendedName{iName} = [extendedName{iName} ext{iName}];
    thisExt = ext{iName};

    % Prepend a period to the extension if one was not specified.
    if ~strcmp(thisExt(1),'.')
        thisExt = ['.' thisExt];
    end % if
    
    extendedName{iName} = [extendedName{iName} thisExt];

end % for

% Return as string if requested and if possible.
if returnAsStr == 1 & returnLength == 1
    extendedName = extendedName{1};
end % if

