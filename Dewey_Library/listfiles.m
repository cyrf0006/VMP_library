function FileList = listfiles(varargin)
%
% listfiles.m--Makes a list of files and returns it as a CELL array.
%
% Note: listfiles.m formerly returned character arrays; some scripts using
% listfiles.m will fail until modified to account for this. Quick fix:
% after running FileList = listfiles('*m'); convert to character array with
% FileList = strvcat(FileList{:})
%
% LISTFILES by itself returns a list of files in the current directory.
%
% LISTFILES(FileSpecifier) returns a list of files according to the string
% found in FileSpecifier. FileSpecifier may include the wildcard character
% "*" and/or path information.
%
% NEW! October, 2003: Can now use '?', '[' and ']' wildcard characters on
% UNIX systems. e.g., FileList = listfiles('cast??.mat'),
% FileList = listfiles('cast[01][0123456789]?.mat')
%
% By default, listfiles.m returns fully-qualified filenames. The optional
% NameType flag can override this behaviour:  
%    LISTFILES(FileSpecifier,'-full') returns the full filenames (default);
%    LISTFILES(FileSpecifier,'-base') returns filenames with no path information.
%
% By default, listfiles.m ignores directories. The FileType flag can be used 
% to override this behaviour:
%    LISTFILES(FileSpecifier,'-file') returns only files (default);
%    LISTFILES(FileSpecifier,'-dir') returns only directories;
%    LISTFILES(FileSpecifier,'-all') returns both files and directories.
%
% If expecting to find only one file or directory matching the file
% specifier, use the '-1' flag. In this case, the variable 'FileList' will
% be returned as a string, rather than a cell array. If more than one file
% is found, an error occurs.
%
% Syntax: FileList = listfiles(<FileSpecifier>,<NameType>,<FileType>,<'-1'>);
%
% e.g., FileList = listfiles('D:\MATLAB\toolbox\matlab\graph2d\*plot*.*')
% e.g., FileList = listfiles('/usr/local/s*','-dir','-base')
% e.g., FileList = listfiles('/usr/local','-dir','-1')
% e.g., FileList = listfiles('/usr/local/s[a-z][ai]*','-all')
%
%%Extended example (cut and paste into your function):
%fileList = listfiles(fullfile(inPath,'*.mat'));
% 
%numFiles = length(fileList);
%       
%for fileCount = 1:numFiles
%
%    disp([mfilename '.m--File ' num2str(fileCount) ' of ' num2str(numFiles)]);
%    thisFile = fileList{fileCount};
%
%end % for

% Developed in Matlab 6.5.0.180913a (R13) on SUN OS 5.8.
% Kevin Bartlett(kpb@hawaii.edu), 2003/12/18, 11:04
%------------------------------------------------------------------------------

% Handle input arguments.

% ...2004-01-16: for consistency with other programs, I am switching to
% using normal "flag" format, with a hypen prepended to each flag (e.g.,
% '-base', instead of 'base'). For backwards compatibility, program will
% still recognise non-hyphened flags.

% ...Default values.
NameType = '-full';
FileType = '-file';
returnSingleFile = 0;

if nargin == 0,
   FileSpecifier = '*';
elseif nargin == 1,
   FileSpecifier = varargin{1};
elseif nargin == 2 | nargin == 3 | nargin == 4
   
   FileSpecifier = varargin{1};
   options = lower(varargin(2:end));
   
   if any(~ismember(options,{'full' 'base' 'file' 'dir' 'all' '-1' '-full' '-base' '-file' '-dir' '-all'}))
      error([mfilename 'm--Unrecognised option.']);
   end % if
   
   if ismember('base',options) | ismember('-base',options)
      NameType = '-base';
   end % if
   
   if ismember('dir',options) | ismember('-dir',options)
      FileType = '-dir';
   elseif ismember('all',options) | ismember('-all',options)
      FileType = '-all';
   end % if

   % ...The hyphen-flag format was used from the start with this flag; no
   % backwards-compatibility need for program to recognise both '-1' and
   % '1' as the same flag.
   if ismember('-1',options)
      returnSingleFile = 1;
   end % if

else   
   error([mfilename 'm--Wrong number of input arguments.']);
end % if

% 2004-01-15: on second thought, this next bit wasn't a good idea.
% listfiles('/home/talpa/') should return empty list. To see contents of
% the directory, user must specify '/home/talpa/*'. 

% % If FileSpecifier represents an existing directory with no file name or
% % wild cards appended to it, then append a '*' wildcard to it.
% if exist(FileSpecifier,'dir') == 7
%    FileSpecifier = fullfile(FileSpecifier,'*');
% end % if

% If path not specified, default to current directory.
[FilePath,name,extension,ver] = fileparts(FileSpecifier);

if isempty(FilePath),
   FilePath = pwd;
   FileSpecifier = fullfile(FilePath,[name,extension,ver]);
end % if

% If the '?' wildcard or '[' or ']' characters have been included in the
% file specifier string, listfiles can use them, providing that this is a
% Unix system.
unixWildCards = '?[]';
doUseUnixCommand = 0;

if any(ismember(FileSpecifier,unixWildCards))
   
   if ~isunix
      disp([mfilename '.m--The wildcard characters ' unixWildCards ' can only be used on UNIX systems.']);
   end % if
   
   doUseUnixCommand = 1;
   
end % if

% If FileSpecifier represents an existing directory with no file name or
% wild cards appended to it, then return an empty list, unless the 'dir' or
% 'all' option is specified, in which case, return the directory name.
if exist(FileSpecifier,'dir') == 7
   
   if strcmp(FileType,'-dir') | strcmp(FileType,'-all') 
      FileList = {};
      FileList{1} = FileSpecifier;
   elseif strcmp(FileType,'-file')
      FileList = {};
   end % if
   
else
   
   % Get list of files conforming to file specifier string. The Unix ls
   % command is not used unless required (and available).
   if doUseUnixCommand == 1
      FileList = listfiles_unix(FileSpecifier,FileType,NameType);
   else
      FileList = listfiles_matlab(FileSpecifier,FilePath,FileType,NameType);   
   end % if
   
   % Sort the file list to get consistent file order.
   FileList = sort(FileList);
   
end % if

% If only a single file/directory expected, return it as a character array.
if returnSingleFile == 1
   
   if length(FileList)>1
      error([mfilename '.m--More than one file or directory found.']);
   elseif length(FileList)==1
      FileList = FileList{1};
   else
      FileList = [];
   end % if
      
end % if

%------------------------------------------------------------------------------
function FileList = listfiles_unix(FileSpecifier,FileType,NameType)
%
% listfiles_unix.m--Uses unix command to get file list for listfiles.m.
%
% Syntax: FileList = listfiles_unix(FileSpecifier,FileType,NameType);

% Developed in Matlab 6.5.1.199709 (R13) Service Pack 1 on SUN OS 5.8.
% Kevin Bartlett(kpb@hawaii.edu), 2004/01/15
%------------------------------------------------------------------------------
[unixStatus,unixResult] = unix(['ls -1d ' FileSpecifier]);

% The ls command executed under unix can fail in two ways. If a
% wild-card combination is not matched, unixStatus will be 1 and
% unixResult will be 'ls: No match.'. If, on the other hand,
% FileSpecifier is an actual filename of a non-existent file (say, for
% example, '/home/talpa/bartlett/silly.mat'), unixStatus will be 0, as
% if there is no problem, but unixResult will be
% '/home/talpa/bartlett/silly.mat not found'. Look for either of these
% failure modes and return an empty file list if encountered.
if unixStatus ~= 0
   FileList = {};      
   return;
else
   
   if strcmp(deblank(unixResult(end-9:end)),'not found') == 1
      FileList = {};
      return;
   else
      FileList = unixResult;
   end % if
   
end % if

% 'ls' returns a \n delimited string of file names. Split the string
% into separate names.
%FileList = split(FileList,sprintf('\n'));
FileList = textscan(unixResult,'%s','delimiter','\n');
FileList = FileList{1};

% Depending on value of FileType, we may want to discard either files or
% directories from the list.
isDirectory = NaN * ones(length(FileList),1);

for FileCount = 1:length(FileList)
   
   if exist(FileList{FileCount},'file')~=2         
      isDirectory(FileCount) = 1;
   else
      isDirectory(FileCount) = 0;
   end % if
   
end % for   

if strcmp(FileType,'-file')      
   FileList = FileList(isDirectory == 0);
elseif strcmp(FileType,'-dir')
   FileList = FileList(isDirectory == 1);
end % if      

% Strip off path information if basenames requested rather than full
% filenames. 
NumFiles = length(FileList);

if strcmp(NameType,'-base')
   
   for FileCount = 1:NumFiles,
      [dummy,baseName,dummy,dummy] = fileparts(FileList{FileCount});
      FileList{FileCount} = baseName;
   end % for
   
end % if   

% Output as a column of filenames to be consistent with output using
% Matlab's "dir" command.
FileList = FileList';

%------------------------------------------------------------------------------
function FileList = listfiles_matlab(FileSpecifier,FilePath,FileType,NameType)
%
% listfiles_matlab.m--Uses Matlab command to get file list for listfiles.m.
%
% Syntax: FileList = listfiles_matlab(FileSpecifier,FilePath,FileType,NameType);

% Developed in Matlab 6.5.1.199709 (R13) Service Pack 1 on SUN OS 5.8.
% Kevin Bartlett(kpb@hawaii.edu), 2004/01/15
%------------------------------------------------------------------------------

eval(['DirStruct = dir(''' FileSpecifier ''');'])
FileList = char(DirStruct.name);

% Exit if no files.
if isempty(FileList)
   FileList = {};
   return;
end % if

IsDir = {};
[IsDir{1:length(DirStruct),1}] = deal(DirStruct.isdir);
IsDir = cat(1,IsDir{:});

% Do not keep '.' or '..' directories.
matchIndex = [strmatch('.',cellstr(FileList),'exact') strmatch('..',cellstr(FileList),'exact')];
keepIndex = setdiff([1:size(FileList,1)],matchIndex);
FileList = FileList(keepIndex,:);
IsDir = IsDir(keepIndex);

% Exit if no files left.
if isempty(FileList)
   FileList = {};
   return;
end % if

% Depending on value of FileType, we may want to discard either files or
% directories from the list.

if strcmp(FileType,'-file')      
   FileList = FileList(IsDir == 0,:);      
elseif strcmp(FileType,'-dir')
   FileList = FileList(IsDir == 1,:);            
else
   % FileType is 'all'; do not discard any.      
end % if      

% Convert to cell array.
if isempty(FileList)
   FileList = {};
else
   FileList = cellstr(FileList);
end % if

% If the user has specified full filenames, prepend path information to the
% list.
NumFiles = length(FileList);

if strcmp(NameType,'-full')
   
   NewFileList = cell(NumFiles,1);
   
   for FileCount = 1:NumFiles,
      %NewFileList{FileCount} = fullfile(FilePath,deblank(FileList(FileCount,:)));
      NewFileList{FileCount} = fullfile(FilePath,FileList{FileCount});
   end % for
   
   FileList = NewFileList;
   clear NewFileList;
   
end % if   

