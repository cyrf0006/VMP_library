function [] = mkpackage(FileSpec,varargin)
%
% mkpackage.m--Makes a software package by checking the dependencies of the
% specified m-file(s).
%
% M-files are specified with the input argument FileSpec. FileSpec can be a
% fully-qualified filename (e.g., '/home/talpa/bartlett/myfunc.m'), the
% filename of an m-file in the current directory (e.g., 'myfunc.m'), a
% file-specifier string containing the '*' wildcard (e.g., 'myf*.m') or the
% name of a function on Matlab's search path (e.g., 'myfunc').
%
% MKPACKAGE(FileSpec) copies the specified m-files to a directory named
% "mkpackageDir", created by mkpackage.m in the current working directory.
% The dependent files are copied to a ".../mkpackageDir/private/"
% sub-directory.
%
% MKPACKAGE(FileSpec,'-file') also makes a package, but does so by appending
% all the code for the dependent m-files to the specified m-file. This only
% works when FileSpec corresponds to a single m-file. The original m-file
% will NOT be overwritten.
%
% Note that the user-specified m-file(s) must be functional in its current
% directory location for a complete package to be created. e.g., if you 
% make the request mkpackage('myprogram.m'), and myprogram.m makes use of
% the function myfunction.m, then myfunction.m must be on the Matlab path at 
% the time mkpackage.m is run or it will not appear in the list of dependent 
% functions and be included in the package.
%
% N.B., the package created by mkpackage.m will frequently contain extra
% files in addition to the ones actually needed; this is due to a bug in
% Matlab's depfun.m function.
% 
% Syntax: mkpackage(FileSpec,<'-file'>)
%
% e.g.,   mkpackage('*.m')
%
% e.g.,   mkpackage('/home/talpa/bartlett/matlab/FileOps/newmf.m','-file')

% Developed in Matlab 6.5.0.180913a (R13) on SUN OS 5.8.
% Kevin Bartlett(kpb@hawaii.edu), 2003/01/30, 16:34
%------------------------------------------------------------------------------

% Development note: this copy of mkpackage.m is itself a package,
% containing the functions split.m, listfiles.m and groupneighbours.m.
% These have been renamed mkpackage_split.m, mkpackage_listfiles.m and
% mkpackage_groupneighbours.m within this file in order to avoid errors
% when running mkpackage on these functions themselves. If split, listfiles
% or groupneighbours is modified in the future, they will have to be
% changed manually within mkpackage.m as well.

% Handle input arguments.
FileFlag = 0;

if nargin == 2
   
   if strcmp(lower(varargin{1}),'-file') == 1
      FileFlag = 1;
   else
      error([mfilename '.m--Unrecognised value for FileFlag.'])
   end % if
   
end % if

% List the files specified by the user.
FileList = mkpackage_listfiles(FileSpec,'full');
FileList = strvcat(FileList{:});

% ...If no files were found using FileSpec, the user may have specified an
% m-file on Matlab's search path.
if isempty(FileList)
   FileList = which(FileSpec);      
end % if

if isempty(FileList)
   error([mfilename '.m--Must specify m-file(s) on Matlab''s search path.'])
end % if

% Creating a package inside an m-file is only allowed when a single m-file
% has been specified by the user.
if FileFlag == 1 & size(FileList,1)>1
   error([mfilename '.m--Cannot create a single-file package for multiple m-files.'])   
end % if

% Want a list of just the basenames for use later.
NumFiles = size(FileList,1);
BaseList = [];

for FileCount = 1:NumFiles
   
   ThisFile = deblank(FileList(FileCount,:));
   [dummy,ThisBaseName,dummy,dummy] = fileparts(ThisFile);
   BaseList = strvcat(BaseList,ThisBaseName);
   
end % for

% Make the new directory to contain the package.
mkpackageDirName = 'mkpackageDir';
mkpackageDir = fullfile(pwd,mkpackageDirName);

if exist(mkpackageDir,'dir')==7
   
    %error([mfilename '.m--Directory ' mkpackageDir ' already exists. Delete it and start over.'])

    disp([mfilename '.m--Directory ' mkpackageDir ' already exists.']);
    resp = input('Delete it and start over? Y/N:   ','s');
   
    if ~isempty(strmatch(lower(resp),'yes'))   
        
        [rmdir_success,rmdir_mssg,rmdir_mssgID] = rmdir(mkpackageDir,'s');
        
        if rmdir_success == 0
            error([mfilename '.m--Failed to remove existing directory: ' mkpackageDir '.']);
        end % if
        
    else
        disp([mfilename '.m--Execution halted; package not created.']);
        return;
    end % if
   
end % if

% Mkdir uses only 2 output arguments on linux version of matlab.
%[mkdir_success,mkdir_mssg,mkdir_mssgID] = mkdir(pwd,mkpackageDirName);
[mkdir_success,mkdir_mssg] = mkdir(pwd,mkpackageDirName);

if mkdir_success ~= 1
   error([mfilename '.m--Error creating directory: ' mkdir_mssg]);
end % if

% Make a "private" directory for holding the dependent files.
privateDirName = 'private';
privateDir = fullfile(mkpackageDir,privateDirName);

[mkdir_success,mkdir_mssg] = mkdir(mkpackageDir,privateDirName);

if mkdir_success ~= 1
   error([mfilename '.m--Error creating directory: ' mkdir_mssg]);
end % if

%[SourceDir,dummy,dummy,dummy] = fileparts(FileSpec);

% Copy all the specified m-files into the mkpackage directory.
for FileCount = 1:NumFiles
   
   RootFile = deblank(FileList(FileCount,:));
   [SourceDir,RootBaseName,RootExt,dummy] = fileparts(RootFile);
   RootFileName = [RootBaseName RootExt];
   RootFullName = fullfile(SourceDir,RootFileName);
   
   % Copy the current file itself to the package directory.
   % N.B., should really use copyfile.m, but it isn't working on my
   % computer.
   %[copyfile_success,copyfile_mssg,copyfile_mssgID] = copyfile(CurrFile,mkpackageDir);
   
   % Okay, it is working now (UVic, 2005-06-20).
   %[status,result] = unix(['cp ' RootFullName ' ' mkpackageDir]);   
   [copyfile_success,copyfile_mssg,copyfile_mssgID] = copyfile(RootFullName,mkpackageDir);
   
   %if status ~= 0
   if copyfile_success ~= 1
      error([mfilename '.m--Error copying file: ' result]);   
   end % if
   
end % for

% If making a single-file package, open the file for writing.
OutFid = NaN;

if FileFlag == 1

   RootFile = deblank(FileList(FileCount,:));
   [dummy,RootBaseName,RootExt,dummy] = fileparts(RootFile);
   RootFileName = [RootBaseName RootExt];

   OutFile = fullfile(mkpackageDir,RootFileName);
   OutFid = fopen(OutFile,'at');
   
   if OutFid < 0
      error([mfilename '.m--Failed to open output file for writing.']);   
   end % if
   
   fseek(OutFid,0,'eof');
   
end % if

% Find the dependencies for the specified m-files.
for FileCount = 1:NumFiles
   
   disp([mfilename '.m--Making package; m-file ' num2str(FileCount) ' of ' num2str(NumFiles)]);
   RootFile = deblank(FileList(FileCount,:));
   [dummy,RootBaseName,RootExt,dummy] = fileparts(RootFile);
   RootFileName = [RootBaseName RootExt];
   
   % Find the dependents of this m-file.
   DepList = mydeps(RootFile);
   
   % ...Only want to copy files/code once.
   AlreadyCopied = BaseList;
   
   for DepCount = 1:length(DepList)
      
      CurrDepFile = DepList{DepCount};
      [dummy,BaseName,ext,dummy] = fileparts(CurrDepFile);
      CurrDepFileName = [BaseName ext];
      
      % Only copy if not copied before.
      if isempty(strmatch(BaseName,AlreadyCopied,'exact'))
         
         if FileFlag == 0
            
            %[status,result] = unix(['cp ' CurrDepFile ' ' privateDir]);
            [copyfile_success,copyfile_mssg,copyfile_mssgID] = copyfile(CurrDepFile,privateDir);

            %if status ~= 0
            if copyfile_success ~= 1
               error([mfilename '.m--Error copying file: ' result]);   
            end % if
            
         else
            
            % Open the current dependent file for reading.
            InFid = fopen(CurrDepFile,'rt');
   
            if InFid < 0
               error([mfilename '.m--Failed to open file ' CurrDepFile ' for reading.']);   
            end % if
   
            fprintf(OutFid,'\n\n%s\n','%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');            

            % Append the code from the dependent file to the output file.
            while ~feof(InFid)            
               tline = fgetl(InFid);
               fprintf(OutFid,'%s\n',tline);
            end % while
            
            fclose(InFid);
               
         end % if
         
         AlreadyCopied = strvcat(AlreadyCopied,BaseName);
         
      end % if
                  
   end % for each dependency
   
end % for each specified m-file.

if FileFlag == 1
   fclose(OutFid);
end % if

% Remove the "private" directory if it is empty.
PrivDirFiles = dir(privateDir);
NumPrivFiles = length(PrivDirFiles)-2; % (dir returns the '.' and '..' directories).

if NumPrivFiles < 1
   
   % Should use Matlab's rmdir.m, but it doesn't exist on the version of
   % Matlab on the laptop.
   %[status,result] = unix(['rmdir ' privateDir]);

   % 2005-06-20 Okay, works now (UVic, release 14).
   [rmdir_success,rmdir_mssg,rmdir_mssgID] = rmdir(privateDir);
   
end % if

disp([mfilename '.m--Package created.']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function groupMembers = mkpackage_groupneighbours(IndexIn)
%
% groupneighbours.m--Takes a vector of numbers (usually an index to another
% vector) and groups them into contiguous segments.
%
% Groupneighbours.m returns groupMembers, a cell array, each cell of which
% contains a group of next-door neighbours. A vector element with no other
% element adjacent to it will appear in its own cell.  
%
% Syntax: groupMembers = groupneighbours(IndexIn);
%
% e.g., groupMembers = groupneighbours([20 30 40 41 42 43 50 60 61 70 71 72 80])

% Kevin Bartlett (bartlettk@dfo-mpo.gc.ca) 11/1999
%------------------------------------------------------------------------------
% Test for development: groupMembers = groupneighbours([20 30 40 41 42 43 50 52 55 57 60 61 70 71 72 80 81])

% Modification, 2003-03-14: program very slow on large vectors. Speed up by
% looping over groups of neighbours, rather than individual points.

if isempty(IndexIn),
   groupMembers = [];
   return;
end % if

% Convert to column vector if necessary.
indexSize = size(IndexIn);

if ~any(indexSize==1)
   error([mfilename '.m--IndexIn must be a vector.'])
end % if

if indexSize(1)>indexSize(2)
   IndexIn = IndexIn';
   wasFlipped = 1;
else
   wasFlipped = 0;
end % if

% Sort the vector of integers.
IndexIn = sort(IndexIn);

% Test that there are no repeated values in the input vector.
if length(unique(IndexIn)) ~= length(IndexIn),
   error([mfilename '.m--IndexIn must not have repeated values.'])
end % if

% Find which points in the vector are NOT adjacent to each other. These
% points mark the divisions ("fences") between groups of neighbours.
FenceIndex = find(diff(IndexIn)>1);

% The start of each group follows each division.
GroupStartIndex = FenceIndex + 1;

% ...The first element in the input vector is always the start of its group.
GroupStartIndex = [1 GroupStartIndex];

% The end of each group is given by the fence index, with the addition of
% the last point in the input vector, which is always the end of its group.
GroupEndIndex = [FenceIndex length(IndexIn)];

% Initialise output variable.
NumGroups = length(GroupStartIndex);
groupMembers = cell(1,NumGroups);

% Find the members of each group.
for GroupCount = 1:NumGroups
   groupMembers{GroupCount} = IndexIn([GroupStartIndex(GroupCount):GroupEndIndex(GroupCount)]);
end % for

% Flip the index to have the same orientation as the input index.
if wasFlipped == 1
   groupMembers = groupMembers';
end % if


function FileList = mkpackage_listfiles(varargin)
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
FileList = mkpackage_split(FileList,sprintf('\n'));

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
% listfiles_unix.m--Uses Matlab command to get file list for listfiles.m.
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stringList = mkpackage_split(inStr,varargin)
%
% split.m--Mimics the behaviour of Perl's split function. Split.m splits up
% the input string "inStr" by some delimiter. 
%
% SPLIT(STR)--splits the input string by whitespace characters.
%
% SPLIT(STR,DELIMITER)--splits the input string by the delimiter string
% given by DELIMITER.
%
% SPLIT(STR,<'-d',DELIMITER>,<'-g',GLUECHAR>)--splits the input string by
% the delimiter string given by DELIMITER. Portions of the input string
% enclosed by the character given by GLUECHAR are not split, even if they
% contain the delimiter character. 
%
% SPLIT(...,<'-i',doIgnore>)--By default, split.m ignores multiple adjacent
% instances of the delimiter, treating them as a single delimiter
% character. This prevents 'Hi  there' from being split into 4 strings (2
% of them empty). Specify this behaviour explicitly with a value of
% doIgnore of 1; override it with a value of 0.
% 
% Note that specifying the delimiter to be ' ' will not give the same
% behaviour as the default whitespace delimiter, as in this case only
% spaces will be considered, leaving other whitespace characters (such as
% tab) alone. Explicitly specify a whitespace delimiter with a NaN.
%
% Strings are returned as separate cells in the cell array "stringList".
%
% e.g.,   stringList = split('This string has 6 space-delimited words')
% e.g.,   stringList = split('This#string#has#6###pound-delimited#words','#')
% e.g.,   stringList = split('#This#string#has#13###pound-delimited#words,#but#3#are#empty','-d','#','-i',0)
% e.g.,   str = sprintf('%s\t%s','Hi','there'); stringList = split(str,' ')
% e.g.,   str = sprintf('%s\t%s','Hi','there'); stringList = split(str,NaN)
% e.g.,   stringList = split('1001 "concrete block" anchor -1000 2 1 1','-g','"')

% Developed in Matlab 6.1.0.450 (R12.1) on SUN OS 5.8.
% Kevin Bartlett(bartlett@soest.hawaii.edu), 2001/10/22, 15:13
%------------------------------------------------------------------------------

% Handle input arguments.

% ...Default values:
delimiter = NaN; % default value denotes whitespace
doIgnore = 1;
glueChar = '';

if nargin == 1,
   delimiter = NaN;
elseif nargin == 2,
    
    delimiter = varargin{1};
    
    if length(varargin)>1
        varargin = varargin(2:end);
    else
        varargin = {};
    end % if

end % if

if rem(length(varargin),2) ~= 0
    error([mfilename '.m--Incorrect number of input arguments']);
end % if
    
for argCount = 1:2:length(varargin)
    
    thisFlag = varargin{argCount};
    thisVal = varargin{argCount+1};
    
    if strcmp(thisFlag,'-g')
        glueChar = thisVal;
    elseif strcmp(thisFlag,'-i')
       doIgnore = thisVal;
    elseif strcmp(thisFlag,'-d')
        delimiter = thisVal;
    else
        error([mfilename '.m--Unrecognised flag: ''' thisFlag '''']);
    end % if
    
end % for

% Test input arguments.
IsBadDelimiter = (length(delimiter)~=1) | (~isstr(delimiter) & ~isnan(delimiter));

if IsBadDelimiter == 1,
   error([mfilename '.m--Bad ''delimiter'' value. Delimiter can be NaN or a single character.'])
end % if

IsBadDoIgnore = (length(doIgnore)~=1) | (~ismember(doIgnore,[0 1]));

if IsBadDoIgnore == 1,
   error([mfilename '.m--Bad ''doIgnore'' value. doIgnore must be 1 or 0.'])
end % if

% "Glue" characters must come in pairs to work.
glueCharIndex = strfind(inStr,glueChar);

if rem(length(glueCharIndex),2)~=0
    error([mfilename '.m--''Glue'' characters must come in pairs.'])
end % if

TotalIndex = 1:length(inStr);

% Split the input string by the delimiter.
if isnan(delimiter), % (if splitting on whitespace characters)
   delimiterIndex = find(isspace(inStr));
else                 
   delimiterIndex = strfind(inStr,delimiter);      
end % if

% Handle case in which the input string contains only delimiters.
if length(delimiterIndex)==length(inStr)
   stringList = [];
   return;
end % if

% Handle case in which no the input string contains no delimiters.
if length(delimiterIndex)==0,
   stringList = {inStr};
   return;
end % if

% Delimiters that appear within a "glued" string should be ignored.
for glueCount = 1:2:length(glueCharIndex)
    delimiterIndex = setdiff(delimiterIndex,glueCharIndex(glueCount):glueCharIndex(glueCount+1));
end % for

% ...If ignoring multiple instances of the delimiter, splitting is fairly straightforward:
if doIgnore == 1,

   % Group the good string data into chunks.
   KeepIndex = find(~ismember(TotalIndex,delimiterIndex));
   GroupedKeepIndex = mkpackage_groupneighbours(KeepIndex);
   
   % Make each chunk of good string data into a separate string. 
   for GroupCount = 1:length(GroupedKeepIndex),
      stringList{GroupCount} = inStr(GroupedKeepIndex{GroupCount}); 
   end % for
   
% ...If NOT ignoring multiple instances of the delimiter, splitting is more complicated:
else

   % Find the non-empty strings.
   KeepIndex = find(~ismember(TotalIndex,delimiterIndex));

   % Group the good string data into chunks.
   GroupedKeepIndex = mkpackage_groupneighbours(KeepIndex);

   % Find the starting index of each chunk of string data.
   ChunkStartIndex = [];

   for ChunkCount = 1:length(GroupedKeepIndex),
      CurrChunk = GroupedKeepIndex{ChunkCount};
      ChunkStartIndex(ChunkCount) = CurrChunk(1);
   end % for

   % Make an index of all points following a delimiter. All of these points mark the beginning
   % of an output string, but some of the output strings may be empty.
   PostdelimiterIndex = delimiterIndex + 1;

   % If the first character of the input string is a delimiter, then the first output string
   % is an empty one. Otherwise, the first output string is the first good string chunk.
   OutCount = 1;
   ChunkCount = 0;

   if delimiterIndex(1)==1,    
      stringList{OutCount} = []; 
   else
      ChunkCount = ChunkCount + 1;
      stringList{OutCount} = inStr(GroupedKeepIndex{ChunkCount});      
   end % if

   % Go through the index of output string beginnings and assemble the output string list.
   for PostDelimiterCount = 1:length(PostdelimiterIndex),

      OutCount = OutCount + 1;
      CurrIndex = PostdelimiterIndex(PostDelimiterCount);

      % If the current string beginning is also the start of one of the non-empty strings,
      % insert the string into the output string list. Otherwise, insert an empty string.
      if ismember(CurrIndex,ChunkStartIndex)==1,         
         ChunkCount = ChunkCount + 1;
         stringList{OutCount} = inStr(GroupedKeepIndex{ChunkCount});      
      else
         stringList{OutCount} = []; 
      end % if

   end % for each output string beginning.

end % if ignoring multiple instances

% If "glue" characters used to keep sections of strings together, we don't
% want them to appear in the output.
if length(glueCharIndex)>0
    
    for strCount = 1:length(stringList)
        stringList{strCount} = strrep(stringList{strCount},glueChar,'');
    end % for
    
end % if

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DepList = mydeps(RootName)
%
% mydeps.m--Locates dependent, user-defined functions of a specified "root"
% function, allowing  the developer to assemble a self-contained software
% package for distribution.
%
% The root function cannot contain a relative path.
%
% The root function being explored may evaluate strings in order to call
% other functions; mydeps.m is unable to detect such dependencies.
%
% N.B., In addition to authentic, dependent functions, mydeps.m will
% usually also return extra files that would never be called if the root
% function were actually evaluated (see "help dep_fun").
%
% Syntax: DepList = mydeps(RootName)
%
% e.g.,   DepList = mydeps('designmooring')

% Developed in Matlab 6.1.0.450 (R12.1) on Win98. Kevin
% Bartlett(kpb@hawaii.edu), 2002/12/03, 16:12
%------------------------------------------------------------------------------

% Crash problem under java session. If "matlab -nojvm" session is used, you see 
% the following error message:
% % Warning: File: /usr/local/matlab/toolbox/matlab/iofun/aviinfo.m Line: 311 Column: 1
% % Variable 'length' has been previously used as a function name.
% % (Type "warning off MATLAB:mir_warning_variable_used_as_function" to suppress this warning.)
% % In /home/talpa/bartlett/matlab/FileOps/mydeps.m (sniff) at line 955
% %  In /home/talpa/bartlett/matlab/FileOps/mydeps.m (depfun) at line 370
% %  In /home/talpa/bartlett/matlab/FileOps/mydeps.m at line 75
% 
% Note, however, that "warning off MATLAB:mir_warning_variable_used_as_function" 
% does not solve the problem of Matlab crashing. ==> Don't run mydeps.m
% under a java session unless you're feeling lucky.
%
% ...2005-06-08: Now on Matlab Release 14, find that nojvm session results
% in error message. ==> use java session instead. Also, '-nographics' flag
% no longer supported for depfun.m, but new '-quiet' flag makes use of
% modified "mod_depfun.m" no longer necessary.
JavaInfo = version('-java');
JavaEnabled = isempty(strfind(JavaInfo,'not enabled'));

% if JavaEnabled
%    error([mfilename '.m--This program will likely crash Matlab if run under a java session. Use a "matlab -nojvm" session instead.']);
% end % if

if JavaEnabled == 0
   error([mfilename '.m--This program will likely crash Matlab if not run under a java session. Don''t use a "matlab -nojvm" session.']);
end % if

% % Running mydeps.m on laptop running Matlab release 12.1 results in DepList
% % being empty, even for m-files known to have dependencies. Get around this
% % by using Matlab's own depfun function instead of my modified version
% % (modified version results in less rubbish echoed to screen).
% if str2num(version('-release'))<13
%    %error([mfilename '.m--This program will not work under this version of Matlab. Must be at least release 13.']);
%    UseOrigDepFun = 1;
% else
%    UseOrigDepFun = 0;   
% end % if
UseOrigDepFun = 1;

% Looking for only at those directories containing non-MathWorks code.
MatPath = matlabroot;

if exist(RootName,'file')~=2
   error([mfilename '.m--M-file does not exist.']);
end % if

FullName = which(RootName);
[RootPath,dummy,dummy,dummy] = fileparts(FullName);
PrivDir = fullfile(RootPath,'private','');

if exist(PrivDir,'dir')==7
   warning([mfilename '.m--"private" directory exists: contents not examined.']);
end % if   

% Search for dependencies.
DepList = [];

% % Prefer to use my modified depfun.m, but use original if necessary.
% if UseOrigDepFun == 0
%    [tracelist,builtins,classnames,probfiles,probsymbols,evalstrs,calledfrom,jclass] = mod_depfun(FullName,'-nographics');
% else
%    [tracelist,builtins,classnames,probfiles,probsymbols,evalstrs,calledfrom,jclass] = depfun(FullName,'-nographics');
% end % if
 [tracelist,builtins,classnames,probfiles,probsymbols,evalstrs,calledfrom,jclass] = depfun(FullName,'-quiet');
%keyboard

% Exclude standard Matlab routines from consideration (only interested in custom
% functions).
MatchIndex = strmatch(MatPath,tracelist);
NonMatchIndex = setdiff([1:length(tracelist)],MatchIndex);
CurrList = strvcat(tracelist{NonMatchIndex});

% Exclude the root function itself from the list of dependencies.
RootIndex = strmatch(FullName,CurrList,'exact');
KeepIndex = setdiff(1:size(CurrList,1),RootIndex);
CurrList = CurrList(KeepIndex,:);

DepList = strvcat(DepList,CurrList);  

% Output a cell array.
if isempty(DepList)
   DepList = {};
else
   DepList = unique(cellstr(DepList));
end % if

