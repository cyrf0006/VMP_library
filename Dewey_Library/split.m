function stringList = split(inStr,varargin)
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
   GroupedKeepIndex = groupneighbours(KeepIndex);
   
   % Make each chunk of good string data into a separate string. 
   for GroupCount = 1:length(GroupedKeepIndex),
      stringList{GroupCount} = inStr(GroupedKeepIndex{GroupCount}); 
   end % for
   
% ...If NOT ignoring multiple instances of the delimiter, splitting is more complicated:
else

   % Find the non-empty strings.
   KeepIndex = find(~ismember(TotalIndex,delimiterIndex));

   % Group the good string data into chunks.
   GroupedKeepIndex = groupneighbours(KeepIndex);

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
