%% inifile_with_instring
% Legacy function replaced by "setupstr"
%%
% <latex>\index{Depreciated!inifile\_with\_instring}</latex>
%
%%% Description
% This is a legacy function.  Please replace with the function 'setupstr'.
%
% Complete documentation for this function can be found in the ODAS Library
% Manual v3.0.  Alternatively, the documentation can be viewed by editing
% this function with the command: 
%
%    >> edit inifile_with_instring
%

% inifile_with_instring
% read data from a standard ini string
%
% |varargout = inifile_with_instring(varargin)|
%
% * [in] varargin - see Usage
%
% * [out] retval varargout - see Usage
% 
%
% *Usage:*
%
%    - [readsett,result] = INIFILE_WITH_INSTRING(mystring,'readstring',keys)
%        Reads the values of the 'keys' where 'readsett' is a cell array of
%        strings and/or numeric values of the keys. If any of the keys
%        is not found, the default value is returned (if given in the
%        5-th column of the keys parameter). 'result' is a cell array of
%        strings - one for each key read; empty if OK, error/warning
%        string if error; in both cases an empty string is returned in
%        readsett{i} for the i-th key if error. 'mystring' is a string that
%        contains keys/values/sections in 'ini' format.
%
%    - [keys,sections,subsections] = INIFILE_WITH_INSTRING(mystring,'readallfromstring')
%        Reads entire 'ini' style string and returns all the sections, subsections
%        and keys found.
%
% Notes on the keys cell array given as an input parameter:
% Cell array of STRINGS; either 3, 4, or 5 columns. 
% Each row has the same number of columns. The columns are:
%
%        - 'section':      section name string (the root is considered if
%                        empty)
%        - 'subsection':   subsection name string (the root is considered
%                        if empty)
%        - 'key':          name of the field to write/read from (given as
%                        a string).
%        - value:          (optional) STRING or NUMERIC value (scalar or
%                        matrix) to be written to the
%                        ini file in the case of 'write' operation OR
%                        conversion CHAR for read operation:
%                        'i' for integer, 'd' for double, 's' or
%                        '' or not given for string (default).
%        - defaultValue:   (optional) string or numeric value (scalar or
%                        matrix) that is returned when the key is not
%                        found or an empty value is found
%                        when reading ('read' operation).
%                        If the defaultValue is not given and the key
%                        is not found, an empty value is returned.
%                        It MUST be in the format as given by the
%                        value, e.g. if the value = 'i' it must be
%                        given as an integer etc.
%
% *NOTES:*
%
% All but the parameter value are NOT case sensitive. Given keys and values
% will be trimmed (leading and trailing spaces will be removed).
% Any duplicates (section, subsection, and keys) are ignored. Empty
% section and/or subsection can be given as an empty string, '',
% but NOT as an empty matrix, [].
% 
%
% *EXAMPLES:*
% 
%  For the following examples, let's assume a string called 'myinistring'
%  contains the following:
%
%  [mysection]
%  mykey=Beethoven
%
%  [myothersection]
%  myotherkey=3.14159265
%  
%
%  To read individual keys:
%  myvalue = inifile_with_instring(myinistring, 'readstring', ...
%               {'mysection','','mykey'});
%
%  myvalue is a 1x1 cell and contains 'Beethoven'
%
%  To read out multiple values:
%  readKeys = {'mysection','','mykey'; 'myothersection','','myotherkey';};
%  myvalues = inifile_with_instring(myinistring,'readstring',readKeys);
%
%  myvalues is a 2x1 cell and contains:
%    'Beethoven'
%    '3.14159265'
%
%  Or, we can just read the entire ini string with:
%  [keys,sections,subsections] = inifile_with_instring(myinistring,
%               'readallfromstring');
%
%  keys, a 2x4 cell, will contain:
%    'mysection'          ''    'mykey'         'Beethoven' 
%    'myothersection'     ''    'myotherkey'    '3.14159265'
%    
%  sections, 2x1 cell, will contain:
%    'mysection'
%    'myothersection'
  
%
% 
% ACKNOWLEDGEMENTS: based on INIFILE - original author:
% 
% (c) Primoz Cermelj, Slovenia
% Contact: primoz.cermelj@gmail.com
% 2009-10-23 AWS added parsing of ini-style string instead of reading from a file
% * 2011-09-01 (AWS) added documentation tags for matlab publishing
% * 2012-04-30 (WID) modified parsing algorithm to allow comments mid-line
% * 2012-09-03 (WID) added a warning statement - it is now depreciated

function varargout = inifile_with_instring(varargin)

global NL_CHAR;

warning('ODAS:depreciatedFunction:inifile_with_instring', ...
    ['The function "inifile_with_instring" has been depreciated. \n' ...
     'Please use the alternative function, "setupstr".  Note this' ...
     'function is used differently.  Type "help setupstr" for more '...
     'information.  Examples of this function in use can be found '...
     'throughout the ODAS library.']);

% Checks the input arguments
if nargin == 0
    disp('INIFILE_WITH_INSTRING v1.0.0');
    disp('Copyright (c) 2003-2007 Primoz Cermelj');
    disp('This is FREE SOFTWARE');
    disp('Type <help inifile> to get more help on its usage');
    return
elseif nargin < 2
    error('Not enough input arguments');
end

fileName = varargin{1};
myString = varargin{1}; %is used by functions that read from string
operation = varargin{2};

if (strcmpi(operation,'read')) | (strcmpi(operation,'deletekeys'))
    if nargin < 3
        error('Not enough input arguments.');
    end
    if ~exist(fileName)
        error(['File ' fileName ' does not exist.']);
    end
    keys = varargin{3};
    [m,n] = size(keys);
    if n < 3
        error('Keys argument must have at least 3 columns for read operation');
    end
    for ii=1:m
        if isempty(keys(ii,3)) | ~ischar(keys{ii,3})
            error('Empty or non-char keys are not allowed.');
        end
    end
elseif (strcmpi(operation,'write'))
    if nargin < 3
        error('Not enough input arguments');
    end
    keys = varargin{3};
    if nargin < 4 || isempty(varargin{4})
        style = 'plain';
    else
        style = varargin{4};
        if ~(strcmpi(style,'plain') | strcmpi(style,'tabbed')) | ~ischar(style)
            error('Unsupported style given or style not given as a string');
        end
    end
    [m,n] = size(keys);
    if n < 4
        error('Keys argument requires 4 columns for write operation');
    end        
    for ii=1:m
        if isempty(keys(ii,3)) | ~ischar(keys{ii,3})
            error('Empty or non-char keys are not allowed.');
        end
    end
elseif (strcmpi(operation,'readstring'))    
    
    if isempty(myString)
        error('ini string cannot be empty');
    end

    if nargin < 3
        error('Not enough input arguments.');
    end
    keys = varargin{3};
    [m,n] = size(keys);
    if n < 3
        error('Keys argument must have at least 3 columns for read operation');
    end
    for ii=1:m
        if isempty(keys(ii,3)) | ~ischar(keys{ii,3})
            error('Empty or non-char keys are not allowed.');
        end
    end
elseif (strcmpi(operation,'readall'))   
    %
elseif (strcmpi(operation,'readallfromstring'))
    %
elseif (~strcmpi(operation,'new'))
    error(['Unknown inifile operation: ''' operation '''']);
end
if nargin >= 3
    for ii=1:m
        for jj=1:3
            if ~ischar(keys{ii,jj})
                error('All cells from the first 3 columns must be given as strings, even the empty ones.');
            end
        end
    end
end


% Sets the new-line character/string
if ispc
    NL_CHAR = '\r\n';
elseif isunix
    NL_CHAR = '\n';
else
    NL_CHAR = '\r';
end

readsett = [];
result = [];


%----------------------------
% CREATES a new, empty file (rewrites an existing one)
%----------------------------
if strcmpi(operation,'new')
    fh = fopen(fileName,'w');
    if fh == -1
        error(['File: ''' fileName ''' can not be (re)created']);
    end
    fclose(fh);
    return

%----------------------------
% READS the whole data (all keys)
%----------------------------    
elseif (strcmpi(operation,'readall'))
    [keys,sections,subsections] = readallkeys(fileName);
    varargout(1) = {keys};
    varargout(2) = {sections};
    varargout(3) = {subsections};
    return    
    
%----------------------------
% READS the whole data (all keys) from string
%----------------------------    
elseif (strcmpi(operation,'readallfromstring'))
    [keys,sections,subsections] = readallkeysfromstring(myString);
    varargout(1) = {keys};
    varargout(2) = {sections};
    varargout(3) = {subsections};
    return    

%----------------------------
% READS key-value pairs out
%----------------------------    
elseif (strcmpi(operation,'read'))
    result = cell(m,1);
    if n >= 4
        conversionOp = keys(:,4);       % conversion operation: 'i', 'd', or 's' ('') - for each key to be read
    else
        conversionOp = cellstrings(m,1);
    end
    if n < 5
        defaultValues = cellstrings(m,1);
    else
        defaultValues = keys(:,5);
    end
    readsett = defaultValues;
    keysIn = keys(:,1:3);
    [secsExist,subsecsExist,keysExist,readValues,so,eo] = findkeys(fileName,keysIn);    
    ind = find(keysExist);
    % For those keys that exist but have empty values, replace them with
    % the default values
    if ~isempty(ind)
        ind_empty = zeros(size(ind));
        for kk = 1:size(ind,1)
            ind_empty(kk) = isempty(readValues{ind(kk)});
        end
        ind(find(ind_empty)) = [];
        readsett(ind) = readValues(ind);
    end 
    % Now, go through all the keys and do the conversion if the conversion
    % char is given
    for ii=1:m
        if ~isempty(conversionOp{ii}) & ~strcmpi(conversionOp{ii},'s')
            if strcmpi(conversionOp{ii},'i') | strcmpi(conversionOp{ii},'d')
                if ~isnumeric(readsett{ii})
                    readsett{ii} = str2num(readsett{ii});
                end
                if strcmpi(conversionOp{ii},'i')
                    readsett{ii} = round(readsett{ii});
                end
                if isempty(readsett{ii})
                    result{ii} = [num2str(ii) '-th key ' keysIn{ii,3} 'or given defaultValue could not be converted using ''' conversionOp{ii} ''' conversion'];
                end
            else
                error(['Invalid conversion char given: ' conversionOp{ii}]);
            end
        end
    end
    varargout(1) = {readsett};
    varargout(2) = {result};
    return
    
%Read from string instead of file
elseif (strcmpi(operation,'readstring'))
    result = cell(m,1);
    if n >= 4
        conversionOp = keys(:,4);       % conversion operation: 'i', 'd', or 's' ('') - for each key to be read
    else
        conversionOp = cellstrings(m,1);
    end
    if n < 5
        defaultValues = cellstrings(m,1);
    else
        defaultValues = keys(:,5);
    end
    readsett = defaultValues;
    keysIn = keys(:,1:3);
    [secsExist,subsecsExist,keysExist,readValues,so,eo] = findkeysInString(myString,keysIn);    
    ind = find(keysExist);
    % For those keys that exist but have empty values, replace them with
    % the default values
    if ~isempty(ind)
        ind_empty = zeros(size(ind));
        for kk = 1:size(ind,1)
            ind_empty(kk) = isempty(readValues{ind(kk)});
        end
        ind(find(ind_empty)) = [];
        readsett(ind) = readValues(ind);
    end 
    % Now, go through all the keys and do the conversion if the conversion
    % char is given
    for ii=1:m
        if ~isempty(conversionOp{ii}) & ~strcmpi(conversionOp{ii},'s')
            if strcmpi(conversionOp{ii},'i') | strcmpi(conversionOp{ii},'d')
                if ~isnumeric(readsett{ii})
                    readsett{ii} = str2num(readsett{ii});
                end
                if strcmpi(conversionOp{ii},'i')
                    readsett{ii} = round(readsett{ii});
                end
                if isempty(readsett{ii})
                    result{ii} = [num2str(ii) '-th key ' keysIn{ii,3} 'or given defaultValue could not be converted using ''' conversionOp{ii} ''' conversion'];
                end
            else
                error(['Invalid conversion char given: ' conversionOp{ii}]);
            end
        end
    end
    varargout(1) = {readsett};
    varargout(2) = {result};
    return
    
%----------------------------
% WRITES key-value pairs to an existing or non-existing
% file (file can even be empty)
%----------------------------  
elseif (strcmpi(operation,'write'))
    if m < 1
        error('At least one key is needed when writing keys');
    end
    if ~exist(fileName)
        inifile(fileName,'new');
    end
    for ii=1:m  % go through ALL the keys and convert them to strings
        keys{ii,4} = n2s(keys{ii,4});
    end
    writekeys(fileName,keys,style);
    return
    
%----------------------------
% DELETES key-value pairs out
%----------------------------        
elseif (strcmpi(operation,'deletekeys'))
    deletekeys(fileName,keys);
    
    
    
else
    error('Unknown operation for INIFILE.');
end      




%--------------------------------------------------
%%%%%%%%%%%%% SUBFUNCTIONS SECTION %%%%%%%%%%%%%%%%
%--------------------------------------------------

% ------------------------------------
%> @brief parse ini string for keys given by 'keysIn'
%>
%> [secsExist,subSecsExist,keysExist,values,startOffsets,endOffsets] = findkeysInString(myString,keysIn)
%>
%> This function parses a string for keys as given by keysIn. keysIn is a cell
%> array of strings having 3 columns; section, subsection and key in each row.
%> section and/or subsection can be empty (root section or root subsection)
%> but the key can not be empty. The startOffsets and endOffsets are start and
%> end bytes that each key occuppies, respectively. If any of the keys doesn't exist,
%> startOffset and endOffset for this key are the same. A special case is
%> when the key that doesn't exist also corresponds to a non-existing
%> section and non-existing subsection. In such a case, the startOffset and
%> endOffset have values of -1.
%>
%> @param myString string to search
%> @param keysIn cell array of keys to find
%> 
%> @retval secsExist set if section exists and is not empty
%> @retval subSecsExists set if sub subsection exists and is not empty
%> @retval keyExist set if key that we are looking for exists
%> @retval keysLocated set if the key's position is located
%> @retval values read values of keys (strings)
%> @retval startOffsets start location of the keys
%> @retval endOffsets end location of the keys
%> 
function [secsExist,subSecsExist,keysExist,values,startOffsets,endOffsets] = findkeysInString(myString,keysIn)

nKeys = size(keysIn,1);         % number of keys
nKeysLocated = 0;               % number of keys located
secsExist = zeros(nKeys,1);     % if section exists (and is non-empty)
subSecsExist = zeros(nKeys,1);  % if subsection...
keysExist = zeros(nKeys,1);     % if key that we are looking for exists
keysLocated = keysExist;        % if the key's position (existing or non-existing) is LOCATED
values = cellstrings(nKeys,1);  % read values of keys (strings)
startOffsets = -ones(nKeys,1);  % start byte-position of the keys
endOffsets = -ones(nKeys,1);    % end byte-position of the keys

keyInd = find(strcmpi(keysIn(:,1),''));  % key indices having [] section (root section)

line = [];
lineN = 0;                      % line number
currSection = '';
currSubSection = '';

remstr = myString;

try
    %--- Searching for the keys - their values and start and end locations in bytes
    while 1
       
        %pos1 = ftell(fh);

        [line,remstr] = strtok(remstr,sprintf('\n'));
        pos1 = length(line)+2;

        if isempty(remstr)               % end of file, exit
            line = [];
            break
        end
        lineN = lineN + 1;
        [status,readValue,readKey] = processiniline(line);        
        if (status == 1)            % (new) section found
            % Keys that were found as belonging to any previous section
            % are now assumed as located (because another
            % section is found here which could even be a repeated one)
            keyInd = find( ~keysLocated & strcmpi(keysIn(:,1),currSection) );  
            if length(keyInd)
                keysLocated(keyInd) = 1;
                nKeysLocated = nKeysLocated + length(keyInd);
            end
            currSection = readValue;
            currSubSection = '';
            % Indices to non-located keys belonging to current section
            keyInd = find( ~keysLocated & strcmpi(keysIn(:,1),currSection) );  
            if ~isempty(keyInd)
                secsExist(keyInd) = 1;
            end
            pos2 = pos1;
            startOffsets(keyInd) = pos2+1;
            endOffsets(keyInd) = pos2+1;
        elseif (status == 2)        % (new) subsection found
            % Keys that were found as belonging to any PREVIOUS section
            % and/or subsection are now assumed as located (because another
            % subsection is found here which could even be a repeated one)
            keyInd = find( ~keysLocated & strcmpi(keysIn(:,1),currSection) & ~keysLocated & strcmpi(keysIn(:,2),currSubSection));
            if length(keyInd)
                keysLocated(keyInd) = 1;
                nKeysLocated = nKeysLocated + length(keyInd);
            end
            currSubSection = readValue;
            % Indices to non-located keys belonging to current section and subsection at the same time
            keyInd = find( ~keysLocated & strcmpi(keysIn(:,1),currSection) & ~keysLocated & strcmpi(keysIn(:,2),currSubSection));
            if ~isempty(keyInd)
                subSecsExist(keyInd) = 1;
            end
            pos2 = pos1;
            startOffsets(keyInd) = pos2+1;
            endOffsets(keyInd) = pos2+1;
        elseif (status == 3)        % key found
            if isempty(keyInd)
                continue            % no keys from 'keys' - from section-subsection par currently in
            end
            currKey = readValue;
            pos2 = pos1;       % the last-byte position of the read key  - the total sum of chars read so far
            for ii=1:length(keyInd)
               if strcmpi( keysIn(keyInd(ii),3),readKey ) & ~keysLocated(keyInd(ii))
                   keysExist(keyInd(ii)) = 1;
                   startOffsets(keyInd(ii)) = pos1+1;
                   endOffsets(keyInd(ii)) = pos2;
                   values{keyInd(ii)} = currKey;
                   keysLocated(keyInd(ii)) = 1;
                   nKeysLocated = nKeysLocated + 1;
               else
                   if ~keysLocated(keyInd(ii))
                       startOffsets(keyInd(ii)) = pos2+1;
                       endOffsets(keyInd(ii)) = pos2+1;
                   end
               end
            end
            if nKeysLocated >= nKeys  % if all the keys are located stop the searching
                break
            end
        else                          % general text found (even empty line(s))
            if (status == -1)
                error(['unknown string found at line ' num2str(lineN)]);
            end
        end       
    %--- End of searching
    end    
catch
    error(['Error parsing the file for keys: ' fileName ': ' lasterr]);
end


% ------------------------------------
%> @brief parse ini file for keys given by 'keysIn'
%>
%> [secsExist,subSecsExist,keysExist,values,startOffsets,endOffsets] = findkeysInString(fileName,keysIn)
%>
%> This function parses a file for keys as given by keysIn. keysIn is a cell
%> array of strings having 3 columns; section, subsection and key in each row.
%> section and/or subsection can be empty (root section or root subsection)
%> but the key can not be empty. The startOffsets and endOffsets are start and
%> end bytes that each key occuppies, respectively. If any of the keys doesn't exist,
%> startOffset and endOffset for this key are the same. A special case is
%> when the key that doesn't exist also corresponds to a non-existing
%> section and non-existing subsection. In such a case, the startOffset and
%> endOffset have values of -1.
%>
%> @param fileName file to search
%> @param keysIn cell array of keys to find
%> 
%> @retval secsExist set if section exists and is not empty
%> @retval subSecsExists set if sub subsection exists and is not empty
%> @retval keyExist set if key that we are looking for exists
%> @retval keysLocated set if the key's position is located
%> @retval values read values of keys (strings)
%> @retval startOffsets start location of the keys
%> @retval endOffsets end location of the keys
%> 
%>
function [secsExist,subSecsExist,keysExist,values,startOffsets,endOffsets] = findkeys(fileName,keysIn)

nKeys = size(keysIn,1);         % number of keys
nKeysLocated = 0;               % number of keys located
secsExist = zeros(nKeys,1);     % if section exists (and is non-empty)
subSecsExist = zeros(nKeys,1);  % if subsection...
keysExist = zeros(nKeys,1);     % if key that we are looking for exists
keysLocated = keysExist;        % if the key's position (existing or non-existing) is LOCATED
values = cellstrings(nKeys,1);  % read values of keys (strings)
startOffsets = -ones(nKeys,1);  % start byte-position of the keys
endOffsets = -ones(nKeys,1);    % end byte-position of the keys

keyInd = find(strcmpi(keysIn(:,1),''));  % key indices having [] section (root section)

line = [];
lineN = 0;                      % line number
currSection = '';
currSubSection = '';

fh = fopen(fileName,'r');
if fh == -1
    error(['File: ''' fileName ''' does not exist or can not be opened.']);
end

try
    %--- Searching for the keys - their values and start and end locations in bytes
    while 1
       
        pos1 = ftell(fh);
        line = fgetl(fh);
        if line == -1               % end of file, exit
            line = [];
            break
        end
        lineN = lineN + 1;
        [status,readValue,readKey] = processiniline(line);        
        if (status == 1)            % (new) section found
            % Keys that were found as belonging to any previous section
            % are now assumed as located (because another
            % section is found here which could even be a repeated one)
            keyInd = find( ~keysLocated & strcmpi(keysIn(:,1),currSection) );  
            if length(keyInd)
                keysLocated(keyInd) = 1;
                nKeysLocated = nKeysLocated + length(keyInd);
            end
            currSection = readValue;
            currSubSection = '';
            % Indices to non-located keys belonging to current section
            keyInd = find( ~keysLocated & strcmpi(keysIn(:,1),currSection) );  
            if ~isempty(keyInd)
                secsExist(keyInd) = 1;
            end
            pos2 = ftell(fh);
            startOffsets(keyInd) = pos2+1;
            endOffsets(keyInd) = pos2+1;
        elseif (status == 2)        % (new) subsection found
            % Keys that were found as belonging to any PREVIOUS section
            % and/or subsection are now assumed as located (because another
            % subsection is found here which could even be a repeated one)
            keyInd = find( ~keysLocated & strcmpi(keysIn(:,1),currSection) & ~keysLocated & strcmpi(keysIn(:,2),currSubSection));
            if length(keyInd)
                keysLocated(keyInd) = 1;
                nKeysLocated = nKeysLocated + length(keyInd);
            end
            currSubSection = readValue;
            % Indices to non-located keys belonging to current section and subsection at the same time
            keyInd = find( ~keysLocated & strcmpi(keysIn(:,1),currSection) & ~keysLocated & strcmpi(keysIn(:,2),currSubSection));
            if ~isempty(keyInd)
                subSecsExist(keyInd) = 1;
            end
            pos2 = ftell(fh);
            startOffsets(keyInd) = pos2+1;
            endOffsets(keyInd) = pos2+1;
        elseif (status == 3)        % key found
            if isempty(keyInd)
                continue            % no keys from 'keys' - from section-subsection par currently in
            end
            currKey = readValue;
            pos2 = ftell(fh);       % the last-byte position of the read key  - the total sum of chars read so far
            for ii=1:length(keyInd)
               if strcmpi( keysIn(keyInd(ii),3),readKey ) & ~keysLocated(keyInd(ii))
                   keysExist(keyInd(ii)) = 1;
                   startOffsets(keyInd(ii)) = pos1+1;
                   endOffsets(keyInd(ii)) = pos2;
                   values{keyInd(ii)} = currKey;
                   keysLocated(keyInd(ii)) = 1;
                   nKeysLocated = nKeysLocated + 1;
               else
                   if ~keysLocated(keyInd(ii))
                       startOffsets(keyInd(ii)) = pos2+1;
                       endOffsets(keyInd(ii)) = pos2+1;
                   end
               end
            end
            if nKeysLocated >= nKeys  % if all the keys are located stop the searching
                break
            end
        else                          % general text found (even empty line(s))
            if (status == -1)
                error(['unknown string found at line ' num2str(lineN)]);
            end
        end       
    %--- End of searching
    end    
    fclose(fh);
catch
    fclose(fh);
    error(['Error parsing the file for keys: ' fileName ': ' lasterr]);
end
%------------------------------------




%------------------------------------
function writekeys(fileName,keys,style)
% Writes keys to the section and subsection pair
% If any of the keys doesn't exist, a new key is added to
% the end of the section-subsection pair otherwise the key is updated (changed).
% Keys is a 4-column cell array of strings.

global NL_CHAR;

RETURN = sprintf('\r');
NEWLINE = sprintf('\n');

[m,n] = size(keys);
if n < 4
    error('Keys to be written are given in an invalid format.');
end

% Get keys position first using findkeys
keysIn = keys;    
[secsExist,subSecsExist,keysExist,readValues,so,eo] = findkeys(fileName,keys(:,1:3));

% Read the whole file's contents out
fh = fopen(fileName,'r');
if fh == -1
    error(['File: ''' fileName ''' does not exist or can not be opened.']);
end
try
    dataout = fread(fh,'char=>char')';
catch
    fclose(fh);
    error(lasterr);
end
fclose(fh);

%--- Rewriting the file -> writing the refined contents
fh = fopen(fileName,'w');
if fh == -1
    error(['File: ''' fileName ''' does not exist or can not be opened.']);
end
try
    tab1 = [];
    if strcmpi(style,'tabbed')
        tab1 = sprintf('\t');
    end
    % Proper sorting of keys is cruical at this point in order to avoid
    % inproper key-writing.
    
    % Find keys with -1 offsets - keys with non-existing section AND
    % subsection - keys that will be added to the end of the file   
    fs = length(dataout);       % file size in bytes
    nAddedKeys = 0;
    ind = find(so==-1);
    if ~isempty(ind)        
        so(ind) = (fs+10);      % make sure these keys will come to the end when sorting
        eo(ind) = (fs+10);
        nAddedKeys = length(ind);
    end
    
    % Sort keys according to start- and end-offsets
    [dummy,ind] = sort(so,1);
    so = so(ind);
    eo = eo(ind);
    keysIn = keysIn(ind,:);
    keysExist = keysExist(ind);
    secsExist = secsExist(ind);
    subSecsExist = subSecsExist(ind);
    readValues = readValues(ind);
    values = keysIn(:,4);   

    % Find keys with equal start offset (so) and additionally sort them
    % (locally). These are non-existing keys, including the ones whose
    % section and subsection will also be added.
    nKeys = size(so,1);
    fullInd = 1:nKeys;
    ii = 1;
    while ii < nKeys
        ind = find(so==so(ii));
        if ~isempty(ind) && length(ind) > 1
            n = length(ind);
            from = ind(1);
            to = ind(end);
            tmpKeys = keysIn( ind,: );
            [tmpKeys,ind2] = sortrows( lower(tmpKeys) );
            fullInd(from:to) = ind(ind2);
            ii = ii + n;
        else
            ii = ii + 1;
        end
    end
    
    % Final (re)sorting
    so = so(fullInd);
    eo = eo(fullInd);
    keysIn = keysIn(fullInd,:);
    keysExist = keysExist(fullInd);
    secsExist = secsExist(fullInd);
    subSecsExist = subSecsExist(fullInd);
    readValues = readValues(fullInd);
    values = keysIn(:,4);    
    
    % Refined data - datain
    datain = [];
    
    for ii=1:nKeys      % go through all the keys, existing and non-existing ones
        if ii==1
            from = 1;   % from byte-offset of original data (dataout)
        else
            from = eo(ii-1);
            if keysExist(ii-1)
                from = from + 1;
            end
        end
        to = min(so(ii)-1,fs);  % to byte-offset of original data (dataout)
        
        if ~isempty(dataout)
            datain = [datain dataout(from:to)];    % the lines before the key
        end
        
        if length(datain) & (~(datain(end)==RETURN | datain(end)==NEWLINE))
            datain = [datain, sprintf(NL_CHAR)];
        end

        tab = [];
        if ~keysExist(ii) 
            if ~secsExist(ii) && ~isempty(keysIn(ii,1))
                if ~isempty(keysIn{ii,1})
                    datain = [datain sprintf(['%s' NL_CHAR],['[' keysIn{ii,1} ']'])];
                end                
                % Key-indices with the same section as this, ii-th key (even empty sections are considered)
                ind = find( strcmpi( keysIn(:,1), keysIn(ii,1)) );
                % This section exists at all keys  corresponding to the same section from know on (even the empty ones)
                secsExist(ind) = 1;
            end
            if ~subSecsExist(ii) && ~isempty(keysIn(ii,2))
                if ~isempty( keysIn{ii,2})
                    if secsExist(ii); tab = tab1;  end;
                    datain = [datain sprintf(['%s' NL_CHAR],[tab '{' keysIn{ii,2} '}'])];
                end
                % Key-indices with the same section AND subsection as this, ii-th key
                % (even empty sections and subsections are considered)
                ind = find( strcmpi( keysIn(:,1), keysIn(ii,1)) & strcmpi( keysIn(:,2), keysIn(ii,2)) );
                % This subsection exists at all keys corresponding to the
                % same section and subsection from know on (even the empty ones)
                subSecsExist(ind) = 1;
            end
        end
        if secsExist(ii) & (~isempty(keysIn{ii,1})); tab = tab1;  end;
        if subSecsExist(ii) & (~isempty(keysIn{ii,2})); tab = [tab tab1];  end;
        datain = [datain sprintf(['%s' NL_CHAR],[tab keysIn{ii,3} ' = ' values{ii}])];
    end
    from = eo(ii);
    if keysExist(ii)
        from = from + 1;
    end
    to = length(dataout);
    if from < to
        datain = [datain dataout(from:to)];
    end
    fwrite(fh,datain,'char');
catch
    fclose(fh);
    error(['Error writing keys to file: ''' fileName ''' : ' lasterr]);
end
fclose(fh);
%------------------------------------



%------------------------------------
function deletekeys(fileName,keys)
% Deletes keys and their values out; keys must have at least 3 columns:
% section, subsection, and the key

[m,n] = size(keys);
if n < 3
    error('Keys to be deleted are given in an invalid format.');
end

% Get keys position first
keysIn = keys;    
[secsExist,subSecsExist,keysExist,readValues,so,eo] = findkeys(fileName,keys(:,1:3));

% Read the whole file's contents out
fh = fopen(fileName,'r');
if fh == -1
    error(['File: ''' fileName ''' does not exist or can not be opened.']);
end
try
    dataout = fread(fh,'char=>char')';
catch
    fclose(fh);
    rethrow(lasterror);
end
fclose(fh);

%--- Rewriting the file -> writing the refined content
fh = fopen(fileName,'w');
if fh == -1
    error(['File: ''' fileName ''' does not exist or can not be opened.']);
end
try
    ind = find(keysExist);
    nExistingKeys = length(ind);
    datain = dataout;
    
    if nExistingKeys
        % Filtering - retain only the existing keys...
        fs = length(dataout);       % file size in bytes
        so = so(ind);
        eo = eo(ind);
        keysIn = keysIn(ind,:);
        % ...and sorting
        [so,ind] = sort(so);
        eo = eo(ind);
        keysIn = keysIn(ind,:);
        
        % Refined data - datain
        datain = [];
        
        for ii=1:nExistingKeys  % go through all the existing keys
            if ii==1
                from = 1;   % from byte-offset of original data (dataout)
            else
                from = eo(ii-1)+1;
            end
            to = so(ii)-1;  % to byte-offset of original data (dataout)
            
            if ~isempty(dataout)
                datain = [datain dataout(from:to)];    % the lines before the key
            end       
        end
        from = eo(ii)+1;
        to = length(dataout);
        if from < to
            datain = [datain dataout(from:to)];
        end
    end
    
    fwrite(fh,datain,'char');
catch
    fclose(fh);
    error(['Error deleting keys from file: ''' fileName ''' : ' lasterr]);
end
fclose(fh);
%------------------------------------




%------------------------------------
function [keys,sections,subsections] = readallkeys(fileName)
% Reads all the keys out as well as the sections and subsections

keys = [];
sections = [];
subsections = [];
% Read the whole file's contents out
try
    dataout = textread(fileName,'%s','delimiter','\n');
catch
    error(['File: ''' fileName ''' does not exist or can not be opened.']);
end
nLines = size(dataout,1);

% Go through all the lines and construct the keys variable
keys = cell(nLines,4);
sections = cell(nLines,1);
subsections = cell(nLines,2);
keyN = 0;
secN = 0;
subsecN = 0;
secStr = '';
subsecStr = '';
for ii=1:nLines
    [status,value,key] = processiniline(dataout{ii});
    if status == 1
        secN = secN + 1;
        secStr = value;
        sections(secN) = {secStr};
    elseif status == 2
        subsecN = subsecN + 1;
        subsecStr = value;
        subsections(subsecN,:) = {secStr,subsecStr};
    elseif status == 3
        keyN = keyN + 1;
        keys(keyN,:) = {secStr,subsecStr,key,value};
    end
end
keys(keyN+1:end,:) = [];
sections(secN+1:end,:) = [];
subsections(subsecN+1:end,:) = [];
%------------------------------------

%------------------------------------
function [keys,sections,subsections] = readallkeysfromstring(mystring)
% Reads all the keys out as well as the sections and subsections

keys = [];
sections = [];
subsections = [];
% Read the whole file's contents out
try
    dataout = strread(mystring,'%s','delimiter','\n');
catch
    error(['File: ''' fileName ''' does not exist or can not be opened.']);
end
nLines = size(dataout,1);

% Go through all the lines and construct the keys variable
keys = cell(nLines,4);
sections = cell(nLines,1);
subsections = cell(nLines,2);
keyN = 0;
secN = 0;
subsecN = 0;
secStr = '';
subsecStr = '';
for ii=1:nLines
    [status,value,key] = processiniline(dataout{ii});
    if status == 1
        secN = secN + 1;
        secStr = value;
        sections(secN) = {secStr};
    elseif status == 2
        subsecN = subsecN + 1;
        subsecStr = value;
        subsections(subsecN,:) = {secStr,subsecStr};
    elseif status == 3
        keyN = keyN + 1;
        keys(keyN,:) = {secStr,subsecStr,key,value};
    end
end
keys(keyN+1:end,:) = [];
sections(secN+1:end,:) = [];
subsections(subsecN+1:end,:) = [];
%------------------------------------



%------------------------------------
function [status,value,key] = processiniline(line)
% Processes a line read from the ini file and
% returns the following values:
%   - status:  -1   => unknown string found
%               0   => empty line found
%               1   => section found
%               2   => subsection found
%               3   => key-value pair found
%               4   => comment line found (starting with ;)
%   - value:    value-string of a key, section, subsection, comment, or unknown string
%   - key:      key as string

status = 0;
value = [];
key = [];

% WID - strip comments from line.  Results in removing all comments from ini 
%       file so only use when comments are not required.
[a,b] = regexp(line, '^(.*?);.*', 'tokens');
if ~isempty(b), line = a{1}{1}; end

line = strim(line);                         % removes any leading and trailing spaces

if isempty(line)                            % empty line
    return
end

if strcmpi(line(1),';')                     % comment found
    status = 4;
    value = line(2:end);
elseif (line(1) == '[') & (line(end) == ']') & (length(line) >= 3)  % section found
    value = lower(line(2:end-1));
    status = 1;
elseif (line(1) == '{') &...                % subsection found
       (line(end) == '}') & (length(line) >= 3)
    value = lower(line(2:end-1));
    status = 2;
else                                        % either key-value pair or unknown string
    pos = findstr(line,'=');
    if ~isempty(pos)                        % key-value pair found
        status = 3;
        key = lower(line(1:pos-1));
        value = line(pos+1:end);
        key = strim(key);                   % removes any leading and trailing spaces
        value = strim(value);               % removes any leading and trailing spaces
        if isempty(key)                     % empty keys are not allowed
            status = 0;
            key = [];
            value = [];
        end
    else                                    % unknown string found
        status = -1;
        value = line;
    end
end


%------------------------------------
function outstr = strim(str)
% Removes leading and trailing spaces (spaces, tabs, endlines,...)
% from the str string.
if isnumeric(str);
    outstr = str;
    return
end
ind = find( ~isspace(str) );        % indices of the non-space characters in the str    
if isempty(ind)
    outstr = [];        
else
    outstr = str( ind(1):ind(end) );
end



%------------------------------------
function cs = cellstrings(m,n)
% Creates a m x n cell array of empty strings - ''
cs = cell(m,n);
cs(:) = {''};


%------------------------------------
function y = n2s(x)
% Converts numeric matrix to string representation.
% Example: x given as [1 2;3 4] returns y = '1,2;3;4'
if ischar(x) | isempty(x)
    y = x;
    return
end
[m,n] = size(x);
y = [num2str(x(1,:),'%15.6g')];
for ii=2:m
    y = [y ';' num2str(x(ii,:),'%15.6g')];
end
