%% setupstr
% Read attributes from an ODAS configuration string.
%%
% <latex>\index{Type A!setupstr}</latex>
%
%%% Syntax
%
%   varargout = setupstr( varargin )
%
% * [varargin] see Usage
% * [varargout] see Usage
%
%%% Description
%
% This function will extract requested attributes from an input string. The
% string should be an ODAS configuration string - an INI-like string used to 
% store values used by the ODAS library.  This string can be obtained by reading
% a configuration file (.cfg) or by extracting the relevant parts of a ODAS v6
% or higher data file.
%
%%% Usage
%
% The setupstr function is called with up to four arguments.  Depending on the
% number and type of arguments used, different results are returned.
%
% The first argument is the configuration string, either extracted from a data
% file or loaded from a configuration file.  When this is the only argument, the
% setupstr function returns an indexed structure containing the contents of the
% configuration string.  This indexed structure can be used in place of the 
% configuration string for subsequent function calls thereby reducing the time
% required to make those calls.
%
%   - obj = SETUPSTR('config_string')
%      Parse the string 'config_string' and return 'obj'.  The returned value 
%      is a data structure containing values within the config_string.  The 
%      returned 'obj' is used in subsequent calls to the function.
%
%   - S = SETUPSTR( obj, 'section' )
%      Find sections within 'obj' that match the value 'section'.  See Notes
%      for additional information regarding search queries.  The returned value
%      'S' is a cell array containing the matching section identifiers in string
%      format.
%
%   - V = SETUPSTR( obj, 'section', 'parameter' )
%      Find parameter values within 'obj' that are located in 'section' and 
%      have the parameter name 'parameter'.  See Notes for additional 
%      information regarding search queries.  The returned value 'V' is a cell 
%      array containing the matching string values.
%
%   - [S, P, V] = SETUPSTR( obj, 'section', 'parameter', 'value')
%      Find values within 'obj' that are located in 'section', have the 
%      parameter name 'parameter', and the parameter value 'value'.  See Notes
%      for additional information regarding search queries.  The returned values
%      are a tuple of cell arrays containing the matching sections, parameter 
%      names, and parameter values.
%
% * [config_string] Configuration string.
% * [obj]     Either the configuration string or a previously returned 
%             structure containing indexed values from a configuration string.
% * [section] Search query for a matching section identifier.
% * [parameter] Search query for a matching parameter name.
% * [value]   Search query for a matching parameter value.
% * []
% * [obj] Structure containing indexed values from a configuration file. Passed
%         to subsequent calls to this function to speed the search.
% * [S]   Matching sections as an array of cells.
% * [P]   Matching parameter names as an array of cells.
% * [V]   Matching parameter values as an array of cells.
%
%%% Notes
%
% All search queries are interpreted as regular expressions.  If empty, the 
% query will match anything by using the regular expression '.*'.  Each search 
% query is prepended with '^' and appended with '$' to ensures exact matches.
%
% Parameters that precede the first section are said to be in the 'root' 
% section.  To access these parameters, use 'root' as the section identifier.
%
%%% Examples
%
%    >> cfg = setupstr( setupfilestr );
%    >> value = setupstr( cfg, 'P', 'coef0' )
%
% From the variable 'setupfilestr', which contains a configuration file, query
% the value of 'coef0' for the pressure channel.
%
%    >> sections = setupstr( setupfilestr, '' )
%
% Find all sections within a configuration file.
%
%    >> cfg = setupstr( setupfilestr );
%    >> [S,P,V] = setupstr( cfg, '', 'id.*', '' )
%
% Find all sections that have an 'id' value.  Allow for 'id_even' and 'id_odd'.
%
% Note: the above example will match any key starting with 'id'.  For a more
% precise match, 'id.*' could be changed to 'id(_(even|odd))?'.

% Version History:
% 
% * 2012-04-30 (WID) initial function
% * 2012-10-24 (WID) Updated documentation
% * 2012-11-22 (WID) added support for arg1 being a configuration string or a
%                    structure.
% * 2012-11-22 (WID) updated documentation to match the manual
% * 2012-02-26 (WID) updated documentation - incorrect pressure example

function varargout = setupstr(arg1, arg2, arg3, arg4)

if nargin == 0,
  display_help();
  varargout{1} = [];      % default return value prevents error during help.
end


if nargin == 1
  varargout{1} = parse_setup_file (arg1);
end


% Assume a section expression was provided, return matching section.
if nargin == 2
  if isa( arg1, 'char' )
	arg1 = parse_setup_file( arg1 );
  end
  out = {};
  for i = section_indexes(arg1, arg2)
    out{end+1} = arg1(i).k{1};
  end
  varargout{1} = out;
end


% Assume a section and key were provided, return all matching values.
if nargin == 3
  if isa( arg1, 'char' )
	arg1 = parse_setup_file( arg1 );
  end
  out = {};
  for i = section_indexes(arg1, arg2)
    for ii = key_indexes(arg1, i, arg3)
      out{end+1} = arg1(i).v(ii).v;
    end
  end
  varargout{1} = out;
end


% Assume a section, key, and value were provided - return all matching truples.
if nargin == 4
  if isa( arg1, 'char' )
	arg1 = parse_setup_file( arg1 );
  end
  if isempty(arg4), value = '.*'; else value = ['^' arg4 '$']; end
  sections = [];
  keys = [];
  values = [];
  for i = section_indexes(arg1, arg2)
    for ii = key_indexes(arg1, i, arg3)
      if ~isempty(regexpi(char(arg1(i).v(ii).v), value, 'match'))
        sections{end+1} = arg1(i).k{1};
        keys{end+1} = arg1(i).v(ii).k;
        values{end+1} = arg1(i).v(ii).v;
      end
    end    
  end
  varargout{1} = sections;
  varargout{2} = keys;
  varargout{3} = values;
end

end



function display_help ()
disp(' ');
disp('SETUPSTR:  parse / query ODAS setup strings');
disp(' ');
disp('  OBJ = setupstr (setup_string)');
disp('    Parse setup_string into an array of struct elements.  The returned');
disp('    array should be passed into subsequent function calls.');
disp(' ');
disp('  S = setupstr (struct, section)');
disp('    Search for matching sections in the struct array.  Standard regex');
disp('    wildcards are accepted.  A cell array of matching sections is');
disp('    returned.');
disp(' ');
disp('  V = setupstr (struct, section, parameter)');
disp('    Search for matching section/parameter pairs in the struct array and');
disp('    return a cell array of associated values.  Standard regex wildcards');
disp('    are accepted.');
disp(' ');
disp('  [S, K, V] = setupstr (struct, section, parameter, value)');
disp('    Search for section, parameter, and value entries.  Standard regex');
disp('    wildcards are accepted.  Cell arrays with the matching sections,');
disp('    parameter, and values are returned.');
end



function out = parse_setup_file (setupFile)

out = [];
section = 'root';
s_index = 1;
index = 1;
dest = 1;

lines = strtrim(regexp(setupFile, '\n', 'split'));

for i=1:length(lines)
  
  % remove comments
  [a,b] = regexp(lines{i}, '^(.*?);.*', 'tokens');
  if ~isempty(b), 
	  if ~isempty(a{1}), lines{i} = a{1}{1}; end
  end
  
  % find sections
  [a,b] = regexp(lines{i}, '^\s*\[\s*(.+)\s*\]', 'tokens');
  if ~isempty(b)
    section = lower(a{1}{1});
    index = 1;               % reset the index..
    dest = 1;
    s_index = s_index + 1;
    continue;
  end
  
  % find assignments
  expr = '^\s*(.+?)\s*=\s*(.+?)\s*$';
  [a,b] = regexp(lines{i}, expr, 'tokens');
  if ~isempty(b)
    
    % check if the key is equal to "name" - if so create a new section 
    % using the value of "name" if it doesn't match the current section.
    % Use 'dest' to ensure the first element in the array is default name.
    if strcmpi(a{1}{1}, 'name'),
      out(s_index).k(1) = lower(a{1}(2));
      dest = 2;
    end
    
    % we found an entry - time to add it.
    out(s_index).k(dest) = {section};
    out(s_index).v(index).k = lower(a{1}{1});
    
    % When we add the entry value, we want to add the evaluated version of the
    % entry if available.  For example, "4 + 7" should result in the string "11"
    % being added as the value.
    try
      temp = sprintf('%d',eval( a{1}{2} ));
    catch
      temp = a{1}{2};
    end
    out(s_index).v(index).v = temp;
    
    index = index + 1;
    continue;
  end
end
end


function out = section_indexes (input, expression)
% return list of index values of struct elements where the section matches
% the provided expression.

if isempty(expression)
  sect = '.*'; 
else
  if iscell(expression), expression = expression{1}; end
  sect = ['^' expression '$'];
end

out = [];
for i=1:length(input)
  for ii=1:length(input(i).k)
    result = regexpi(input(i).k{ii}, sect, 'match');
    if ~isempty(result), out = [out i]; break; end
  end
end

end


function out = key_indexes (input, section, expression)
% return list of index values at the given section where the key matches
% the provided expression.

if isempty(expression)
  key = '.*'; 
else
  if iscell(expression), expression = expression{1}; end
  key = ['^' expression '$'];
end

out = [];
for i=1:length(input(section).v)
  result = regexpi(input(section).v(i).k, key, 'match');
  if ~isempty(result), out = [out i]; end
end

end


