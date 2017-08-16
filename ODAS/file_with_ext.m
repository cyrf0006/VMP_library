%% file_with_ext
% Find file with default extension
%%
% <latex>\index{Type A!file\_with\_ext}</latex>
%
%%% Syntax
%   [P,N,E,fullName] = file_with_ext( file, extensions, errorStr )
%
% * [file] Name of file with or without extension.
% * [extensions] Cell array of acceptable extensions. ex, {'.p' '.P' ''}
% * [errorStr] Error string used if error is to be triggered.  (optional)
% * []
% * [P] Full path to the requested file.
% * [N] Name of requested file.
% * [E] Extension of requested file.
% * [fullName] Full path and file name with path seperators. (optional)
%
%%% Description
% Find file with a default extension.  Convenience function that allows users 
% to input file names without an extension.
%
% This function searches for a file with the name 'file' and an extension from
% the array 'extensions'.  When found, it returns the file components [P,N,E]
% along with the fully qualified name, 'fullName'.  If not found the function 
% triggers the error 'errorStr'.  Should 'errorStr' be empty then no error is 
% triggered and empty values are returned for [P,N,E,fullName].
%
% The return value 'fullName' is provided to simplify addressing the resulting
% file.  Alternatively, this value can be constructed from the [P,N,E] values.
%
% The 'extensions' array lists the acceptable file extensions for the requested
% file.  Each extension should consist of a string value prepended with a 
% period (.).  To allow for the extension to be declared within the 'file' 
% variable, prepend the array with an empty string.  See the example for 
% details.
%
% Be sure to include both possible case values within 'extensions' - for 
% example, {'.p', '.P'}.  Some operating systems are case sensitive while others
% are not.
%
%%% Examples
%
%    >> file_err = 'No valid input file was found..';
%    >> [P,N,E] = file_with_ext( file, {'' '.p' '.P'}, file_err );
%
% Typical usage where 'file' is a string provided by the user.
%
%    >> [P,N,E] = file_with_ext( file, {'' '.p' '.P'} );
%    >> if isempty(N), %processing error...;  end
%
% Errors are optional.  If 'errorStr' is not provided, this function will not 
% trigger an error if the file is not found.
%
%    >> [P,N,E,fName] = file_with_ext( fName, 
%                                     {'' '.p' '.P'},
%                                     ['Unable to find file: ' fName] );
%    >> dos(['mv ' fName ' /dev/null']);
%
% Record the full matching file identifier in 'fName'.  This example 
% demonstrates how the function can be utilized in a single call.  The 
% subsequent dos command will effectively delete the file on a UNIX system.  
% Using the full name and path ensures the wrong file does not get deleted.

% Version History
%
% * 2012-04-30 (WID) initial
% * 2012-05-14 (WID) updated documentation
% * 2012-07-19 (WID) fixed iteration over cell array + support path inputs
% * 2012-08-14 (WID) Documentation update + some simplified code
% * 2012-09-09 (WID) Bug fix for testing if files are present.
% * 2012-11-05 (WID) Documentation update.
  
function [P,N,E,fullName] = file_with_ext(file, extensions, error_str)

if nargin < 2
  error('Invalid number of arguments.');
end

if nargin == 2
  error_str = [];
end

% Open the 2-byte binary data file for reading
[P,N,E]=fileparts(file);

% Test for the following extensions.
for ext = extensions
  ext = char(ext);
  fullName = fullfile(P, [N E ext]);
  
  if ~isempty(dir(fullName)),
	  [P N E] = fileparts(fullName);
	  return
  end
end

% No extensions were found.  Return an error if requested.
P = ''; N = ''; E = ''; fullName = '';
if ~isempty(error_str), error(error_str); end

end
