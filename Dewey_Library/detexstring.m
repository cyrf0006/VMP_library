function [escapedStr] = detexstring(inStr)
%
% detexstring.m--Given a string containing carats and underscores (that
% would otherwise be interpreted as Tex commands for super and subscripts),
% detexstring.m returns a string with these characters escaped, suitable
% for use as input to title.m, xlabel.m, etc.
% 
% The alternative to this is to set the text object's "Interpreter" to
% "none", but this will prevent the use of Greek and other special
% characters.  
%
% Syntax: escapedStr = detexstring(inStr)
%
% e.g.,   figure; filename = 'dummy_file_name.mat';
%         title(detexstring(['\Delta from file ' filename ' (underscores escaped)']));
%         xlabel(['\Delta from file ' filename ' (underscores NOT escaped)']);

% Developed in Matlab 6.5.0.180913a (R13) on SUN OS 5.8.
% Kevin Bartlett(kpb@hawaii.edu), 2003/11/17, 14:11
%------------------------------------------------------------------------------

escapedStr = strrep(inStr,'_','\_');
escapedStr = strrep(escapedStr,'^','\^');

