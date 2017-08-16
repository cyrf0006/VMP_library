%% fix_underscore
% Escape underscore in supplied string for LaTeX
%%
% <latex>\index{Type B!fix\_underscore}</latex>
%
%%% Syntax
%   stringOut = fix_underscore( stringIn )
%
% * [stringIn] Input string.
% * []
% * [stringOut] Copy of the input string with all underscore characters 
%         prepended with a forward slash '\'.
%
%%% Description
% Escape the underscore character '_' into '\_'  within a string so that Matlab 
% functions do not interpret it as a subscript command.  This function
% allows file names to contain an underscore within plot text.

% Version History
%
% * 2000-07-01 (RGL)
% * 2011-09-01 (AWS) added documentation tags for matlab publishing
% * 2012-05-14 (WID) documentation update
% * 2012-11-05 (WID) documentation update

function string_out=fix_underscore(string_in)

l = length (string_in); 
string_out = string_in;
n = find(string_out == '_');
if ~isempty(n)
	for k = 1:length(n)
		string_out = [string_out(1:n(k)-1) '\_' string_out(n(k)+1:l)];
		n = find(string_out == '_');
		l = length(string_out); % get new length
	end
end % got all of the under_scores
