%% get_latest_file
% Get the name of the latest ODAS binary data file in current folder.
%%
% <latex>\index{Type A!get\_latest\_file}</latex>
%
%%% Syntax
%   testString = get_latest_file()
%
% * [testString] Name of the latest ODAS binary data file in the current 
%          directory. Empty if no ODAS binary data files exist.
%
%%% Description
% Retrieve the name of the latest ODAS binary data file in the local
% directory. It is useful for near real-time data processing and plotting
% because such operations are frequently conducted on the latest data file
% recorded in the local directory.  
%

% Version History
% * 2004-06-06 (RGL) initial
% * 2011-09-01 (AWS) added documentation tags for matlab publishing
% * 2012-04-23 (WID) Updated to be case insensitive.  Simplified loops.
% * 2012-07-19 (WID) simplified loop

function     latest_file = get_latest_file()

% Set initial value to zero so an empty string is returned if no files exist.
newest.datenum = 0;
newest.name = [];

for file = [dir('*.p')' dir('*.P')']
  if file.datenum > newest.datenum, newest = file; end
end

latest_file = newest.name;
