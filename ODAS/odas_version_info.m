%% odas_version_info - EXCLUDED
% Return the version of the current ODAS library.
%%
% <latex>\index{Functions!odas\_version\_info}</latex>
%
%%% Syntax
%   version = odas\_version\_info( )
%
% * []
% * [version] Numeric value of the version.
%
%%% Description
% Function called by other functions to determin the current version of
% the ODAS Libarary.  Used when checking for outdated data files.

% Version History
%
% * 2015-06-03 (WID) initial version 4.0

function version = odas_version_info( )

version = 4.02;

end

