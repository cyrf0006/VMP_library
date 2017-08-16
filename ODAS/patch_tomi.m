%% patch_tomi
% Legacy function, replace with patch_odas. This function call is redirected to
% patch_odas.
%%
% <latex>\index{Type A!patch\_tomi}</latex>
%
%%% Syntax
%   [bad_record, fix_manually] = patch_tomi( file );
%
%%% Description
%
% Calls to patch_tomi can be directly replaced with patch_odas.  The calling 
% structure is the exact same so all that is required is for you to change the
% function name.
%

function  [bad_record, fix_manually] = patch_tomi(file)
warning('"patch_tomi" has been depreciated.  Please use "patch_odas" instead.');
[bad_record, fix_manually] = patch_odas(file);
