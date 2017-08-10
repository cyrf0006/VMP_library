%% save_rolf
% Legacy function replaced by "save_odas"
%%
% <latex>\index{Depreciated!save\_rolf}</latex>
%
%%% Description
% This is a legacy function, please change your functions/scripts to use the
% function 'save_odas' instead. This function calls 'save_odas' using:
%
%    >> success_flag = save_odas( file_name, vector_name, vector )
%

function success_flag = save_rolf(file_name, vector_name, vector)
warning('The function "save_rolf" has been depreciated.  Please use "save_odas".');
success_flag = save_odas(file_name, vector_name, vector);
