%% csd_rolf
% Legacy function replaced by csd_odas
%%
% <latex>\index{Depreciated!csd\_rolf}</latex>
%
%%% Description
% This is a legacy function.  Please replace with the function 'csd_odas'.
%
% This function calls 'csd_odas' using:
%
%    [Cxy, Fs] = csd_odas( x, y, nFFT, rate, [], nFFT/2, 'linear' )
%

function [Cxy, F] = csd_rolf(x, y, n_fft, rate)
warning('The function "csd_rolf" has been depreciated.  Please use "csd_odas".');
[Cxy, F] = csd_odas(x,y,n_fft,rate,[],n_fft/2,'linear');
