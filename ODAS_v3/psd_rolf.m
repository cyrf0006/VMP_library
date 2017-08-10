%% psd_rolf
% Calculate auto-spectrum. This is a legacy function, please use csd_odas.
%%
% <latex>\index{Depreciated!psd\_rolf}</latex>
%
%%% Syntax
%   [Px, Fs] = psd_rolf( x, n_fft, rate )
%
%%% Description
%
% Calls to psd_rolf can be replaced with csd_odas.  The calling structure
% for csd_odas differs from psd_rolf so some minor changes to your code will be 
% required.
%
% The three arguments to psd_rolf; x, n_fft, and rate; can all be passed
% into csd_odas but their location changes.  Assuming the above three arguments
% were used to call psd_rolf, one would call csd_odas as follows;
%
%    [Px, Fs] = csd_odas( x, x, n_fft, rate, [], n_fft/2, 'linear');


function [Px,Fs]=psd_rolf(x,n_fft,rate)
warning('The function "psd_rolf" has been depreciated.  Please use "csd_odas".');
[Px, Fs] = csd_odas(x,x,n_fft,rate,[],n_fft/2,'linear');
