%% csd_odas - EXCLUDED
% Legacy function replaced with csd_matrix_odas
%%
% <latex>\index{Depreciated!csd\_odas}</latex>
%
%%% Description
% Depreciated function replaced with csd_matrix_odas. The csd_matrix_odas
% function generates identicle results but operates more quickly due to a
% more efficient algorithm.  Please see the csd_matrix_function for
% information on the calling parameters and return values.

% *Version History:*
%
% * 2009-04-01 (RGL) Original function
% * 2011-08-23 (RGL) Modified to accept the case of y = [ ] for implicit
%   auto-spectrum. Forced Cxx and Cyy to [], for auto-spectrum
% * 2011-09-01 (AWS) added documentation tags for matlab publishing
% * 2012-07-04 (RGL) added parabolic and cubic detrending feature
% * 2015-10-22 (WID) Depreciated - now calls csd_matrix_odas.m.

function [Cxy, F, Cxx, Cyy] = csd_odas(x, y, n_fft, rate, Window, over_lap, msg)

[Cxy,F,Cxx,Cyy] = csd_matrix_odas(x,y,n_fft,rate,Window,over_lap,msg);

