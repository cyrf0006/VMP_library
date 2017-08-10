function [y, No] = Dspike(x, N, tresh)
    
% function I=Dspike(x, N, tresh)
% 
% Use to despike any signal. Fred's own version of despike.m 
% 
% - x: Initial signal to Despike
% - N: Smoothing window for running mean
% - tresh: Treshold for the ratio between x and smoothed signal.
%
% - y: Corrected signal where spikes are replace by value
% - No: Number of corrected points  
% interpolated between good point.    
%
% Author: Frederic Cyr, 2011/02/07
%
% ---------------------------------------------------------- %
    
    
    
I = find(~isnan(x)==1); % only consider non-NaNs
init = x(I); % initial signal, no Nans

smooth_sig = runmean(init,N);
ind = 1:length(smooth_sig);

J = find(init./smooth_sig > tresh);

init(J)=[];
IND = ind;
ind(J) = [];

corrected = interp1(ind, init, IND);

y = x;

y(I) = corrected;

if isempty(J)==1
    No=0;
else
    No = length(J);
end


    