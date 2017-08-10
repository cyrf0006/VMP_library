function Trbr = getrbr(filename, time_vec)


% I want to make it that way... after QO! 
% function [trbr, Trbr, Prbr] = getrbr(filename, varargin)
% function [trbr Trbr] = getrbr(filename)
%
% This functions is called in T_MN080.m or T_MRIKI.m to get data
% from an RBR thermistor.
    
t0 = time_vec(1);
tf = time_vec(end);
dt = time_vec(2) - time_vec(1);
    
    
load(filename)
tstart = RBR.starttime;
t1 = datenum([tstart(7:10) '-' tstart(4:5) '-' tstart(1:2) tstart(11:end)]);
dt_rbr = RBR.sampleperiod; % in sec.
dt_rbr = dt_rbr/60/60/24;
t2 = t1+size(RBR.data, 1)*dt_rbr;
trbr = [t1:dt_rbr:t2];
T_rbr = RBR.data(:,1)';

% to required resolution
for i = 1:length(time_vec)
    I = find(trbr >= time_vec(i)-dt/2 & trbr <= time_vec(i)+dt/2);
    Trbr(i) = nanmean(T_rbr(I));
end

