%% despike
% Remove short-duration spikes from a signal.
%%
% <latex>\index{Functions!despike}</latex>
%
%%% Syntax
%   [y, spike, pass_count, ratio] = despike( dv, thresh, smooth, Fs, N, ...)
%
% * [dv] Signal to be despiked.
% * [thresh] Threshold value for the ratio of the instantaneous rectified
%       signal to its smoothed version. A value of 8 is a good starting
%       point for a VMP. 
% * [smooth] The cut-off frequency of the first-order Butterworth filter
%       that is used to smooth the rectified input signal. The time scale
%       of smoothing is approximately 1/(2*smooth). A value of 0.5 is a
%       good starting point for a VMP. 
% * [Fs] Sampling rate (Hz). 
% * [N] Spike removal scale. A total of 1.5*N data points are removed. N/2
%        points are removed before a spike, and N points are removed after
%        a spike. The replaced data is an average of the adjacent
%        neighbouthood. Averaging uses ~Fs/(4*smooth) points from each side
%        of a spike. 
% * ['-debug'] Optional string that provides a plot of the data after every pass.
% * ['-single_pass'] Optional string to force the function to make only one
%       pass.
% * []
% * [y] Despiked signal.
% * [spike] Indices to the spikes.
% * [pass_count] The number of iterations used to find and remove the spikes.
% * [ratio] The fraction of the data replaced with a local mean.
%
%%% Description
% This function removes spikes in shear probe and micro-conductivity
% signals, such as those generated by a collision with plankton and other
% detritus. It identifies spikes by comparing the instantaneous rectified
% signal against its local psuedo standard deviation. The instantaneous
% rectified signal is obtained by:
%
% # high-pass filtering the input signal, $\texttt{dv}$, at 0.5 Hz with a
% zero-phase, single-pole, Butterworth filter, and
% # rectifying this signal by taking its absolute value.
%
% The pseudo standard deviation is calculated by smoothing the rectified
% signal with a low-pass, single-pole, zero-phase, Butterworth filter with
% a cut-off frequency of $\texttt{smooth}\ \si{\hertz}$.
%
% The function applies itself iteratively. After identifying and replacing
% spikes, it tries again to identify and replace spikes, until no more
% spikes are detected. The maximun number of iterations is capped at 10.
%
% You can restrict despiking to a single pass, which duplicates the
% behaviour of the function in previous versions, using the
% $\texttt{'-single\_pass'}$ string option. You can visualize the iterative
% operation of the function using the $\texttt{'-debug'}$ string option.
%%% Examples
%
%    >> [y,spike,pass_count,ratio] = despike(sh1,8,0.5,fs_fast,round(0.04*fs_fast))
%
% @image @images/despike @Despike function example @Despike function
% applied with despike(sh1, 7, 0.1, 512, 30). The smoothing scale, the
% inverse of 0.1 Hz, is long because this profiler was moving at only
% $\SI{0.17}{\m\per\s}$. The units for shear should be $\si{\per\s}$ and
% not $\si{\per\square\s}$. 
%
% @image @images/despike2 @Despike example with zoomed in plot. @The same plot as the
% previous figure but with a close-up view.
%
% The parameters for this function should be chosen with care and
% thoroughly tested using the $\texttt{'-debug'}$ option.
%
% The number of points
% around a spike that are removed is best expressed in terms of a duration
% because the ringing of the shear probe after a collision is not speed
% dependent. For example, $\texttt{round(0.04*fs\_fast)}$ will remove 
% $\SI{40}{\milli\s}$ of data after a spike and half that amount before a spike.
% 
% The parameter $\texttt{smooth}$ is based on the notion that turbulence
% occurs over a path-length that exceeds some minimum, perhaps
% $\SI{0.5}{\m}$. Therefore, the smoothing time-scale and the cut-off
% frequency, $\texttt{smooth}$, of the filter are speed dependent. The
% value of $\texttt{smooth}$ is proportional to speed. 
%
% The threshold should
% be chosen so that an extrama is indeed anomalous with respect to a
% neighbourhood defined by the parameter $\texttt{smooth}$. A value of 8
% usually gives good result.
%
% 

% *Version History:*
%
% * 1995-07-27 (FW) initial
% * 1995-08-22 (FW) revised
% * 1995-08-28 (FW) revised
% * 1997-07-14 (FW) revised
% * 1998-08-31 (RGL) revised
% * 2000-05-03 (RGL) added unipulse feature
% * 2000-06-07 (RGL) cured problem with spikes near ends. Improved memory efficiency
% * 2000-06-21 (FW) finally cured problem with spikes near ends.
% * 2002-08-29 (IG) modified spike to correspond to correct location of spikes in y
% * 2009-03-25 (RGL) Set default rate to 512 Hz.
% * 2011-09-01 (AWS) added documentation tags for matlab publishing
% * 2012-05-07 (WID) update to documentation + removed outdated matlab
%                    functions
% * 2012-11-05 (WID) update to documentation
% * 2013-02-26 (WID) merge in changes from Rolf wrt input parameter parsing
% * 2013-10-28 (WID) replacement algorithm replaced.  Accept N as duration
%                    in seconds.
% * 2013-12-12 (WID) added option to return the removed points.  Other
%                    algorithm modifications.
% * 2013-12-23 (WID) small segments of valid data no longer removed in R
% * 2015-10-28 (RGL) Documentation corrections.
% * 2015-11-05 (RGL) Removed input options. Now requires all input
%                    parametrrs. Added a debug optioon to plot graphs of
%                    the rectified, with and without smoothing, and the
%                    results of debugging. Despiking is now down
%                    iteratively untill all spikes are removed. Added
%                    option for a single pass.
% * 2015-11-06 (WID) Reworked code to make it more managable and slightly
%                    faster. Results unchanged.
% * 2015-11-09 (RGL) Added the fraction of removed data to the return
%                    values.
% * 2015-11-18 (RGL) Documenttion corrections.

function [y, spike, pass_count, despike_fraction] = despike (dv, thresh, smooth, Fs, N, varargin)

% check for input parameters
if (nargin < 5), error('Too few input arguments to despike'), end

max_pass_count = 10;
generate_plot = false;
for arg = varargin
    if strcmpi(char(arg), '-debug')
        generate_plot = true;
    end
    if strcmpi(char(arg), '-single_pass')
        max_pass_count = 1;
    end
end

% Iterations of the despike function start here.
spike = [];
pass_count = 0;
y = dv;
despike_fraction = 0;

while pass_count < max_pass_count    
    [y,spike_pass,dv_LP,dv_HP] = single_despike(y, thresh, smooth, Fs, N, varargin);
    
    if isempty(spike_pass),
        despike_fraction = length(find(y ~= dv)) / length(y);
        return;
    end
    
    spike = union( spike', spike_pass' );
    pass_count = pass_count + 1;

    if generate_plot, make_figure(); end
end

    function [dv,spike,dv_LP,dv_HP] = single_despike(dv, thresh, smooth, Fs, N, varargin)
        
        N = N / 2;
        
        % Force into a column vector if required
        dv = dv(:);
        
        % Zero padding alleviates problems when spikes are at the
        % beginning or the end of the vector, and avoids problems
        % caused by filter transients of the smoothing filter
        % and rectifying filters. (Fab)
        len = length(dv);
        padLen = min(length(dv), 2*floor(Fs/smooth));    % padLen must not exceed the length of the input vector
        range = (1+padLen):(len+padLen); % these two values mark the start and end of the input
        % vector in the padded vector (used later)
        
        dv = [flipud(dv(1:padLen)); dv; flipud(dv(len-padLen+1:len))];
        
        % high pass filter & rectify
        [b,a] = butter(1, 0.5/(Fs/2), 'high');
        dv_HP = abs(filtfilt(b,a, dv));
        
        % smoothing coefficients
        [b_LP,a_LP] = butter(1, smooth/(Fs/2));
        
        % find spikes
        s = warning; % cache warning state
        warning off  % avoid divide-by-zero warnings (fab)
        dv_LP = filtfilt(b_LP,a_LP,dv_HP);
        spike = find( dv_HP ./ dv_LP > thresh);
        warning(s);  % reset warning state
        
        % Save interm variables for generating plots.
        dv_LP = dv_LP(range);
        dv_HP = dv_HP(range);
        
        % ignore spikes detected in the padding
        spike(spike<range(1) | spike>range(end)) = [];
        
        if isempty(spike) % get out if no spikes are found (RGL)
            dv = dv(range);   % remove the padding
            return
        end
        
        good_points = true(size(dv));
        
        R = round(Fs /(4*smooth)); % arbitrarily decided by Rolf 2015-11-04
        %R = round(Fs/35);  % R for REGION - area size to average when calculating
        % replacement values. Calculated from the sampling rate.
        
        
        %%%% First pass, mark invalid data.
        for s = spike'
            idx = round( max(1,s-N) : min(length(good_points),s+2*N) );
            good_points(idx) = false;
        end
        
        
        %%%% Second pass, find start / stop points.
        start_points = find(diff([true; good_points]) == -1);
        end_points   = find(diff(good_points) == 1);
        start_stop = [start_points end_points];
        
        
        %%%% Third pass, calculate and replace bad points.
        for index = start_stop'
            start = index(1);
            stop  = index(2);
            
            % Find points within region R that are valid and use to
            % calcualte the start value.
            idx = max(start-R,1):start-1;
            points = dv( find(good_points(idx)) + idx(1) - 1 );
            start_value = sum(points) / length(points);
            
            % Find points within region R that are valid and use to
            % calcualte the last value.
            idx = start+1 : min(start+R,length(good_points));
            points = dv( find(good_points(idx)) + idx(1) - 1 );
            stop_value = sum(points) / length(points);
            
            dv(start:stop) = (start_value + stop_value) / 2;
        end
        
        
        % remove the padding
        dv = dv(range);
        spike = spike - padLen;
    end



    function make_figure()
        figure(1000)
        clf
        
        t_junk = (0:length(dv)-1)'/Fs; % Used for plotting
        
        subplot(2,1,1)
        title_string = [...
            'pass-count = ' num2str(pass_count) ...
            ', thresh = ' num2str(thresh) ...
            ', smooth = ' num2str(smooth) ...
            ', Fs = ' num2str(Fs) ...
            ', N = ' num2str(2*N) ...
            ', spikes = ' num2str(length(spike_pass))];
        h = semilogy(...
            t_junk, dv_HP, ...
            t_junk, dv_LP, ...
            (spike_pass - 1)/Fs, dv_HP(spike_pass), '*'); grid on
        set(h(2),'linewidth',3)
        xlabel('\it t \rm [s]')
        y_lim = [0.01 max(abs(dv))];
        set(gca,'ylim', y_lim)
        legend('|sh|', '|sh| - LP', 'Spike', 'location', 'northeast')
        title(title_string)
        
        subplot(2,1,2)
        plot(t_junk, [dv y]); grid on
        xlabel('\it t \rm [s]')
        legend('Dirty', 'Clean', 'location','northeast')
        pause
    end

end