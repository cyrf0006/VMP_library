function x_cal = calibrate_micro(x, x_ref, W, fs, FS, sep_length, fit_order);

% usage:  x_cal = calibrate_micro(x, x_ref, W, fs, FS, sep_length, fit_order);
%
% x:          Microscale vector that needs to be calibrated.
% x_ref:      Finescale vector against which calibration is done.
% W:          Finescale fall-speed of the VMP
% fs_fast:    Sampling frequency of the microscale data.
% fs_slow:    Sampling frequency of the finescale data.
% sep_length: Separation length between the micro-scale sensors and 
%             the fine-scale sensors [m]. 
%             sep_length = 0.325 for the temperature sensors (found by trial an error)
%             sep_length = 0.550 for the salinity sensors.
% fit_order:  Order of the polynomial fit use to do the calibration
%             (normally 1 or 2)
%
%
% Author: Daniel Bourgault
%
%     modified by F. Cyr - 2009/12/15: 
%             - re-wrote the whole function, using a more compact form, but
%             the method remains the same (removed the conditions on the speed, 
%             worthless bec. used wihtin an identified profile, changed the 
%             expressions of the ratio and x_dec, etc.)
%             - try to add a correction of the distance between
%             probes due to thepitch/roll of the profiler, but
%             realized that the error is    
%             <5mm at an angle of 10 deg.
% -------------------------------------------------------------------- %
    
  ratio = round(fs / FS); % 512/64 = 8
  
  % Decimate the microscale data to match the response of the finescale data. 
  x_dec = mean(reshape(x,ratio,length(x)./ratio))'; % takes tendhe mean value of a matrix 8 X length(x)/8 over the 8 lines

  % Find the lag between the finescale data and the decimated microscale data. 
  for i=1:length(x_ref);
    
    lag = round(sep_length*FS/W(i)); % lag in time step: 1sec. = 64 time steps... time is d/W and lag is time*freq...
      
    if lag > 0 & i-lag > 0 %& W(i) >= Wmin & W(i) <= Wmax %lag represents the advance that has x on x_ref. (i-lag is negative at the beginning of the profile)
        x_dec_lag(i) = x_dec(i-lag);
    else
        x_dec_lag(i) = NaN;
    end
    
  end
  x_dec_lag = x_dec_lag';
  
  % Find the index of the non-nan data points. 
  k_good = find(isnan(x_dec_lag) == 0 & isnan(x_ref) == 0);
  
  % Look at the correlationm coefficient. Should be very high.
  corrcoef(x_dec_lag(k_good),x_ref(k_good));
  
  % Calibrate with a polynomial fit.
  p = polyfit(x_dec_lag(k_good),x_ref(k_good),fit_order);
  x_cal = polyval(p,x);
      
  % Examine visually the relationship between the finescale data (x_ref),
  % the lagged and decimated microscale data and the fit. 
%   plot(x_ref,x_dec_lag,'.');
%   hold on;
%   plot(x_cal,x,'r');
%   hold off;

