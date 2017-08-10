function corrected_shear = shear_preanalysis(shear, w, pitch, roll, p, sec)

% function corrected_shear = shear_preanalysis(shear, w, ax, ay, p, sec)
% 
% This function clean and prepare shear data for spectral analysis
%
% Where
%   - shear: is the raw shear profile
%   - w: is the falling speed
%   - pitch: the X-angle of the profiler
%   - roll: the Y-angle of the profile
%   - p: is the pressure profile
%   - sec: is the running-mean length (higher it is, more you loose data
%           but more you are sure to be in free falling)
%
% ** You can uncomment the visual inspection plot and try different running
% mean time.
%
%
% Author: Frederic Cyr - 2009/12/17
%   Contributions from D. Bourgault
%
% Modifications:
%                - F. cyr (2010/05/25): inverse 1st and 2nd corrections,
%                despike was not working well when applied after 1st...
%% -------------------------------------------------------------- %%

% --------------------------- %
% fisrt correction (despike.m)
% --------------------------- %

% Parameters for despiking

thresh_spike = 5; smooth_spike = 0.5; N_spike = 5; fs=512;
%thresh_spike = 2; smooth_spike = 0.05; N_spike = 5; fs=512;


% Despiking
[s1, spike] = despike(shear, thresh_spike, smooth_spike, fs, N_spike); % does only despiking 

% Pop-up if despiking method looks bizzare!
if length(spike)>length(shear)/4 % if more than 25% of the profile as been despiked (could be changed...)
    plot(shear, p);
    hold on
    plot(s1, p, 'r')
    set(gca, 'ydir', 'reverse');
    hold off
    legend('original', 'despiked')
    disp('Problem with despiking???');    
    disp('Press any key to continue');    
    pause
end


% $$$ % --------------------------------------------------- %
% $$$ % second correction (home-made method on acceleration)
% $$$ % --------------------------------------------------- %
% $$$ 
% $$$ ddwddz=diff(w); %falling speed derivative (acceleration of the probe)
% $$$ ddwddz_m = runmean(abs(ddwddz), sec*512); %running mean on ?sec.
% $$$ 
% $$$ 
% $$$ % $$$ % --- test --- %
% $$$ % $$$ % $$$ Fs=512; % sample frequency
% $$$ % $$$ Fc = .45*Fs; %cutoff frequency
% $$$ % $$$ Fc2 = .40*Fs; %cutoff frequency
% $$$ % $$$ [B,A]=butter(4,0.001); 
% $$$ % $$$ signal_filt=filtfilt(B,A,ddwddz_m);
% $$$ % $$$ [B,A]=butter(4,0.0001); 
% $$$ % $$$ signal_filt2=filtfilt(B,A,ddwddz_m);
% $$$ % $$$ % ------------ %
% $$$ 
% $$$ %I=find(ddwddz_m>0.00005); % Threshold for the smoothed acceleration (based on visual inspection)
% $$$ I=find(ddwddz_m>0.0005); % Threshold for the smoothed acceleration (based on visual inspection)
% $$$ 
% $$$ s2=s1; 
% $$$ 
% $$$ p2=p;
% $$$ s2(I)=NaN; %flag data wich doesnt satisfies the preceeding condition
% $$$ p2(I)=NaN;
% $$$ 
% $$$ % $$$  % Visual inspection of the home-made correction
% $$$ % $$$  figure(1)
% $$$ % $$$  clf
% $$$ % $$$  plot(s1, p);
% $$$ % $$$  hold on
% $$$ % $$$  plot(s2, p2, 'r')
% $$$ % $$$  set(gca, 'ydir', 'reverse');
% $$$ % $$$  hold off
% $$$ % $$$  % pause(0.2)
% $$$ 
% $$$ % $$$ 
% $$$ % $$$ % --- another test --- %
% $$$ % $$$ ddwddz_m2 = runmean(ddwddz, 5*512); %running mean on ?sec.
% $$$ % $$$ ddwddz_m1 = runmean(ddwddz, 10*512); %running mean on ?sec.
% $$$ % $$$ 
% $$$ % $$$ dif = abs(ddwddz_m2 - ddwddz_m1); 
% $$$ % $$$ 
% $$$ % $$$ I = find(dif>5e-6 & p(1:end-1)>20);
% $$$ % $$$ s2(I)=NaN; %flag data wich doesnt satisfies the preceeding condition
% $$$ % $$$ p2(I)=NaN;
% $$$ % $$$ % ------------ %

% --------------------------- %
% third correction (pitch/roll)
% --------------------------- %

s2 = s1; % uncomment to skip 2nd correction

% Flag data for which pitch or roll > 10  % written by D. Bourgault
%ax = ang2acc(pitch); ay = ang2acc(roll);
J = find(abs(pitch) > 15 | abs(roll) > 15);
s2(J)= NaN;


% returning corrected shear
corrected_shear = s2;
