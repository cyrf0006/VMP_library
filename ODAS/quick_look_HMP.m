%% quick_look_HMP
% Generation of plots and figures from profiles recorded with a horizontal 
% profiler.
%%
% <latex>\index{Type B!quick\_look\_HMP}</latex>
%
%%% Syntax
%   quick_look_HMP( fname, P_start, P_end )
%
% * [fname] Name of the file to be processed.  The file extension is optional.
% * [P_start] Start point of the data section used for spectral estimates.  This
%             integer value is a depth value.
% * [P_end] End point of the data section used for spectral estimates.
%
%%% Description
%
% A function for the evaluation of a horizontal profiler using an actual
% profile in the ocean. It is usually used during sea trials of an instrument.
% It produces six figures of various channels including spectra of the shear
% probes, thermistors and micro-conductivity sensors. Clients are encouraged to
% examine this file to see one method of data processing. It is nearly
% equivalent to quick_look, so please look at quick_look for guidance on how to
% use this function. The data for all figures below are courtesy of Tim Boyd.
%
% See quick_look.m for more information.
%
% @image @images/quick_look_HMP_2.pdf @quick_look_HMP - example #1 @The first figure produced by 
% quick_look_HMP. It shows the time series of pressure for the entire file 
% (upper panel) and the pressure for the selected time range (lower panel).
%
% @image @images/quick_look_HMP_1.pdf @quick_look_HMP - example #2 @The second figure produced by 
% quick_look_HMP. It shows the pitch, roll, and 100 times the rate of change of 
% pressure (vertical velocity) for the selected time range.
%
% @image @images/quick_look_HMP_3.pdf @quick_look_HMP - example #3 @The third figure produced by 
% quick_look_HMP. It shows the temperature for the selected time range using 
% only the nominal linear conversion to physical units.
%
% @image @images/quick_look_HMP_4.pdf @quick_look_HMP - example #4 @The fourth figure produced by 
% quick_look_HMP. It shows the gradient of the microstructure signals for the 
% selected time range. BP indicates low-pass filtered at 30 Hz.
%
% @image @images/quick_look_HMP_5.pdf @quick_look_HMP - example #5 @The fifth figure produced by 
% quick_look_HMP. It shows the spectra of acceleration and all microstructure 
% signals for the selected time range. Shear probe 1 and thermistor 1 were not 
% installed. Note the propeller vibrations, evident in the acceleration spectra.
%
% @image @images/quick_look_HMP_6.pdf @quick_look_HMP - example #6 @The sixth figure produced by 
% quick_look_HMP. It shows the wavenumber spectra of the two shear probes. The 
% thin green and red lines are the spectra of shear. The thick lines are the 
% spectra after removal of acceleration-coherent vibrations and correction for 
% the spatial response of the shear probes. The dotted curves are polynomial 
% fits to determine the upper limit of spectral integration to obtain the 
% variance of shear. The thin black lines are the Nasmyth spectrum for the rate 
% of dissipation obtained from the shear variance. The triangles indicate the 
% upper limit of spectral integration. Shear probe 1 and thermistor 1 were not 
% installed. Note that the propeller vibrations have been removed from the shear
% probe spectrum (thin versus thick green curves).

% *Revision history:*
% * 2009-03-31 (RGL) customization for the Remus-600 AUV
% * 2010-04-28 (AWS) changes related to odas v6
% * 2011-09-01 (AWS) added documentation tags for matlab publishing
% * 2012-11-09 (WID) documentation update
% * 2012-11-09 (WID) removed inifile_withstr and save_rolf commands.  Use 
%                    fig2pdf where applicable.
% * 2013-02-26 (WID) code for converting channels taken from VMP function

function quick_look_HMP(fname, t_start, t_end)

fft_length = 2; % length in seconds
HP_cut = 0.5; % Frequency (Hz) of high-pass filter cut-off for shear probes.
LP_CUT = 20; % Frequency (Hz) of low-pass filter for shear signals. 
% It is used only for profile display and is not used in spectral calculations. 
sh1_sens      = 0.1304; % shear probe SN474
sh2_sens      = 0.1286; % SN423
P_diff_gain   = 20.5; % Gain of presure pre-emphasis ~20.5
T1_diff_gain  = 1.0; % T1 differentiator gain ~1.0
T2_diff_gain  = 1.0; % T2 differentiator gain ~1.0
Sh1_diff_gain = 1.03; % Shear Channel 1 differentiator gain ~1.0
Sh2_diff_gain = 1.025; % Shear Channel 2 differentiator gain ~1.0
ADC_FS = 5; % Full-scale voltage of A-to-D converter
ADC_bits = 16; % Number of bits in A-to-D converter
speed = 1.0; % This is the speed of the AUV which is temporarily set to a fixed value.
% ____________________________________________________________
% File opening etc.
% Check file name, check for a .mat file, and open the .p file if no *.mat
% file
% _________________________________________________________________

% Ask the user for a file name if not provided as an input variable



we_have_old_mat_file = 0;% if mat-file already exists, assume it is in physical units.

% Look for the .p and .mat files, deal with problems
[P N E] = file_with_ext( fname, ...
                         {'.p', '.P', ''}, ...
                         ['Could not find file: ' fname]);

if exist([N '.mat'],'file') % Then the mat-file already exist and we will
    % use it and go straight to the plotting.
    disp(['Loading file: ' N '.mat']);
    load([N '.mat'])
    we_have_old_mat_file = 1;
else
    if isempty(read_odas([N E])),
        error('Could not convert data into a mat-file');
    end
    disp(['Loading file  = ' N '.mat']);
    load([N '.mat']); % use the one we just converted
end

if header_version >= 6,
    cfg = setupstr(setupfilestr);
    if ~isempty(setupstr(cfg, '', 'xmp'))
        error('please run quick_look_XMP.m');
    end
end

if header_version >= 6
    tmp = setupstr(cfg, 'shear1', 'sens');
    sh1_sens = str2double(tmp{1});
    tmp = setupstr(cfg, 'shear1', 'diff_gain');
    Sh1_diff_gain = str2double(tmp{1});

    tmp = setupstr(cfg, 'shear2', 'sens');
    sh2_sens = str2double(tmp{1});
    tmp = setupstr(cfg, 'shear2', 'diff_gain');
    Sh2_diff_gain = str2double(tmp{1});
    
    tmp = setupstr(cfg, 'ucond1', 'diff_gain');
    C1_diff_gain = str2double(tmp{1});
    tmp = setupstr(cfg, 'ucond2', 'diff_gain');
    C2_diff_gain = str2double(tmp{1});
end

% __________________________________________________________________
%
% This is where we get some information about the data sampling

%header = header(:); header = reshape(header,64, length(header)/64)';
if (~we_have_old_mat_file), %We do not have a file that is already in physical

    t_slow = (0:length(P)-1)'/fs_slow; % make time vectors for plotting
    t_fast = (0:length(sh1)-1)'/fs_fast;
        
    Year  = header(2,4);
    Month = header(2,5);
    Day   = header(2,6);
    Hour  = header(2,7);
    Minute= header(2,8);
    if header_version >= 6,
        % header_version >= 6 % We have a new style data file that has a configuration string
        % Deconvolve the channels with pre-emphasis
        process_list = {{ 'T1_dT1', 'T1_hres'}, ...
                        { 'T2_dT2', 'T2_hres'}, ...
                        { 'C1_dC1', 'C1_hres'}, ...
                        { 'C2_dC2', 'C2_hres'}};
        
        for item = process_list,
            set = item{1};
            
            % Skip item if it is not a variable - start again with the next item
            if ~exist(set{1}, 'var'), continue, end
            
            % Deconvole the item.  tmp used because it simplifies the eval
            % statement
            tmp = deconvolve( set{1}, [], eval(set{1}), fs_fast, ...
                              setupfilestr, header_version );
            eval([set{2} ' = tmp;']);
        end
        
        % Pressure is special
        P_hres  = deconvolve('P_dP', P, P_dP, fs_slow, setupfilestr, header_version);
        % end of deconvolution
        %
        % Now convert to physical units
        convert_list = {'T1_hres', 'T2_hres', 'C1_hres', 'C2_hres', 'P', 'P_hres', ...
            'PTV', 'sbt', 'sbc', 'Ax', 'Ay', 'Az', 'Rx', 'Ry', 'Rz', 'Mx', 'My', 'Mz', ...
            'U', 'Ux', 'Uy', 'V', 'W', 'V_Bat', 'Incl_X', 'Incl_Y', 'Incl_T','sh1', 'sh2'};
        name_list = convert_list;
        if exist('T1','var'), name_list{1} = 'T1'; else name_list{1} = 'T1_dT1';end;% use coefficients for non-differentiated channel
        if exist('T2','var'), name_list{2} = 'T2'; else name_list{2} = 'T2_dT2';end
        name_list{3} = 'C1_dC1';
        name_list{4} = 'C2_dC2';
        name_list{6} = 'p';
        
        % name_list gives the names of the channels that contain the conversion
        % coefficients. For example, "P" is a variable and is found in the
        % configuration file, but P_hres is not in the coefficient file. The
        % coefficients for P_hres are in the section with the name "P".
        % Similarly for "T1_hres", etc. Most of the variables in the convert_list
        % have names that are identical to those in the configuration file.
        
        % Note that sh1 and sh2 (shear probes) still must be divided by speed
        % squared in order to complete the conversion to physical units.
        for ii = 1:length(convert_list),
            item = convert_list{ii};
            name = name_list{ii};
            
            if ~exist(item, 'var'), continue, end
            tmp = convert_odas( eval(item), name, 'string', setupfilestr, header_version );
            eval([name ' = tmp;']);
        end
    else %% header_version 1 (old style)
    %_____________________________________
    % Deconvolve the channels with pre-emphasis
        P_hres = deconvolve('pressure', P, P_dP, fs_slow, P_diff_gain);
        T1_hres = deconvolve('temperature', [], T1_dT1, fs_fast, T1_diff_gain);
        T2_hres = deconvolve('temperature', [], T2_dT2, fs_fast, T2_diff_gain);
        T1_hres = (ADC_FS / 2^ADC_bits)*T1_hres*15 + 15; % +/-1V for +/-15C
        T2_hres = (ADC_FS / 2^ADC_bits)*T2_hres*15 + 15; % +/-1V for +/-15C

        if exist('C1_dC1','var'), 
            C1_hres = deconvolve('conductivity', [], C1_dC1, fs_fast, 0.99);
            C1_hres = (ADC_FS / 2^ADC_bits)*C1_hres*30 + 30; % +/-1V for +/-30 mS/cm.
        end
        if exist('C2_dC2','var'), 
            C2_hres = deconvolve('conductivity', [], C2_dC2, fs_fast, 0.99);
            C2_hres = (ADC_FS / 2^ADC_bits)*C2_hres*30 + 30; % +/-1V for +/-30 mS/cm.
        end
        %___________________________________________________
        % Convert to physical units
        
        P_hres = convert_odas(P_hres, 'Pres', 'file', 'setup.txt');
        P      = convert_odas(P,      'Pres', 'file', 'setup.txt');
        
        if exist('sbt','var'), sbt = convert_odas (sbt,'SBT1E', 'file', 'setup.txt'); end
        if exist('sbc','var'), sbc = convert_odas (sbc,'SBC1E', 'file', 'setup.txt'); end
        
        Ax  = convert_odas (Ax,'Pitch', 'file', 'setup.txt');
        Ay  = convert_odas (Ay,'Roll' , 'file', 'setup.txt');
        if exist('Az','var'), Az  = convert_odas (Az,'az'   , 'file', 'setup.txt'); end
        
        if exist('Mx','var'), Mx = convert_odas (Mx,'Mx', 'file', 'setup.txt'); end
        if exist('My','var'), My = convert_odas (My,'My', 'file', 'setup.txt'); end
        if exist('Mz','var'), Mz = convert_odas (Mz,'Mz', 'file', 'setup.txt'); end
        
    end
        
    % Generate pressure and vertical speed vectors
    % 0.5 Hz seems to give the best smoothing for normal profiling. Adjust,
    % as needed.
    W_slow = gradient(P_hres, 1/fs_slow);
    [b,a] =butter(4,0.5/fs_slow/2);
    W_slow = filtfilt(b,a,W_slow);
    P_fast = interp1(t_slow, P_hres,t_fast);
    W_fast = interp1(t_slow, W_slow,t_fast);
    
    % Shear Probes
    if (header_version >= 6) %% header_version 6 and up
        tmp = setupstr(cfg, 'shear1', 'adc_fs');
        ADC_FS = str2double(tmp{1});
        tmp = setupstr(cfg, 'shear1', 'adc_bits');
        ADC_bits = str2double(tmp{1});
    end
    sh1 = sh1.*(ADC_FS / 2^ADC_bits)./(2*sqrt(2)*Sh1_diff_gain*sh1_sens*speed.^2);

    if (header_version >= 6) %% header_version 6 and up
        tmp = setupstr(cfg, 'shear2', 'adc_fs');
        ADC_FS = str2double(tmp{1});
        tmp = setupstr(cfg, 'shear2', 'adc_bits');
        ADC_bits = str2double(tmp{1});
    end
    sh2 = sh2.*(ADC_FS / 2^ADC_bits)./(2*sqrt(2)*Sh2_diff_gain*sh2_sens*speed.^2);
    
    % Thermistors
    gradT1 = gradient(T1_hres,1/fs_fast);
    gradT2 = gradient(T2_hres,1/fs_fast);
    gradT1 = gradT1 ./ speed; % in units of degree C per meter
    gradT2 = gradT2 ./ speed;
    
    % Micro-conductivity
    if exist('C1_dC1','var'), 
        gradC1 = gradient(C1_hres,1/fs_fast);
        gradC1 = gradC1 ./ speed; % in units of mS/cm per meter
    end
    if exist('C2_dC2','var'),
        gradC2 = gradient(C2_hres,1/fs_fast);
        gradC2 = gradC2 ./ speed;
    end

% Save the converted variables to the mat-file
eval (['save ' N ' fs_fast fs_slow Year t_slow t_fast Month Day Hour Minute -append -v6'])
flag = save_odas([N '.mat'], 'P_hres',  P_hres);  if (flag < 0), error(['Could not write P_hres to '  N '.mat']), end
flag = save_odas([N '.mat'], 'P',       P);       if (flag < 0), error(['Could not write P to '  N '.mat']), end
flag = save_odas([N '.mat'], 'T1_hres', T1_hres); if (flag < 0), error(['Could not write T1_hres to ' N '.mat']), end
flag = save_odas([N '.mat'], 'T2_hres', T2_hres); if (flag < 0), error(['Could not write T2_hres to ' N '.mat']), end
if exist('sbt','var');
    flag = save_odas([N '.mat'], 'sbt', sbt); 
end
if exist('sbc','var');
    flag = save_odas([N '.mat'], 'sbc', sbc); 
    if (flag < 0), error(['Could not write sbc to ' N '.mat']), end, 
end
flag = save_odas([N '.mat'], 'Ax', Ax); if (flag < 0), error(['Could not write Ax to ' N '.mat']), end
flag = save_odas([N '.mat'], 'Ay', Ay); if (flag < 0), error(['Could not write Ay to ' N '.mat']), end
flag = save_odas([N '.mat'], 'Az', Az); if (flag < 0), error(['Could not write Az to ' N '.mat']), end

flag = save_odas([N '.mat'], 'W_slow', W_slow); if (flag < 0), error(['Could not write W_slow to ' N '.mat']), end
flag = save_odas([N '.mat'], 'W_fast', W_fast); if (flag < 0), error(['Could not write W_fast to ' N '.mat']), end
flag = save_odas([N '.mat'], 'P_fast', P_fast); if (flag < 0), error(['Could not write P_fast to ' N '.mat']), end

flag = save_odas([N '.mat'], 'sh1', sh1); if (flag < 0), error(['Could not write sh1 to ' N '.mat']), end
flag = save_odas([N '.mat'], 'sh2', sh2); if (flag < 0), error(['Could not write sh2 to ' N '.mat']), end

flag = save_odas([N '.mat'], 'gradT1', gradT1); if (flag < 0), error(['Could not write gradT1 to ' N '.mat']), end
flag = save_odas([N '.mat'], 'gradT2', gradT2); if (flag < 0), error(['Could not write gradT2 to ' N '.mat']), end

if exist('gradC1','var')
    flag = save_odas([N '.mat'], 'gradC1', gradC1); 
    if (flag < 0), error(['Could not write gradC1 to ' N '.mat']), end
    flag = save_odas([N '.mat'], 'C1_hres', C1_hres); 
    if (flag < 0), error(['Could not write C1_hres to ' N '.mat']), end
end
if exist('gradC2','var')
    flag = save_odas([N '.mat'], 'gradC2', gradC2); 
    if (flag < 0), error(['Could not write gradC2 to ' N '.mat']), end
    flag = save_odas([N '.mat'], 'C2_hres', C2_hres); 
    if (flag < 0), error(['Could not write C2_hres to ' N '.mat']), end
end
if exist('Mx','var')
    flag = save_odas([N '.mat'], 'Mx', Mx); 
    if (flag < 0), error(['Could not write Mx to ' N '.mat']), end
end
if exist('My','var')
    flag = save_odas([N '.mat'], 'My', My); 
    if (flag < 0), error(['Could not write My to ' N '.mat']), end
end
if exist('Mz','var')
    flag = save_odas([N '.mat'], 'Mz', Mz); 
    if (flag < 0), error(['Could not write My to ' N '.mat']), end
end
end

% This is the end of the section for converting to physical units. This
% section was skipped if the mat-file already existed. If the mat-file is
% corrupted or has a bad conversion, then erase it and run this function
% again

% _______________________________________________________________________
% Make figures
% Pitch, Roll Fall-Rate

n = find(t_fast> t_start & t_fast < t_end); % Creat index for section of profile to be processed.
m = find(t_slow> t_start & t_slow < t_end);

if isempty(n)
    error('No profile detected')
    return
end

if(Minute < 10)
    Minute_string = ['0' num2str(Minute)];
else
    Minute_string = num2str(Minute);
end

date_string = [num2str(Year) '\_' num2str(Month) '\_' num2str(Day) ...
    ', ' num2str(Hour) ':' Minute_string];

%---------------------------------------------------------------
figure(10); clf
subplot(2,1,1)
h = plot(t_slow, P_hres);grid
ylabel('\it P \rm [dBar]','fontsize',16)
set(h(1), 'linewidth', 2,   'color','b')
title([ fix_underscore(N) ';  ' date_string], 'fontsize', 16)
set(gca,'fontsize', 16)

subplot(2,1,2)
h = plot(t_slow(m), P_hres(m)); grid
set(gca,'fontsize', 16)
set(h(1), 'linewidth', 2,   'color','b')
xlabel('\it t \rm [s]','fontsize',16)
ylabel('\it P \rm [dBar]','fontsize',16)
orient landscape

%file_name = [N '_' num2str(t_start) '_Fig_10'];
%fig2pdf( gcf, file_name );


figure (1); clf
if exist('Mx','var')
    h = plot(asin([Ax(n),Ay(n)]./9.81)*180/pi, t_fast(n), [W_slow(m) 1.8*angle(Mx(m)+i*My(m))/pi], t_slow(m)); grid
    legend('Pitch', 'Roll', 'W', 'Compass/100','location','eastoutside') 
else
    [b,a] = butter(4,2/(fs_fast/2)); % low-pass to eliminate vibrations
%     h = plot(asin([Ax(n),Ay(n)]./9.81)*180/pi, t_fast(n), ...
%         filtfilt(b,a,asin([Ax(n),Ay(n)])./9.81)*180/pi, t_fast(n),  ...
%         W_slow(m), t_slow(m)); grid
    h = plot(t_fast(n), filtfilt(b,a,asin([Ax(n) Ay(n)]./9.81))*180/pi,  ...
        t_slow(m), 100*W_slow(m)); grid
%     set(h(1), 'linewidth', 0.5, 'color','b')
%     set(h(2), 'linewidth', 0.5, 'color',[0 0.5 0])
    set(h(1), 'linewidth', 2,   'color','b')
    set(h(2), 'linewidth', 2,   'color','g')
    set(h(3), 'linewidth', 2, 'color','r')
    legend('Pitch', 'Roll','100\timesd\itP\rm/d\itt', 'Location', 'eastoutside')
end
set(gca,'fontsize', 16)
xlabel('\it t \rm [s]','fontsize',16)
ylabel('[  ^{\circ} ,  m s^{-1}]','fontsize',16)
title([ fix_underscore(N) ';  ' date_string], 'fontsize', 16)
orient landscape

%file_name = [N '_' num2str(t_start) '_Fig_01'];
%fig2pdf( gcf, file_name );

%-----------------------------------------------------------------------
%Sea-Birds
if (exist('sbt','var') && exist('sbc','var'))
    figure(2); clf
    sbt(m) = despike(sbt(m), 4, 1, fs_slow, 5);
    S = salinity(P_hres,sbt,sbc);
    S_0 = mean(S(m(1:1000))); % Near-surface average salinity
    h = plot([sbt(m) 0.2218+500*(sbc(m)-0.2218) 500*(S(m)-S_0)], t_slow(m));grid
    set(gca, 'ydir', 'rev','xlim',[-3 8],'fontsize',16)
    ylabel('\it P \rm [dBar]','fontsize',16)
    xlabel('[  ^{\circ}C, mS cm^{-1}, PSU]','fontsize',16)
    hh=legend('\itT_{SB}', '0.2218+500(\itC_{SB}\rm - 0.2218)', ['500(\itS\rm - ' num2str(S_0) ')'], 'Location','eastoutside');
    set(hh,'fontsize', 14)
    title([ fix_underscore(N) ';  ' date_string], 'fontsize',16)
    orient landscape

    %file_name = [N '_' num2str(t_start) '_Fig_02'];
    %fig2pdf( gcf, file_name );
end

%----------------------------------------------------------
%Thermistors
% Remeber that sbt_offset is the distance between FP07 and SBE3 in meters
%mm = find(t_slow> P_min+sbt_offset & W_slow > W_min);
% This lines up the pressure range of the Sea-Bird and the thermistors
%mm = mm(1:end-8*64);% trim off bottom 8 seconds (~1m)

figure (3);clf
if exist('sbt','var') % The first part is a rough but reliable way of scaling the thermistors
    sbt_mean = mean(sbt(mm)); % mean actual temperature
    T1_hres_mean = mean(T1_hres(n));
    T2_hres_mean = mean(T2_hres(n));
    T1_offset = T1_hres_mean - sbt_mean; % calculate the approximate offset
    T2_offset = T2_hres_mean - sbt_mean; % calculate the approximate offset
    T1_hres = T1_hres-T1_offset; % This should align the thermistor with the SBT values
    T2_hres = T2_hres-T2_offset; 

    % The next part is a more formal way to calibrate the thermistors
    % against the Sea-Bird. However, nothing is perfect because the
    % fall-rate is not constant. Perhaps a still better way to do this is
    % to re-interpolate the data to evenly spaced depth intervals; that
    % is cumbersome and I leave it to others.
    Ratio = floor(length(n)/length(mm)); % Should be almost exactly 8
    nn = n(1:Ratio*length(mm)); % trim indeces in n to make them a whole number multiple of length(m)
    T1_junk = reshape(T1_hres(nn),Ratio, length(nn)/Ratio);
    T1_junk = mean(T1_junk)';
    
    p_T1 = polyfit(T1_junk, sbt(mm), 3);
    T1_hres = polyval(p_T1, T1_hres);
    p_T2 = polyfit(T2_hres(n), T1_hres(n), 3);
    T2_hres = polyval(p_T2, T2_hres);
    h = plot([sbt(m)], P_hres(m)-sbt_offset, [T1_hres(n) T2_hres(n)], t_fast(n));grid
    legend('T_{SB}', 'T1', 'T2', 'location','eastoutside')
else
%    h = plot(t_fast(n), [T1_hres(n) T2_hres(n)]);grid
    h = plot(t_fast(n), [T2_hres(n)]);grid
%    legend('T1', 'T2', 3)
    legend('T2', 'location', 'eastoutside')
end
set(gca,'fontsize', 16)
xlabel('\it t \rm [s]','fontsize',12)
ylabel('[  ^{\circ}C ]','fontsize',12)
title([ fix_underscore(N) ';  ' date_string], 'fontsize',16)
orient landscape

%file_name = [N '_' num2str(t_start) '_Fig_03'];
%fig2pdf( gcf, file_name );

%--------------------------------------------------------------------------
% Temperature gradient, microconductivity gradient, and shear
figure (4);clf

sh1 = sh1.*(speed.^2); % This removes the contribution from fall-rate
sh1 = sh1 - mean(sh1(n));
sh2 = sh2.*(speed.^2);
sh2 = sh2 - mean(sh2(n));

sh1(n) = despike(sh1(n), 7, 1, fs_fast, 20);
sh2(n) = despike(sh2(n), 7, 1, fs_fast, 20);

% Apply gentle high-pass filter
[b,a] = butter(1,HP_cut/(fs_fast/2),'high'); 
sh1_HP = sh1;
sh2_HP = sh2;
sh1_HP(n) = filtfilt(b,a,sh1(n));
sh2_HP(n) = filtfilt(b,a,sh2(n));

sh1 = sh1 ./ speed.^2;% This scales the signal back into physical units
sh2 = sh2 ./ speed.^2;
sh1_HP = sh1_HP ./ speed.^2;
sh2_HP = sh2_HP ./ speed.^2;

[b,a] = butter(4, LP_CUT/(fs_fast/2));
sh1_BP = sh1_HP;
sh2_BP = sh2_HP;
sh1_BP(n)=filtfilt(b,a,sh1_BP(n)); % High- and low-pass (band-pass) for plotting only
sh2_BP(n)=filtfilt(b,a,sh2_BP(n));

h = plot(t_fast(n), [gradT1(n)-10 gradT2(n)-5  sh1(n) sh2(n)+5 ], ...
    t_fast(n), [sh1_BP(n)+10 sh2_BP(n)+15 ], t_fast(n),T2_hres(n));grid
hh=legend('\partial T_1/\partial z', '\partial T_2/\partial z', '\partial u_1/\partial z',...
    '\partial u_2/\partial z',...
    '\partial u_1/\partial z BP','\partial u_2/\partial z BP','T_2', 'location','eastoutside');

set(hh,'fontsize', 14)
set(gca,'fontsize',16)
xlabel('\it t \rm [s]')
ylabel('[  ^{\circ}C m^{-1},  s^{-1} ]')
title([ fix_underscore(N) ';  ' date_string], 'fontsize',16)
orient landscape

%file_name = [N '_' num2str(t_start) '_Fig_04'];
%fig2pdf( gcf, file_name );

%----------------------------------------------------------------
%spectra
figure(5); clf
fft_num = round(fft_length*fs_fast);% fft length in units of samples

%n = find(((t_fast > t_start) & (t_fast < t_end)) & (W_fast > W_min));
[P_Ax, F]  = psd_rolf( Ax(n), fft_num, fs_fast);
[P_Ay, F]  = psd_rolf( Ay(n), fft_num, fs_fast);
[P_Az, F]  = psd_rolf( Az(n), fft_num, fs_fast);
[P_sh1_HP, F] = psd_rolf(sh1_HP(n), fft_num, fs_fast);
[P_sh2_HP, F] = psd_rolf(sh2_HP(n), fft_num, fs_fast);
[P_gradT1, F] = psd_rolf(gradT1(n), fft_num, fs_fast);
[P_gradT2, F] = psd_rolf(gradT2(n), fft_num, fs_fast);

    epsilon = 1e-10*[1 10 100 1e3 1e4]; % 
    nu = visc35(mean(T1_hres(n)));
    [Pn, kn]=nasmyth(epsilon, nu);
    fn = kn*mean(speed); Pn = Pn/mean(speed);

    h = loglog(F,[P_Ax P_Ay P_Az P_sh1_HP P_sh2_HP P_gradT1 P_gradT2], ...
        fn, Pn, 'k');grid
    hh=legend('A_x', 'A_y', 'A_z', ...
    'du_1/dz-HP', 'du_2/dz-HP', 'dT_1/dz', 'dT_2/dz', 'Nasmyth','location','eastoutside');
    set(hh,'fontsize', 16)

    set(h(1),'linewidth',3,'color','b'); % Ax
    set(h(2),'linewidth',3,'color',[0 0.5 0]); % Ay
    set(h(3),'linewidth',3,'color','r'); % Az
    set(h(4),'linewidth',2,'color','b'); % sh1_HP
    set(h(5),'linewidth',2,'color',[0 0.5 0]); % sh2_HP
    set(h(6),'linewidth',2,'color','c'); % grad_T1
    set(h(7),'linewidth',2,'color','m'); % grad_T2

    for index = 1: length (epsilon)
            junk = find(fn(:,index)>=1); junk=junk(1);
            text(fn(junk,index), Pn(junk,index), ['10^{' num2str(log10(epsilon(index)),2) '}'], ...
         'fontsize', 16, 'fontweight','bold')
    end

% if (exist('gradC1','var') && exist('gradC2','var'))
%     f_upper = max([F  F_C1 F_C2]); 
% elseif (exist('gradC1','var') && ~exist('gradC2','var'))
%     f_upper = max([F F_C1]);
% else
%     f_upper = F;
% end
% f_upper = f_upper(1);
f_upper = F(end);
set(gca, 'XMinorGrid', 'off','YMinorGrid', 'off','fontsize',16)
set(gca, 'ylim', [1e-9 10], 'xlim',[1./fft_length f_upper])
ylabel('[Variance Hz^{-1}]')
xlabel('\it f \rm [Hz]')
title([ fix_underscore(N) ';  ' date_string ' ; ' num2str(t_start) ' - ' num2str(t_end) ' s'], 'fontsize',16)
orient landscape

%file_name = [N '_' num2str(t_start) '_Fig_05'];
%fig2pdf( gcf, file_name );

%_____________________________
figure(6);clf
[clean_UU1, AA, UU1, UA, F] = clean_shear_spec([Ax(n) Ay(n) Az(n)], ...
   sh1_HP(n), fft_num, fs_fast);
[clean_UU2, AA, UU2, UA, F] = clean_shear_spec([Ax(n) Ay(n) Az(n)], ...
   sh2_HP(n), fft_num, fs_fast);
W = mean(speed);
K = F/W;
clean_UU1 = clean_UU1 * W; UU1 = UU1 * W;
clean_UU2 = clean_UU2 * W; UU2 = UU2 * W;
clean_UU1 = clean_UU1 .* (1 + (K / 48).^2); % Wavenumber correction after Macoun & Lueck
clean_UU2 = clean_UU2 .* (1 + (K / 48).^2);
UU1 = UU1 .* (1 + (K / 48).^2);
UU2 = UU2 .* (1 + (K / 48).^2);

% Now find first minimum of spectrum in log-log space. Look only at
% wavenumbers smaller than 150 cpm. This is a little complicated.
junk = find (K < 150);
Index_limit = length(junk);
y1 = log10(clean_UU1(2:Index_limit));
y2 = log10(clean_UU2(2:Index_limit));
x = log10(K(2:Index_limit));
% Note K(1) is always zero and has problems with logarithms. Now we fit a
% polyninomial to the spectrum

order = 4; % Keep the order reasonable. Any value from 4 to 8 could work
if (order<2), order=3; end % at least third order
p1 = polyfit(x,y1,order);
p1_div = p1(1:end-1); p1_div = p1_div .*[order:-1:1]; % poly coefficients of first derivative
div_zeros = roots (p1_div);
z = find (div_zeros >= 1); div_zeros = div_zeros(z); % only look for minimum above log_10(10 cpm)
div_zeros = sort(div_zeros); % put into ascending order
p1_curve = p1_div(1:end-1); p1_curve = p1_curve .* [order-1:-1:1]; % poly coefficients of second derivative

K1_limit = log10(150); % absolute upper limit is 150 cpm
if ~isempty(div_zeros) % check if its a max or min
    for M = length(div_zeros):-1:1
        if (polyval(p1_curve, div_zeros(M)) > 0); % we found one.
            K1_limit = div_zeros(M); % this is a slope zero.
        end % There might be more, so get the one with the lowest wavenumber
    end
end
if (K1_limit<1); K1_limit = 1; end
if (K1_limit>log10(150)); K1_limit = log10(150); end

p2 = polyfit(x,y2,order); 
p2_div = p2(1:end-1); p2_div = p2_div .*[order:-1:1]; % poly coefficients of first derivative
div_zeros = roots (p2_div); div_zeros = sort(div_zeros);
z = find (div_zeros >= 1); div_zeros = div_zeros(z); % only look for minimum above log_10(10 cpm)
div_zeros = sort(div_zeros); % put into ascending order
p2_curve = p2_div(1:end-1); p2_curve = p2_curve .* [order-1:-1:1]; % poly coefficients of second derivative

K2_limit = log10(150); 
if ~isempty(div_zeros) % check if its a max or min
    for M = length(div_zeros):-1:1
        if (polyval(p2_curve, div_zeros(M)) > 0); % we found one.
            K2_limit = div_zeros(M); % this is a slope zero.
        end % There might be more, so get the one with the lowest wavenumber
    end
end
if (K2_limit<1); K2_limit = 1; end
if (K2_limit>log10(150)); K2_limit = log10(150); end
%K1_limit = 1; K2_limit = 1;

Range_1 = find(K <= 10.^K1_limit+0.5); %Integration range for sh1
Range_2 = find(K <= 10.^K2_limit+0.5); %Integration range for sh1

nu = visc35(mean(T1_hres(n)));
e1= 7.5*nu*trapz(K(Range_1), clean_UU1(Range_1));
[phi1,k1] = nasmyth(e1,nu,512);
e2= 7.5*nu*trapz(K(Range_2), clean_UU2(Range_2));
[phi2,k2] = nasmyth(e2,nu,512);
h = loglog(k1, phi1, '--k', K, [clean_UU1 UU1], ...
    10.^(x), 10.^(polyval(p1,x)),':', ...
    k2, phi2, 'k', K, [clean_UU2 UU2], ...
    10.^(x), 10.^(polyval(p2,x)), ':',...    
    K(Range_1(end)), clean_UU1(Range_1(end)), '^', ...
    K(Range_2(end)), clean_UU2(Range_2(end)), '^') ;grid

set(gca,'xlim',[0.5 200], 'ylim', [1e-9 1e0])
set(h(1),'linewidth',1.5,'color',[0 0 0]) % Nasmith 1
set(h(2),'linewidth',4.0,'color',[1 0 0]) % sh1 clean
set(h(3),'linewidth',1.0,'color',[1 0 0]) % sh1 original
set(h(4),'linewidth',3  ,'color',[1 0 1]) % Polyfit 1
set(h(5),'linewidth',1.5,'color',[0 0 0]) % Nasmyth 2
set(h(6),'linewidth',4.0,'color',[0 0.5 0]) % sh2 slean
set(h(7),'linewidth',1.0,'color',[0 0.5 0]) % sh2 original
set(h(8),'linewidth',3  ,'color',[0 1.0 0]) % polyfit 2
set(h(9), 'MarkerFaceColor','m','MarkerSize',14) % sh1 limit
set(h(10),'MarkerFaceColor','g','MarkerSize',14) % sh2 limit
 
xlabel('\itk \rm [cpm]','fontsize',16)
ylabel('\Phi (\itk\rm)  [s^{-2} cpm^{-1}]','fontsize',16)
hh=legend(['\epsilon_1 = ' num2str(e1,2) ' W kg^{-1}'], ...
    '\partial u_1/\partial z_{clean} HP','\partial u_1/\partial z HP', ...
    'p_1 - fit', ...
    ['\epsilon_2 = ' num2str(e2,2) ' W kg^{-1}'], ...
    '\partial u_2/\partial z_{clean} HP','\partial u_2/\partial z HP', ...
    'p_2 - fit',  'orientation','vertical','location','eastoutside');

set(gca, 'XMinorGrid', 'off','YMinorGrid', 'off','fontsize',16)
set(hh,'fontsize',14)
orient landscape
title([ fix_underscore(N) ';  ' date_string ' ; ' ...
    '\itU\rm=' num2str(mean(speed),2) ' m/s ; ' ...
    num2str(t_start) ' - ' num2str(t_end) ' s'], 'fontsize',16)

%file_name = [N '_' num2str(t_start) '_Fig_06'];
%fig2pdf( gcf, file_name );

clear all




