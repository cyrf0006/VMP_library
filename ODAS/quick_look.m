%% quick_look
% Generation of plots and figures from profiles recorded with a vertical 
% profiler.
%%
% <latex>\index{Type B!quick\_look}</latex>
%
%%% Syntax
%   quick_look( fname, P_start, P_end )
%
% * [fname] Name of the file to be processed.  The file extension is optional.
% * [P_start] Start point of the data section used for spectral estimates.  This
%             integer value is a depth value.
% * [P_end] End point of the data section used for spectral estimates.
%
%%% Description
%
% This function is for the rapid generation of plots and figures from profiles
% recorded with a vertical profiler.  Horizontal profilers should use 
% quick_look_HMP. It is intended for the quick evaluation of an instrument 
% during a sea trial. Clients are encouraged to examine this file to see an 
% example of data processing.
%
% If the depth range for calculating spectra is unknown, try some reasonable
% values for P_start and P_end and re-run the function after examing the 
% resulting figures. Most figures do not depend on these values.
%
% When using legacy data files, the local directory should contain a copy of the
% ODAS setup.txt file that contains calibration coefficients. This function will
% use those coefficients for the conversion to physical units.  This is not 
% required for newer v6 data files as the calibration coefficients are embedded
% within the data file.
%
% The function performs a massive number of jobs. It starts by converting the
% binary data file into a Matlab mat-file and then loads the file. It
% systematically converts most channels into physical units. The following
% coefficients control the behavior of the quick_look.m function.  Non-legacy
% ODAS data files (v6 or higher) only have to worry about the first section of
% variables.
%
%    % ************************************************************************
%    % Play with the following variables to optimize the outputs.
%    % ************************************************************************
%    P_min = 3.0;    % Minimum depth for profile start
%    W_min = 0.4;    % Minimum speed for profile start
%    Rise = 0;       % for rising profilers
%    fft_length = 1; % length in seconds
%    diss_length = 3*fft_length;
%    over_lap = diss_length/2;
%    HP_cut = 0.25;  % High-pass filter cut-off frequency (Hz) for shear probe
%    LP_CUT = 25;    % Low-pass filter cut-off frequency (Hz) for shear signals
%    %
%    % ************************************************************************
%    % Required for processing legacy data files.  The following default values 
%    % should not be changed if processing version 6 data files.
%    % ************************************************************************
%    T1_offset =  0;         % Adjust T1 by this amount
%    T2_offset =  0;         % Adjust T2 by this amount
%    sbt_offset = 0.35;      % distance between FP07 and SBE3 in meters 
%                            %(usually 0.3 m or 2.5 m)
%    pump = 0;               % pump flag: 1 = SBT5 pump installed, 0 otherwise.
%    sh1_sens      = 0.0759; % SN725, 2010-11-25
%    sh2_sens      = 0.0719; % shear probe SN735, 2010-11-25
%    P_diff_gain   = 19.9;   % Gain of presure pre-emphasis ~20.5
%    T1_diff_gain  = 0.99;   % T1 differentiator gain ~1.0
%    T2_diff_gain  = 0.99;   % T2 differentiator gain ~1.0
%    C1_diff_gain = 0.4;     % micro-conductivity 
%    Sh1_diff_gain = 1.01;   % Shear Channel 1 differentiator gain ~1.0
%    Sh2_diff_gain = 0.97;   % Shear Channel 2 differentiator gain ~1.0
%    ADC_FS = 4.096;         % Full-scale voltage of A-to-D converter
%    ADC_bits = 16;          % Number of bits in A-to-D converter
%    % ************************************************************************
%
% Some variables are not completely converted to physical units. The thermistor
% and micro-conductivity signals are converted to physical units using the
% nominal sensitivity of these sensors by a linear transformation if you are
% using a data file prior to version 6.  A linear transformation is reversible.
% The values will be close to the true values but not as close as they
% potentially could be. The reason is that these sensors are not calibrated. We
% expect that they will be calibrated using in situ measurements and a
% comparison against the Sea-Bird SBE3F thermometer and the SBE4C conductivity
% cell. The variables are then saved to the mat-file and then reloaded. If the
% mat-file already exists, then the reading of the binary file and its
% conversion are by-passed.
%
% One of the figures shows the results of calibrating the thermistors but this
% result is not saved into the mat-file because the process is non-linear.
% However, it could be saved for future use by modifying the function. If the
% instrument carries micro-conductivity sensors, then calibration is also
% performed and shown in another window. Again, the results are not saved.
% One of the figures shows the frequency spectrum of acceleration, shear,
% temperature and micro-conductivity for the pressure range specified in the
% input arguments to this function. One can re-run this function for different
% pressure ranges to examine the shape of the spectra and to decipher the noise
% limits of the instrument. A wavenumber spectrum of the shear probe signals is
% shown in another figure. This section of the function does extensive
% computations and corrects the spectrum for the spatial resolution of the shear
% probes and for accelerometer-coherent vibrations of the profiler. It also,
% automatically finds the limits of spectral integration to compute the shear
% variance, which then is used to calculate the rate of dissipation of kinetic
% energy.
%
% By default, the function saves plots as figures (.fig) and PDF files (.pdf).
% The figures allow the plots to be reformatted at any time in the future.  The
% PDF files are high quality vector files which can be incorporated into most
% documents.  Example plots generated by quick_look follow:
%
% @image @images/quick_look1.pdf @Quick_look example for pitch and roll. 
% @Quick_look example where pitch and roll are derived from the horizontal 
% accelerometers and the fall-rate from the pressure transducer. If the 
% instrument has a magnetometer then the magnetic orientation is shown in 
% degrees divided by 100.
% 
% @image @images/quick_look2.pdf @Quick_look exmple plot of Sea-Bird data. 
% @Quick_look example plots of data from Sea-Bird sensors.  This figure is blank
% if the instrument does not carry Sea-Bird sensors.
% 
% @image @images/quick_look3.pdf @Quick_look example of polynomial regression 
% of the thermistor against the Sea-Bird thermometer.
% @Quick_look example showing the results of a polynomial regression of the 
% thermistors against the Sea-Bird thermometer. The user has to manually adjust 
% the value of the separation of the Sea-Bird thermometer and the thermistors by
% editing the function. This instrument did not have a second thermistor 
% installed. The polynomial coefficients indicate that the nominal thermistor 
% conversion to physical units had to be shifted by 0.306 degrees Celsius and 
% scaled by 0.97405 or its inverse, in order to agree with the Sea-Bird 
% thermometer.
%
% @image @images/quick_look4.pdf @Quick_look example showing vertical gradients
% of the micro-structure signals.
% @Quick_look example showing the vertical gradients of the micro-structure 
% signals. LP indicates low-pass filtering at 30 Hz. This instrument had only a 
% single thermometer and shear probe.
%
% @image @images/quick_look5.pdf @Quick_look example spectra plot.
% @Quick_look plotted spectra for the pressure range selected by the input 
% arguments to this function. All spectra are in physical units. Thermistor T2 
% and shear probe 2 were not installed, consequently, their spectra show the 
% noise level of the electronics.
%
% @image @images/quick_look6.pdf @Quick_look example spectra plot also showing
% the Nasmyth spectra.
% @Final quick_look plot where the thin (red and green) lines are the wavenumber
% spectra for the selected pressure range. The thick lines are the spectra after 
% correction for acceleration-coherent vibrations and the spatial averaging by 
% the shear probes. The dotted lines are a polynomial fit used to estimate the 
% upper limit of spectral integration for calculating the shear variance 
% (triangle). The solid black lines are the Nasmyth spectra for the estimated 
% rate of dissipation. The second shear probe was not installed, so its spectrum
% is mostly in nano-land below the limits of the figure.

% Version History:
% 2006-04-15 (RGL)
% 2006-05-24 (RGL) added conditional plotting for magnetometer
% 2008-09-27 (RGL) Modified by adding more variables into the top of
%     the program and no assumed sampling rate for spectra.
% 2009-01-16 (RGL) Modiied to produce wavenumber spectrum,
%    specfic depth range for spectra, file must now be specfied 
%    on input to function, and rationalization of the names of certain 
%    variables.
% 2009-08-07 (RGL) fixed nagging problem with the alinement of FP07 and 
%    SBE7 with Sea-Birds. Minimum pressure now correct.
% 2010-01-15 (AWS) support odas v6 and up
% 2010-10-04 (AWS) added tags for Doxygen
% 2012-05-11 (RGL) modified to use for loops for saving and conversion,
%       and check for existence of variables.
% 2012-07-19 (WID) removal of inifile_with_instring.m
% 2012-10-12, RGL some changes for Partrac measurements in Scotish Tidal
%       Channel.
% 2012-10-15, RGL modified to make slides for Scotish presentations.
%       Removed low-pass filtered profiles.
% 2012-10-31, RGL, numerous changes to pretty up the figures. Now we
%       same the dissipation structure for later use.
% 2012-11-09 (WID) made looking for input file case-indifferent.
% 2012-11-09 (WID) documentation update
% 2013-02-26 (WID) merge of changes made by Rolf.

function quick_look(fname, P_start, P_end)

% ******************************************************************************
% Play with the following variables to optimize the outputs.
% ******************************************************************************
P_min = 3.0;    % Minimum depth for profile start
W_min = 0.05;    % Minimum speed for profile start
Rise = 0;       % for rising profilers
fft_length = 1; % length in seconds
diss_length = 3*fft_length;
over_lap = diss_length/2;
HP_cut = 0.25;  % Frequency (Hz) of high-pass filter cut-off for shear probes.
LP_CUT = 25;    % Frequency (Hz) of low-pass filter for shear signals. 

% ******************************************************************************
% Required for processing legacy data files.  The following default values 
% should not be changed if processing version 6 data files.
% ******************************************************************************
T1_offset =  0;         % Adjust T1 by this amount
T2_offset =  0;         % Adjust T2 by this amount
sbt_offset = 0.35;      % distance between FP07 and SBE3 in meters 
                        %(usually 0.3 m or 2.5 m)
pump = 0;               % pump flag = 1 if SBT5 pump is installed, 0 otherwise.
sh1_sens      = 0.0759; % SN725, 2010-11-25
sh2_sens      = 0.0719; % shear probe SN735, 2010-11-25
P_diff_gain   = 19.9;   % Gain of presure pre-emphasis ~20.5
T1_diff_gain  = 0.99;   % T1 differentiator gain ~1.0
T2_diff_gain  = 0.99;   % T2 differentiator gain ~1.0
C1_diff_gain = 0.4;     % micro-conductivity 
Sh1_diff_gain = 1.01;   % Shear Channel 1 differentiator gain ~1.0
Sh2_diff_gain = 0.97;   % Shear Channel 2 differentiator gain ~1.0
ADC_FS = 4.096;         % Full-scale voltage of A-to-D converter
ADC_bits = 16;          % Number of bits in A-to-D converter
% ******************************************************************************


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

% __________________________________________________________________
%
% This is where we get some information about the data sampling

%header = header(:); header = reshape(header,64, length(header)/64)';
if (~we_have_old_mat_file), %We do not have a file that is already in physical
    
    number_of_records = size(header,1); % number of records in this file.
    length_of_record = (header(2,19) - header(2,18))/2; % number of samples in a record
    number_of_fast_columns = header(2,29);
    number_of_slow_columns = header(2,30);
    number_of_columns      = number_of_fast_columns + number_of_slow_columns;
    number_of_rows         = header(2,31);
    length_of_fast_vector  = (length_of_record / number_of_columns)*number_of_records; % for standard practice
    length_of_slow_vector  = length_of_fast_vector / number_of_rows; % for standard practice
    % double sampled channels and ground "fill" channels may have different
    % lengths
    
    t_slow = (0:length_of_slow_vector - 1)'/fs_slow; % make time vectors for plotting
    t_fast = (0:length_of_fast_vector - 1)'/fs_fast;% assumes that we always have variables P and sh1.
    % should be ok for all but weirdly sampled channels
    
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
            'U', 'V', 'W', 'V_Bat', 'Incl_X', 'Incl_Y', 'Incl_T','sh1', 'sh2'};
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
%         for ii = 1:length(convert_list),
%             item = convert_list{ii};
%             name = name_list{ii};
%             
%             if ~exist(item, 'var'), continue, end
%             tmp = convert_odas( eval(item), name, 'string', setupfilestr, header_version );
%             eval([name ' = tmp;']);
%         end
        for k = 1:length(convert_list)
            if exist(convert_list{k},'var')
                eval([convert_list{k} '=' ...
                    'convert_odas(' convert_list{k} ', '''  name_list{k} ''', ' ...
                    '''string''' ' , setupfilestr, header_version);']);
            end
        end
    else
        P_hres = deconvolve('pressure', P, P_dP, fs_slow, P_diff_gain);
        T1_hres = deconvolve('temperature', [], T1_dT1, fs_fast, T1_diff_gain);
        T2_hres = deconvolve('temperature', [], T2_dT2, fs_fast, T2_diff_gain);
        T1_hres = (ADC_FS / 2^ADC_bits)*T1_hres*15 + 15; % +/-1V for +/-15C
        T2_hres = (ADC_FS / 2^ADC_bits)*T2_hres*15 + 15; % +/-1V for +/-15C
        if exist('C1_dC1','var'),
            C1_hres =  despike(C1_dC1, 20, 1, fs_fast, 40);
            C1_hres = deconvolve('conductivity', [], C1_hres, fs_fast, C1_diff_gain);
            C1_hres = (ADC_FS / 2^ADC_bits)*C1_hres*40 + 40; % +/-1V for +/-30 mS/cm.
        end
        if exist('C2_dC2','var'),
            C2_hres =  despike(C2_dC2, 20, 1, fs_fast, 40);
            C2_hres = deconvolve('conductivity', [], C2_hres, fs_fast, C2_diff_gain);
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
        if exist('Az','var'), Az  = convert_odas (Az,'az'   , 'file', 'setup.txt');end

        if exist('Mx','var'), Mx = convert_odas (Mx,'Mx', 'file', 'setup.txt'); end
        if exist('My','var'), My = convert_odas (My,'My', 'file', 'setup.txt'); end
        if exist('Mz','var'), Mz = convert_odas (Mz,'Mz', 'file', 'setup.txt'); end

        if exist('Ux','var'),   Ux =   convert_odas (Ux,   'Ux',  'file', 'setup.txt'); end
        if exist('Uy','var'),   Uy =   convert_odas (Uy,   'Uy',  'file', 'setup.txt'); end
        if exist('Vbat','var'), Vbat = convert_odas (Vbat,'Vbat', 'file', 'setup.txt'); end
    end
    
    % Generate pressure and vertical speed vectors
    % 1.5 Hz seems to give the best smoothing for slow profiling
    W_slow = gradient(P_hres, 1/fs_slow);
    [b,a] =butter(4,1.5/fs_slow/2);
    W_slow = filtfilt(b,a,W_slow);
    P_fast = interp1(t_slow, P_hres,t_fast,'spline','extrap');
    W_fast = interp1(t_slow, W_slow,t_fast,'spline','extrap');
    profile_speed = W_fast; %
    
    if exist('sh1','var'), sh1 = sh1 ./ profile_speed.^2; end
    if exist('sh2','var'), sh2 = sh2 ./ profile_speed.^2; end
    
    % generate gradient signals
    gradient_list        = {'gradT1',  'gradT2',  'gradC1',  'gradC2'};
    gradient_source_list = {'T1_hres', 'T2_hres', 'C1_hres', 'C2_hres'};
    
    for k = 1:length(gradient_list)
        if exist(gradient_source_list{k},'var')
            eval([gradient_list{k} ' = ' gradient_source_list{k} ';']);
            eval([gradient_list{k} '(2:end) = diff(' gradient_source_list{k} ') *fs_fast;']);
            eval([gradient_list{k} '(1) = ' gradient_list{k} '(2);']);
            eval([gradient_list{k} ' = ' gradient_list{k} './ profile_speed;']);
        end
    end
    
    % Save the converted variables to the mat-file
    eval (['save ' N '.mat t_slow t_fast Year Month Day Hour Minute -append -v6'])
    
    var_list = {'P', 'P_hres', 'T1_hres', 'T2_hres', 'sbt', 'sbc', ...
        'Ax', 'Ay', 'Az', 'Rx', 'Ry', 'Rz', 'GC_T', 'APy', 'APz', ...
        'W_fast', 'W_slow', 'P_fast', 'sh1', 'sh2', 'gradT1', 'gradT2', ...
        'gradC1', 'C1_hres', 'gradC2', 'C2_hres', 'Mx', 'My', 'Mz', ...
        'U', 'V', 'W', 'V_Bat', 'V_bat', 'Incl_X', 'Incl_Y', 'Incl_T', 'PTV' };
    
    for k = 1: length(var_list)
        if exist(var_list{k},'var')
            flag = save_odas([N '.mat'], var_list{k},  eval(var_list{k}) );
            if (flag < 0), error(['Could not write ' var_list{k} ' to '  N '.mat']), end
        end
    end
    
    
    % 'end of conversion to physical units' for version 6
    
    % This is the end of the section for converting to physical units. This
    % section was skipped if the mat-file already existed. If the mat-file is
    % corrupted or has a bad conversion, then erase it and run this function
    % again
end


% _______________________________________________________________________
% Make figures
% Pitch, Roll Fall-Rate
n = find(P_fast> P_min & W_fast > W_min); % No up-casts, create index for profile
m = find(P_hres> P_min & W_slow > W_min);

if isempty(n)
    warning('No profile detected');
    return
end

Year_string   = num2str(Year);
Month_string  = sprintf('%02d',Month);
Day_string    = sprintf('%02d',Day);
Hour_string   = sprintf('%02d',Hour);
Minute_string = sprintf('%02d',Minute);

figure_file_name = ['QL_' N '_' num2str(round(P_start)) '_' num2str(round(P_end))] ;
fig_num = 0;
set(0,'Units','pixels') 
screen_size = get(0,'ScreenSize');
portrait_size  = [round(0.45*screen_size(3)) round(0.75*screen_size(4))];
landscape_size = [round(0.65*screen_size(3)) round(0.80*screen_size(4))];
figure_position = [20 40  portrait_size];

title_string = [ fix_underscore(N) '; '  Year_string '\_' Month_string '\_' Day_string ...
    ', ' Hour_string ':' Minute_string ' UTC'];

%
%---------------------------------------------------------------
fig_num = 0;
% Figure 01
fig_num = fig_num + 1;
figure (fig_num); clf
[b,a] =   butter(1,1/(fs_slow/2));% low-pass to show only the gravity signal
Incl_X_LP = filtfilt(b,a,Incl_X);
Incl_delay = 0.57; % Inclinometer delay in seconds due to internal averaging
%                    of 256 points and 482 samples per second.

subplot(1,2,1)
h = plot(Incl_X_LP(m), P_hres(m) - W_slow(m)*Incl_delay);grid
set(h(1), 'linewidth', 2,   'color','b')
set(gca, 'ydir', 'rev')
legend('Roll',1)
set(gca,'fontsize', 16)
ylabel('\it P \rm [dBar]','fontsize',16)
xlabel('[ ^{\circ} ]','fontsize',16)
title(title_string, 'fontsize', 15,'horizontalalignment','left')
Y_lim = get(gca,'ylim');
if (Y_lim(1) < -15), Y_lim(1) = -15; end
if (Y_lim(2) < +15), Y_lim(2) = +15; end
set(gca,'ylim',Y_lim)

subplot(1,2,2)
h = plot(W_slow(m), P_hres(m)); grid
set(h(1), 'linewidth', 2,   'color','b')
set(gca, 'ydir', 'rev')
legend('W',1)
set(gca,'fontsize', 16)
set(gca,'yticklabel',[])
xlabel('[m s^{-1}]','fontsize',15)

%axis tight
orient portrait

figure_name = [N  '_Fig_' sprintf('%02d',fig_num)];
saveas (fig_num, figure_name, 'fig');
fig2pdf(fig_num, figure_name);

%
%_________________________________________________________
% Dumb bell accelerometers
% Figure 02
fig_num = fig_num + 1;
figure (fig_num); clf
    [b,a] =   butter(1,1/(fs_fast/2));% low-pass to show only the gravity signal
plot_data = [Ax Ay];
    
h = plot([plot_data(n,:) filtfilt(b,a,plot_data(n,:))], P_fast(n));grid
set(h(1), 'color','b')
set(h(2), 'color','r')
set(h(3), 'color', 'k','linewidth',1);
set(h(4), 'color','k','linewidth',1)
set(gca, 'ydir', 'rev','fontsize',16)
xlabel('[counts]','fontsize',16)
ylabel(' \itP\rm [dbar]','fontsize',16)
hh=legend('\itA_x', '\itA_y','location','NorthEastOutside');
set(hh,'fontsize', 14)
title( title_string, 'fontsize',15)
axis tight
orient portrait

figure_name = [N  '_Fig_' sprintf('%02d',fig_num)];
saveas (fig_num, figure_name, 'fig');
fig2pdf(fig_num, figure_name);
%
%----------------------------------------------------------
% Figure 03
% Thermistors

mm = find((P_hres> P_min) & W_slow > W_min);
nn = find(P_fast >= P_hres(mm(1)));

fig_num = fig_num +1;
figure (fig_num); clf
legend_string = '\itT\rm_1';
plot_data = T1_hres + T1_offset;
if exist('T2','var')
    plot_data = [ plot_data T2_hres+T2_offset];
    legend_string = {legend_string, '\itT\rm_2'};
end
h = plot(plot_data(n,:), P_fast(n));grid
legend(legend_string,  'location','NorthEastOutside')

set(gca, 'ydir', 'rev', 'fontsize', 16)
ylabel('\it P \rm [dBar]','fontsize',16)
xlabel('[  ^{\circ}C ]','fontsize',16)
title(title_string, 'fontsize',15)
orient landscape

figure_name = [N  '_Fig_' sprintf('%02d',fig_num)];
saveas (fig_num, figure_name, 'fig')
fig2pdf(fig_num, figure_name);

%
%--------------------------------------------------------------------------
% Figure 04
% Temperature gradient, microconductivity gradient, and shear
fig_num = fig_num + 1;
figure(fig_num); clf

sh1(n) = despike(sh1(n), 7, 1, fs_fast, 20);
sh2(n) = despike(sh2(n), 7, 1, fs_fast, 20);
Ax(n)  =  despike(Ax(n),20, 1, fs_fast, 20);
Ay(n)  =  despike(Ay(n),20, 1, fs_fast, 20);

% Apply gentle high-pass filter
[b,a] = butter(1,HP_cut/(fs_fast/2)); % HP at ~0.25 cpm
sh1_HP = sh1;
sh2_HP = sh2;
Ax_HP  = Ax;
Ay_HP  = Ay;

sh1_HP(n) = filtfilt(b,a,sh1(n)) - sh1(n);
sh2_HP(n) = filtfilt(b,a,sh2(n)) - sh2(n);
Ax_HP(n)  = filtfilt(b,a, Ax(n)) -  Ax(n);
Ay_HP(n)  = filtfilt(b,a, Ay(n)) -  Ay(n);

% [b,a] = butter(4, LP_CUT/(fs_fast/2));% for profile plotting only
% sh1_BP    = sh1_HP;
% sh2_BP    = sh2_HP;
% sh1_BP(n) = filtfilt(b,a,sh1_HP(n));
% sh2_BP(n) = filtfilt(b,a,sh2_HP(n));

plot_gradient_wish_list = {'sh1_HP', 'sh2_HP',  ...
    'gradT1', 'gradT2'};
legend_gradient_wish_list = {'sh_1 HP', 'sh_2 HP',  ...
    '\partialT_1/\partialz', '\partialT_2/\partialz'};
plot_gradient_list = [];

plot_scalar_wish_list = {'T1_hres', 'T2_hres'};
legend_scalar_wish_list = {'T_1','T_2'};
plot_scalar_list = [];

index_gradient = 0;
for k = 1:length(plot_gradient_wish_list)
    if exist(plot_gradient_wish_list{k},'var')
        index_gradient = index_gradient + 1;
        plot_gradient_list{index_gradient} = plot_gradient_wish_list{k};
        legend_list{index_gradient} = legend_gradient_wish_list{k};
    end
end
% assume all are fast channels with same sampling rate
X_gradient_data = zeros(length(n),index_gradient);
x_gradient_offset = -10;
LL = index_gradient; % the number of gradient channels

index_scalar = 0;
for k = 1:length(plot_scalar_wish_list)
    if exist(plot_scalar_wish_list{k},'var')
        index_scalar = index_scalar + 1;
        plot_scalar_list{index_scalar} = plot_scalar_wish_list{k};
        legend_list{LL + index_scalar} = legend_scalar_wish_list{k};
    end
end
% assume all are fast channels with same sampling rate
X_scalar_data = zeros(length(n),index_scalar);
x_scalar_offset = 0;

for k = 1: index_gradient
   eval(['X_gradient_data(:,k) = '  plot_gradient_list{k} '(n);'])
   X_gradient_data(:,k) = X_gradient_data(:,k) + x_gradient_offset;
   x_gradient_offset = x_gradient_offset + 10;
end

for k = 1: index_scalar
   eval(['X_scalar_data(:,k) = '  plot_scalar_list{k} '(n);'])
   X_scalar_data(:,k) = X_scalar_data(:,k) + x_scalar_offset;
   x_scalar_offset = x_scalar_offset + 0;
end

if (isempty(X_gradient_data) && isempty(X_scalar_data))
    error ('no vertical profile')
end

X_data = [];
if ~isempty(X_gradient_data)
    X_data = X_gradient_data;
end
if ~isempty(X_scalar_data)
    X_data = [X_data X_scalar_data];
end

plot(X_data, P_fast(n) ); grid
set(gca,'ydir','rev','fontsize',12)
legend(legend_list, 'location','eastoutside')
ylabel('\itP\rm [dBar]','fontsize',14)
title(title_string, 'fontsize',15)

figure_name = [N  '_Fig_' sprintf('%02d',fig_num)];
saveas (fig_num, figure_name, 'fig');
fig2pdf(fig_num, figure_name);

%
% Figure 05
% ----------------------------------------------------------------
% Frequency spectra
fig_num = fig_num + 1;
figure(fig_num); clf

fft_num = round(fft_length*fs_fast);% fft length in units of samples

nn = find((P_fast(n) > P_start) & (P_fast(n) <= P_end) & (W_fast(n) > W_min));
if isempty(nn)
    error(['The pressure range of ' num2str(P_start) ' to ' num2str(P_end) ...
        ' was not found in this profile'])
end
range = n(nn);

P_start = round(P_fast(range(1))); % So the labels give the actual range of data used
P_end   = round(P_fast(range(end)));

info.fft_length    = fft_num;
info.diss_length   = length(range);
info.over_lap      = info.diss_length;
info.fs            = fs_fast;
info.speed         = W_fast(range);
info.T             = T1_hres(range);
info.P             = P_fast(range);
info.fit_2_Nasmyth = 1; % Force a fit to the Nasmyth spectrum
%info.Tolerance   = 0.2;

A = [Ax(range) Ay(range)];

diss = get_diss([sh1_HP(range) sh2_HP(range)], A, info);
W_mean = diss.speed;
F = diss.F;
K = diss.K;
% get_diss returns wavenumber spectra. Convert to frequency spectra by
% dividing by the mean speed. The spectra are corrected for the wavenumber
% response of the shear probe.
P_sh1       = squeeze(diss.sh(1,1,1,:))       / W_mean;
P_sh2       = squeeze(diss.sh(1,2,2,:))       / W_mean;
P_sh1_clean = squeeze(diss.sh_clean(1,1,1,:)) / W_mean;
P_sh2_clean = squeeze(diss.sh_clean(1,2,2,:)) / W_mean;


P_Ax =squeeze(diss.AA(1,1,1,:));
P_Ay =squeeze(diss.AA(1,2,2,:));

P_Ax = 1e-7*P_Ax; % we have a piezo-accelerometer. Scale to bring it on axis
P_Ay = 1e-7*P_Ay;
if exist('gradT1','var'), [P_gradT1, F] = psd_rolf(gradT1(range), fft_num, fs_fast);end
if exist('gradT2','var'), [P_gradT2, F] = psd_rolf(gradT2(range), fft_num, fs_fast);end

epsilon = 1e-10*[1 10 100 1e3 1e4]; %
nu = visc35(mean(T1_hres(range)));
[Pn, kn]=nasmyth(epsilon, nu);
fn = kn*mean(abs(W_fast(range))); Pn = Pn/mean(abs(W_fast(range)));

plot_wish_list = {'P_Ax', 'P_Ay', 'P_Az', 'P_sh1', 'P_sh2', ...
    'P_gradT1', 'P_gradT2'};
legend_wish_list = {'A_x', 'A_y', 'A_z', 'sh_1', 'sh_2', ...
    '\partialT_1/\partialz', '\partialT_2/\partialz'};

index = 0;
plot_list = [];
legend_list = [];
for k = 1:length(plot_wish_list)
    if exist(plot_wish_list{k},'var')
        index = index + 1;
        plot_list{index} = plot_wish_list{k};
        legend_list{index} = legend_wish_list{k};
    end
end
legend_list{index+1} = 'Nasmyth';

num_vectors = index;
Y_data = zeros(size(F,1),num_vectors);
for k = 1: num_vectors
   eval(['Y_data(:,k) = '  plot_list{k} ';'])
end

h = loglog(F, Y_data, fn, Pn, 'k');grid
hh=legend(legend_list,'location','NorthEastOutside');

for index = 1: size(Y_data,2)
    set(h(index), 'linewidth',1.5)
end

set(h(1),'linewidth',3,'color','b'); % Ax
set(h(2),'linewidth',3,'color',[0 0.5 0]); % Ay
set(h(3),'color','r'); % Ax
set(h(4),'color',[0 0.5 0]); % Ay

if exist('Az','var'), set(h(3),'linewidth',3 ,'color','3');end % Az

f_upper = max(F);

y_limit = [1e-10 10]; % limits for y-axis
x_limit = [1./fft_length f_upper];% limits for x-axis
set(gca, 'ylim', y_limit, 'xlim',x_limit)

for index = 1: length (epsilon)
    f_location = find(fn(:,index) > x_limit(1));
    if isempty(f_location)
        f_location = 1;
    else
        f_location = f_location(1);
    end
    h=text(fn(f_location,index), Pn(f_location,index), ...
        ['10^{' num2str(log10(epsilon(index)),2) '}'], ...
        'fontsize', 16, 'fontweight','bold');
end
ylabel('[Variance Hz^{-1}]','fontsize',14')
xlabel('\it f \rm [Hz]','fontsize',14)
new_title_string = title_string;

% new_title_string{size(new_title_string,2)+1} = ...
%     [ num2str(P_start) ' < P < ' num2str(P_end) ' m, f_{HP} = ' num2str(HP_cut) ' Hz'];
new_title_string = {title_string , ' ' , ...
    ['P =   ' num2str(P_start) ' - ' num2str(P_end) ' m' , ...
    ' ; W=' num2str(mean(W_fast(range)),2) ' m s^{-1}' , ...
    ' ; f_{HP} = ' num2str(HP_cut) ' Hz']};
title(new_title_string,'fontsize',15)
orient landscape

figure_name = [N '_' sprintf('%02d',P_start) '_' sprintf('%02d',P_end) ...
               '_Fig_' sprintf('%02d',fig_num)];
saveas (fig_num, figure_name);
fig2pdf(fig_num, figure_name);

%
% Figure 06
%_____________________________
% Wavenumber Spectra
fig_num = fig_num + 1;
figure(fig_num); clf

% Convert frequency spectra to wavenumber spectra
P_sh1_clean = P_sh1_clean * W_mean; P_sh1 = P_sh1 * W_mean;
P_sh2_clean = P_sh2_clean * W_mean; P_sh2 = P_sh2 * W_mean;

F_0 = 25*W_mean; % frequency response of thermistor aft Vachon and Lueck
P_gradT1 = P_gradT1 .* (1 + (F/F_0).^2); % frequency response correction
P_gradT1 = P_gradT1 .* W_mean; % conversion to wavenumbers
P_gradT2 = P_gradT2 .* (1 + (F/F_0).^2); % frequency response correction
P_gradT2 = P_gradT2 .* W_mean; % conversion to wavenumbers

% find index to the limit of integration.

K_max = diss.K_max;
K_max_index_1 = find (K == K_max(1));
K_max_index_2 = find (K == K_max(2));

nu = diss.nu;
e1 = diss.e(1);
phi1 = squeeze(diss.Nasmyth_spec(1,1,:));
phi2 = squeeze(diss.Nasmyth_spec(1,2,:));
%nasmyth(e1,nu,512);
e2 = diss.e(2);
%[phi2,k2] = nasmyth(e2,nu,512);

h = loglog(K, phi1, '--k', K, [P_sh1_clean P_sh1], ...
    K, phi2, 'k', K, [P_sh2_clean P_sh2], ...
    K_max(1), P_sh1_clean(K_max_index_1), '^', ...
    K_max(2), P_sh2_clean(K_max_index_2), '^', ...
    K, [P_gradT1 P_gradT2]) ;grid

set(gca,'xlim',[0.5 200], 'ylim', [1e-8 1e1])
set(h(1),'linewidth',1.5,'color',[0 0 0]) % Nasmith 1
set(h(2),'linewidth',4.0,'color',[1 0 0]) % sh1 clean
set(h(3),'linewidth',1.0,'color',[1 0 0]) % sh1 original
set(h(4),'linewidth',1.5,'color',[0 0 0]) % Nasmyth 2
set(h(5),'linewidth',4.0,'color',[0 0.5 0]) % sh2 slean
set(h(6),'linewidth',1.0,'color',[0 0.5 0]) % sh2 original
set(h(7), 'MarkerFaceColor',[1 0 0],  'MarkerSize',18) % sh1 limit
set(h(8), 'MarkerFaceColor',[0 0.5 0],'MarkerSize',18) % sh2 limit
set(h(9), 'linewidth',1.5,'color', [0.75 0.75 0]); % grad_T1
set(h(10),'linewidth',1.5,'color', [0.75 0 0.75]); % grad_T2

xlabel('\itk \rm [cpm]','fontsize',14)
ylabel('\Phi (\itk\rm)  [s^{-2} cpm^{-1}]','fontsize',14)
e1_string = ['\epsilon_1=' make_scientific(e1,2) 'W kg^{-1}'];
e2_string = ['\epsilon_2=' make_scientific(e2,2) 'W kg^{-1}'];

hh=legend(e1_string, ...
    '\partial u_1/\partial z_{clean} ','\partial u_1/\partial z ', ...
    e2_string, ...
    '\partial u_2/\partial z_{clean}','\partial u_2/\partial z', ...
    ['k_{max} u_1 ' num2str(round(K_max(1)))], ...
    ['k_{max} u_2 ' num2str(round(K_max(2)))], ...
    '\partial T_1/\partial z ', ...
    '\partial T_2/\partial z ', ...
    'orientation','vertical','location','NorthEastOutside');

set(gca, 'XMinorGrid', 'off','YMinorGrid', 'off','fontsize',16)
set(hh,'fontsize',12)
orient landscape
new_title_string = {title_string , ' ' , ...
    ['P =   ' num2str(P_start) ' - ' num2str(P_end) ' m' , ...
    ' ; W=' num2str(mean(W_fast(range)),2) ' m s^{-1}' , ...
    ' ; f_{HP} = ' num2str(HP_cut) ' Hz']};

title(new_title_string,'fontsize',15)

figure_name = [N '_' sprintf('%02d',P_start) '_' sprintf('%02d',P_end) ...
               '_Fig_' sprintf('%02d',fig_num)];
saveas( fig_num, figure_name, 'fig')
fig2pdf(fig_num, figure_name);


%
% Figure 07 
%--------------------------------------------------------------------------
% Velocity profiles
fig_num = fig_num + 1;
figure(fig_num); clf

f_int = 0.25; % in Hz
[b,a] = butter(1,f_int/(fs_fast/2));

u1 = filter(b,a,detrend(sh1(n)));
u2 = filter(b,a,detrend(sh2(n)));
u1 = u1/(2*pi*f_int); % now properly scaled velocity
u2 = u2/(2*pi*f_int);

u2_offset = max(u1 - u2);
u2_offset = ceil(10*u2_offset)/10;
if (u2_offset < 0.1), u2_offset = 0.1; end

h=plot([u1 (u2+u2_offset)], P_fast(n) ); grid
set(gca,'ydir','rev','fontsize',12)
set(h(1),'linewidth',1,'MarkerSize',18,'color',[1 0   0])
set(h(2),'linewidth',1,'MarkerSize',18,'color',[0 0.5 0])
ylabel('\itP\rm [dBar]','fontsize',14)
xlabel('[ m s^{-1}]', 'fontsize',14)
title(title_string, 'fontsize',15)
legend('\itu\rm_1', ['\itu\rm_2 + ' num2str(u2_offset,2)], 'location','northeastoutside')

figure_name = [N  '_Fig_' sprintf('%02d',fig_num)];
saveas (fig_num, figure_name, 'fig')
fig2pdf(fig_num, figure_name);

%
% Figure 08
% Now we calculate a dissipation profile

%n = start_index_fast : end_index_fast; n = n';
%m = start_index_slow : end_index_slow; m = m';

fig_num = fig_num + 1;
figure(fig_num); clf

info.fft_length  = fft_num;
info.diss_length = round(diss_length*fs_fast);
info.over_lap    = round(over_lap   *fs_fast);
info.fs          = fs_fast;
info.speed       = W_fast(n);
info.T           = T1_hres(n);
info.P           = P_fast(n);
info.fit_order   = 4;
info.Tolerance   = 0.3;

diss = get_diss([sh1_HP(n) sh2_HP(n)], [Ax(n) Ay(n)], info);

e1 = diss.e(:,1);
e2 = diss.e(:,2);
P_e = diss.P;

new_name = ['diss_' N];
eval([new_name ' = diss;'])
if ~exist('dissipation.mat', 'file'),
    save('dissipation', new_name, '-v6' );
else
    save('dissipation', new_name, '-v6','-append')
end

h=semilogx([e1 e2],P_e,'.-');grid
% junk = get(gca, 'ylim');
% set(gca,'ylim',[7 junk(2)])

set(h(1),'linewidth',1,'MarkerSize',18,'color',[1 0   0])
set(h(2),'linewidth',1,'MarkerSize',18,'color',[0 0.5 0])

set(gca,'ydir','rev','fontsize',14)
ylabel ('\it P \rm [dBar]','fontsize',14)
xlabel ('\epsilon [W kg^{-1}]','fontsize',14)
legend('\epsilon_1','\epsilon_2','location','NorthEastOutside');

orient landscape
title(title_string,'fontsize',15)

figure_name = [N  '_Fig_' sprintf('%02d',fig_num)];
saveas( fig_num, figure_name, 'fig')
fig2pdf(fig_num, figure_name);

fig_num = fig_num + 1;
figure(fig_num)
plot(diss.K_max, diss.P);grid
set(gca,'ydir','rev')


