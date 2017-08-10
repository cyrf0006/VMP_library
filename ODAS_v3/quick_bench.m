%% quick_bench
% Quick evaluation of a data file collected while the instrument is on a bench.
%%
% <latex>\index{Type B!quick\_bench}</latex>
%
%%% Syntax
%   quick_bench( 'dataFileName', 'serialNumber' )
% 
% * [dataFileName] Name of the file to be processed (extension optional).
% * [serialNumber] Serial number of the instrument as a string.  Used in 
%                  the resulting graphs.
% * []
% * [empty] No return parameters but this function produces two figures.
%
%%% Description
%
% This function generates plots and figures from data collected from an RSI 
% instrument that is being tested.  While the instrument is on a bench, data
% is collected then processed using this function.  Dummy probes should be 
% installed when collecting data.  The resulting graphs allow the user to 
% determine if an instrument is working correctly.
%
% The graphs are primarily used to detect excessive noise within an instrument.
% This helps identify problems such as corroded connections and other faults 
% which would otherwise go unseen.
% 
% For real profiles taken in the ocean, use quick_look.m to verify an instrument
% is working correctly.
% 
%%% Examples:
%
%    >> quick_bench( 'data_001.p', '43' )
%
% Perform a plot of the data file 'data_001.p' for an instrument with serial 
% number 43.  The serial number is not required but when provided, will be 
% added to the resulting figures.
%
% @image @images/quick_bench_1 @The time-series output from the quick_bench 
% function. @Quick_bench tries to plot most of the variables in your raw data 
% file. These include accelerometers, pressure signals, shear probes, 
% thermistors, magnetometers, inclinometers, voltage-output oxygen sensors, 
% and micro-conductivity sensors. For some signals, the function subtracts the 
% mean and this is indicated in the legend on the right-hand side. Only the 
% inclinometer signals are converted into physical units. All others remain in 
% raw units of counts.
%
% @image @images/quick_bench_2 @quick_bench spectra output @Spectra of some of 
% the signals shown in the time-series figure. The instrument should be well 
% cushioned to isolate it from vibrations. Even so, it is nearly impossible to 
% suppress the output from the extremely sensitive accelerometers. Line 
% frequency (50/60 Hz) contamination is also difficult to suppress and may show 
% up as narrow spectral peaks. Dummy probes have been installed in place of the 
% shear probes, thermistor and micro-conductivity sensor. Their spectra can be 
% directly compared against those in your instrument calibration report to check
% if the noise level is close to that observed at RSI before your instrument was
% shipped.

% Revision History:
% 2007-01-01 (RGL) initial version
% 2007-11-05 (RGL) Added conditional plotting of magnetometer and oxygen sensor
% 2011-04-19 (AWS) support odas v6 and up, added tags for Doxygen, dynamically 
%                  build list of vectors and legends.
% 2011-05-04 (RGL) removed most scaling to leave data in raw units.
% 2011-07-26 (RGL) added optional serial number for figure titles
% 2011-09-01 (AWS) added documentation tags for matlab publishing
% 2012-04-11 (WID) changed inifile_with_instring calls to setupstr
% 2012-04-23 (WID) changed plotting functions to provide improved output for 
%                  latex input
% 2012-05-02 (RGL) Corrected the counting of windows for the first
%                  figure.
% 2012-11-07 WID documentation update
% 2013-01-15 WID changed dpres to P_dP for v3.1 files


function quick_bench(fname, SN)

if nargin<2, SN='___'; end
P_diff_gain   = 20.3; % Gain of presure pre-emphasis ~20.5
fft_length_in_seconds = 2;
% ____________________________________________________________
% File opening etc.
% Check file name, check for a .mat file, and open the .p file if no *.mat
% file
% _________________________________________________________________

% Ask the user for a file name if not provided as an input variable


dft_fname = get_latest_file;
if nargin<1
    fname = input(['Enter data file name (default: ' dft_fname '): '], 's');
    if isempty(fname), fname = dft_fname; end
end
% if nargin ~<1, then use the name passed by the function.

% Look for the .p and .mat files, deal with problems
[P,N,E,fname] = file_with_ext( fname,                               ...
                               {'' '.p' '.P'},                      ...
                               ['Unable to find file: ' fname] );

% Do not create a mat-file with a base name equal to the data file. Make a
% temporary mat-file and then wipe it out when finsihed.
if exist([N '_tmp.mat'],'file') 
    % Then the temp mat-file already exist, probably to do a previous 
    % failure and we will wipe it out
    display('file already exists... deleting');
    delete([N '_tmp.mat']);
end

[variable_list, tmpname] = read_odas(fname); % convert to a mat-file
if isempty(variable_list)
    error('Could not convert data into a mat-file');
end

disp(['Loading file  = ' tmpname]);
load(tmpname); % use the one we just converted
delete (tmpname);

%Do we need this??? Yes!!
if header_version >= 6
    cfg = setupstr(setupfilestr);
    tmp = setupstr(cfg, '', 'xmp');
    if ~isempty(tmp)
        error('please run quick_bench_XMP.m');
    end
end

% __________________________________________________________________
%
% This is where we get some information about the data sampling

    t_slow = (0:length(P)-1)'/fs_slow; % make time vectors for plotting
    t_fast = (0:length(sh1)-1)'/fs_fast;% assumes that we always have variables P and sh1.

    Year  = header(2,4);
    Month = header(2,5);
    Day   = header(2,6);
    Hour  = header(2,7);
    Minute= header(2,8);

    if (header_version >= 6) % new version files
        P_hres = deconvolve('P_dP', P, P_dP, fs_slow, setupfilestr, header_version);
        if exist('Incl_X','var'), Incl_X = convert_odas(Incl_X, 'Incl_X', 'string', setupfilestr, header_version);end 
        if exist('Incl_Y','var'), Incl_Y = convert_odas(Incl_Y, 'Incl_Y', 'string', setupfilestr, header_version);end 
        if exist('Incl_T','var'), Incl_T = convert_odas(Incl_T, 'Incl_T', 'string', setupfilestr, header_version);end 
    else % legacy odas

        P_hres = deconvolve('pressure', P, P_dP, fs_slow, P_diff_gain);
        if exist('Incl_X','var'), Incl_X = convert_odas(Incl_X, 'inclx', 'file', 'setup.txt');end 
        if exist('Incl_Y','var'), Incl_Y = convert_odas(Incl_Y, 'incly' ,'file', 'setup.txt');end 
        if exist('Incl_T','var'), Incl_T = convert_odas(Incl_T, 'inclt', 'file', 'setup.txt');end 
        
        %___________________________________________________
    end % legacy version
    
 
% _______________________________________________________________________
% Make figures

n_fft = round(fft_length_in_seconds*fs_fast); % length of fft for spectra in points

Year_string   = num2str(Year);
Month_string  = sprintf('%02d',Month);
Day_string    = sprintf('%02d',Day);
Hour_string   = sprintf('%02d',Hour);
Minute_string = sprintf('%02d',Minute);

figure_file_name = ['QB_' SN '_' N '_'] ;
figure_num = 0;
set(0,'Units','pixels') 
screen_size = get(0,'ScreenSize');
landscape_size = round([0.90*screen_size(4) 0.65*screen_size(4)]);
figure_position = [20 150  landscape_size];

title_string = {[ fix_underscore(N) '; '  Year_string '\_' Month_string '\_' Day_string ...
    ', ' Hour_string ':' Minute_string ' UTC'], ['SN\_' fix_underscore(SN) ', Time Series']};

figure_num = figure_num + 1;
figure_position(1) = figure_position(1)+20;
fig_handle = figure(figure_num);
clf(fig_handle)
%set(fig_handle,'paperorientation','landscape','Position',figure_position);clf
set(fig_handle,'Position',figure_position);

windows = 3;
if exist ('Mx','var')
    windows = windows + 1;
end
if exist('O2','var')
    windows = windows + 1;
end
if exist ('Incl_X','var')
    windows = windows + 1;
end
if (exist ('C1_dC1','var') || exist('C2_dC2','var'))
    windows = windows + 1;
end
if (exist('T1_dT1','var') || exist('T2_dT2','var'))
    windows = windows + 1;
end
window_index = 1;

h(window_index)=subplot(windows,1,window_index);
value = num2str(round(abs(mean(Ax))));
if mean(Ax)<0;
    Ax_string = ['A_x + ' value];
else
    Ax_string = ['A_x - ' value];
end
value = num2str(round(abs(mean(Ay))));
if mean(Ay)<0;
    Ay_string = ['A_y + ' value];
else
    Ay_string = ['A_y - ' value];
end
if exist('Az','var')
    value = num2str(round(abs(mean(Az))));
    if mean(Az)<0;
        Az_string = ['A_z + ' value];
    else
        Az_string = ['A_z - ' value];
    end
end

if ~exist('Az','var')
    plot(t_fast, [Ax-mean(Ax), Ay-mean(Ay) ]); grid
    legend(Ax_string, Ay_string, 'location','eastoutside')
    legend boxoff
else
    plot(t_fast, [Ax-mean(Ax), Ay-mean(Ay) Az-mean(Az) ]); grid
    legend(Ax_string, Ay_string, Az_string, 'location','eastoutside')
    legend boxoff
end

title(title_string, 'fontsize',12)
set(gca,'xticklabel',[])
set(gca,'fontsize',11)

%____________
if (exist('C1_dC1','var') || exist('C2_dC2','var'))
    window_index = window_index + 1;
    h(window_index)=subplot(windows,1,window_index);
    plot_vector = [];
    legend_string = [];
    if exist('C1_dC1','var')
        mean_C1_dC1 = mean(C1_dC1);
        plot_vector = C1_dC1 - mean_C1_dC1;
        value = num2str(round(abs(mean_C1_dC1)));
        if mean_C1_dC1<0;
            C1_string = ['C1\_dC1+' value];
        else
            C1_string = ['C1\_dC1-' value];
        end
        legend_string = [legend_string; C1_string];
        legend boxoff
    end
    if exist('C2_dC2','var')
        mean_C2_dC2 = mean(C2_dC2);
        plot_vector = [plot_vector C2_dC2 - mean_C2_dC2];
        value = num2str(round(abs(mean_C2_dC2)));
        if mean_C2_dC2<0;
            C2_string = ['C2\_dC2+' value];
        else
            C2_string = ['C2\_dC2-' value];
        end
        legend_string = {legend_string, C2_string};
        legend boxoff
    end
   
    plot(t_fast, plot_vector); grid
    legend(legend_string, 'location','eastoutside')
    legend boxoff
    set(gca,'xticklabel',[])
    set(gca,'fontsize',11)
end

%window_index = window_index + 1;
if exist('Mx','var')
    window_index = window_index + 1;
    h(window_index)=subplot(windows,1,window_index);
    plot(t_slow, [Mx, My, Mz]); grid
    legend('M_x', 'M_y', 'M_z', 'location','eastoutside') 
    legend boxoff
    set(gca,'xticklabel',[])
    set(gca,'fontsize',11)
end

if exist('O2','var')
    window_index = window_index + 1;
    h(window_index)=subplot(windows,1,window_index);
    plot(t_slow, O2); grid
    legend('O_2 [Hz]', 'location','eastoutside') 
    legend boxoff
    set(gca,'xticklabel',[])
    set(gca,'fontsize',11)
end

if exist('Incl_X','var')
    window_index = window_index + 1;
    h(window_index)=subplot(windows,1,window_index);
    plot(t_slow, [Incl_X Incl_Y ]); grid
    legend('Incl\_X', 'Incl\_Y',  'location','eastoutside') 
    legend boxoff
    ylabel('^{\circ}')
    set(gca,'xticklabel',[])
    set(gca,'fontsize',11)
end

if (exist('T1_dT1','var') || exist('T2_dT2','var'))
    window_index = window_index + 1;
    h(window_index)=subplot(windows,1,window_index);
    plot_vector = [];
    legend_string = [];
    if exist('T1_dT1','var')
        mean_T1_dT1 = mean(T1_dT1);
        plot_vector = T1_dT1 - mean_T1_dT1;
        value = num2str(round(abs(mean_T1_dT1)));
        if mean_T1_dT1<0;
            T1_string = ['T1\_dT1+' value];
        else
            T1_string = ['T1\_dT1-' value];
        end
        legend_string = [legend_string; T1_string];
    end
    if exist('T2_dT2','var')
        mean_T2_dT2 = mean(T2_dT2);
        plot_vector = [plot_vector T2_dT2 - mean_T2_dT2];
        value = num2str(round(abs(mean_T2_dT2)));
        if mean_T2_dT2<0;
            T2_string = ['T2\_dT2+' value];
        else
            T2_string = ['T2\_dT2-' value];
        end
        legend_string = {legend_string, T2_string};
    end
    
    plot(t_fast, plot_vector); grid
    legend(legend_string, 'location','eastoutside')
    legend boxoff
    set(gca,'xticklabel',[])
    set(gca,'fontsize',11)
end

window_index = window_index + 1;
h(window_index)=subplot(windows,1,window_index);
plot(t_fast, [sh1 sh2]);grid
set(gca,'xticklabel',[])
set(gca,'fontsize',11)
legend('Sh1','Sh2','location','eastoutside')
    legend boxoff
window_index = window_index + 1;

h(window_index)=subplot(windows,1,window_index);
plot(t_slow, [ P P_dP P_hres]);grid
legend('P','P\_dP','P\_hres','location','eastoutside')
    legend boxoff
xlabel('\it t \rm [s]','fontsize',12)
set(gca,'fontsize',11)

position = zeros(window_index,4);
for k = 1:window_index
    position(k,:) = get(h(k), 'position');
end

min_width = min(position(:,3));
position(:,3) = 0.95*min_width;% 100% causes problems on right edge

for k=1:window_index
    set(h(k), 'position',position(k,:))
end

fig_name = [figure_file_name '_Fig_' num2str(figure_num)];
fig2pdf( fig_handle, fig_name );
saveas(fig_handle, fig_name, 'fig')
%----------------------------------------------------------------
%spectra
figure(2)

title_string = {[ fix_underscore(N) '; '  Year_string '\_' Month_string '\_' Day_string ...
    ', ' Hour_string ':' Minute_string ' UTC'], ['SN\_' fix_underscore(SN) ', Spectra']};

figure_num = figure_num + 1;
figure_position(1) = figure_position(1)+50;
fig_handle = figure(figure_num);
%set(fig_handle,'paperorientation','landscape','Position',figure_position);clf
set(fig_handle,'Position',figure_position);clf

[P_Ax, F]  = csd_odas(Ax, Ax, n_fft, fs_fast, [], n_fft/2,'linear');
[P_Ay, F]  = csd_odas(Ay, Ay, n_fft, fs_fast, [], n_fft/2,'linear');
% start off with a given list of vectors and the corresponding legends
% and grow these two arrays depending on the presence of other variables
all_spectra = [P_Ax P_Ay ];
all_legends = {'A_x' 'A_y'};
if exist('Az','var'),
    [P_Az, F]  = csd_odas(Az, Az, n_fft, fs_fast, [], n_fft/2,'linear');
    all_spectra = [all_spectra P_Az];
    all_legends = [all_legends 'A_z' ];
end

[P_sh1, F] = csd_odas(sh1, sh1, n_fft, fs_fast, [], n_fft/2,'linear');
[P_sh2, F] = csd_odas(sh2, sh2, n_fft, fs_fast, [], n_fft/2,'linear');

all_spectra = [all_spectra P_sh1 P_sh2];
all_legends = [all_legends 'Sh1' 'Sh2' ];

if exist('T1_dT1','var')
    [P_T1, F] = csd_odas(T1_dT1, T1_dT1, n_fft, fs_fast, [], n_fft/2,'linear');
    all_spectra = [all_spectra P_T1];
    all_legends = [all_legends 'T1\_dT1'];
end
if exist('T2_dT2','var')
    [P_T2, F] = csd_odas(T2_dT2, T2_dT2, n_fft, fs_fast, [], n_fft/2,'linear');
    all_spectra = [all_spectra P_T2];
    all_legends = [all_legends 'T2\_dT2'];
end
if exist('C1_dC1','var')
    [P_C1, F] = csd_odas(C1_dC1, C1_dC1, n_fft, fs_fast, [], n_fft/2,'linear');
    all_spectra = [all_spectra P_C1];
    all_legends = [all_legends 'C1\_dC1'];
end
if exist('C2_dC2','var')
    [P_C2, F] = csd_odas(C2_dC1, C2_dC2, n_fft, fs_fast, [], n_fft/2,'linear');
    all_spectra = [all_spectra P_C2];
    all_legends = [all_legends 'C2\_dC2'];
end
if exist('P_dP','var')
    nn_fft = round(n_fft*fs_slow/fs_fast);
    [P_P, FF] = csd_odas(P_dP, P_dP, nn_fft, fs_slow, [], nn_fft/2,'linear');
    all_spectra_slow = [P_P];
    all_legends_slow = ['P\_dP'];
end

h = loglog(F,all_spectra,FF, all_spectra_slow);grid
legend([all_legends all_legends_slow], 'location','eastoutside');
title(title_string,'fontsize',12)
set(h(1),'linewidth', 1.5)
set(h(2),'linewidth', 1.5)
if exist('Az','var'), set(h(3),'linewidth', 1.5);end
set(gca,'fontsize',11)
xlim_lower = 0.9/fft_length_in_seconds;
xlim_upper = 1.1*fs_fast/2;
set(gca, 'ylim', [1e-4 1e2], 'xlim',[xlim_lower xlim_upper])
ylabel('[counts Hz^{-1}]')
xlabel('\it f \rm [Hz]')

fig_name = [figure_file_name '_Fig_' num2str(figure_num)];
fig2pdf( fig_handle, fig_name );
saveas(fig_handle, fig_name, 'fig')

