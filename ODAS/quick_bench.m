%% quick_bench
% Quick evaluation of a data file collected while the instrument is on a bench.
%%
% <latex>\index{Functions!quick\_bench}</latex>
%
%%% Syntax
%   quick_bench( 'dataFileName', 'serialNumber' )
%
% * [dataFileName] String. The name of the file to be processed (extension
%        optional). 
% * [serialNumber] Serial number of the instrument as a string, or any
%        other useful information. This string is placed into a line of the
%        title of each figure. 
%
% * []
% * [empty] No return parameters but this function produces two figures.
%
%%% Description
%
% This function generates two figures from data collected with a RSI
% instrument, that you wish to test. The instrument should be on a bench, or
% just standing in a laboratory. Dummy probes should be installed when
% collecting data. Data should be collected for a few minutes. This
% function processes the data to produce time-series and spectra of some of
% the channels in the instrument. The resulting graphs allow the user to
% determine if an instrument is working correctly.  
%
% The graphs are primarily used to detect excessive noise within an instrument.
% This helps identify problems such as corroded connections and other faults
% which would otherwise go unseen. For example, you can compare the spectra
% produced by $\texttt{quick\_bench}$ against the noise spectra in the
% calibration report for your instrument.
%
% For real profiles taken in the ocean, use quick_look.m to verify that
% your instrument is working correctly.
%
%%% Examples:
%
%    >> quick_bench( 'data_001.p', '43' )
%
% Plot the data in file $\texttt{data\_001.p}$ collected with the instrument
% that has the serial number 43.  The serial number is not required, but
% the string will be added to the title of the figures.
%
% @image @images/quick_bench_1 @The time-series output from the
% $\texttt{quick\_bench}$ function. @The time-series output from the
% $\texttt{quick\_bench}$ function. Quick_bench tries to plot most of the
% variables in your raw data file. These include accelerometers, pressure
% signals, shear probes, thermistors, magnetometers, inclinometers,
% voltage-output oxygen sensors, micro-conductivity sensors, and JAC -T,
% -C, -Turbidity and -Chlorophyll sensors. For some signals, the function
% subtracts the mean and this is indicated in the legend on the right-hand
% side. Only the inclinometer signals are converted into physical units.
% All others remain in raw units of counts. 
%
% @image @images/quick_bench_2 @The spectra output from the
% $\texttt{quick\_bench}$ function.
% @Spectra of some of the signals shown in the time-series figure. The
% instrument should be well cushioned to minimize its vibrations. Even so,
% it is nearly impossible to suppress the output from the extremely
% sensitive accelerometers. AC power line frequency (50/60 Hz)
% contamination is also difficult to suppress and may show up as narrow
% spectral peaks. Dummy probes have been installed in place of the shear
% probes, thermistor and micro-conductivity sensors. Their spectra can be
% directly compared against those in your instrument calibration report to
% check if the noise level is close to that observed at RSI before your
% instrument was shipped.


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
% 2012-11-07 (WID) documentation update
% 2013-01-15 (WID) changed dpres to P_dP for v3.1 files
% 2013-06-10 (WID) updated to use modified read_odas
% 2013-09-10 (WID) modified to make use of texstr
% 2014-03-07 (WID) Fixed bug for when 2 uC probes are present
% 2014-04-03 (WID) Allow working with only 1 shear probe by faking sh2
% 2014-12-15 (RGL) updated to make function test the existance of every
%                  variable to be plotted and to generate a time vector
%                  independent of the existance of any particular channel.
% 2015-01-23 (RGL) Added a figure for JAC sensors, if they exist.
% 2015-02-18 (RGL) Mystery change... variable "n_fft_slow" created.
% 2015-03-02 (WID) Removed file check.  Use check within read_odas.
% 2015-10-30 (RGL) Corrected spectrum of P_dP using fs_slow. Documentation
%                  changes. Added saving the figures as both Matlab figures
%                  and as pdf-files. Set the gca fontsize explicitly to 16.

function quick_bench(fname, SN)

%
P_diff_gain   = 20.3; % Gain of presure pre-emphasis ~20.5, for legacy files
fft_length_in_seconds = 2;

% Set the default serial number
if nargin < 2, SN = '___'; end

% Let read_odas sort out the file name.  This is easier and facilitates the
% use of wilecards in the file name.
if nargin < 1
    [variable_list, data] = read_odas(); % convert to a mat-file
else
    [variable_list, data] = read_odas(fname); % convert to a mat-file
end

if isempty(variable_list)
    error(['No data found in data file: ' data.fullName]);
end
for i=fieldnames(data)', i=char(i); eval([i ' = data.' i ';']); end

% When wildcards are used fname is meaningless.  Must use name returned
% from read_odas.
[filepath, filename, fileext] = fileparts( fullPath );

% __________________________________________________________________
%
% This is where we get some information about the data sampling
Year  = header(2,4);
Month = header(2,5);
Day   = header(2,6);
Hour  = header(2,7);
Minute= header(2,8);

if exist('P_dP','var') && exist('P','var')
    P_hres = deconvolve('P_dP', P, P_dP, fs_slow, setupfilestr);
end
if exist('Incl_X','var'), Incl_X = convert_odas(Incl_X, 'Incl_X', setupfilestr);end
if exist('Incl_Y','var'), Incl_Y = convert_odas(Incl_Y, 'Incl_Y', setupfilestr);end
if exist('Incl_T','var'), Incl_T = convert_odas(Incl_T, 'Incl_T', setupfilestr);end


%_____________________________________________________________________
% Make figures

n_fft = round(fft_length_in_seconds*fs_fast); % length of fft for spectra in points
n_fft_slow = round(n_fft*fs_slow/fs_fast); % length of fft for spectra of slow channels

figure_file_name = ['QB_' SN '_' filename '_'] ;
figure_num = 0;
set(0,'Units','pixels')
screen_size = get(0,'ScreenSize');
landscape_size = round([0.90*screen_size(4) 0.65*screen_size(4)]);
figure_position = [20 150  landscape_size];

title_string = { texstr(sprintf('%s; %d_%02d_%02d, %02d:%02d UTC', ...
                                filename, Year, Month, Day, Hour, Minute)),
                 texstr(sprintf('SN_%s, Time Series', SN))};

figure_num = figure_num + 1;
figure_position(1) = figure_position(1)+20;
fig_handle = figure(figure_num);
clf(fig_handle)
%set(fig_handle,'paperorientation','landscape','Position',figure_position);clf
set(fig_handle,'Position',figure_position);

%
% Figure out how many subplots are required
windows = 0;
if exist('Ax','var') || exist('Ay','var') || exist('Az','var')
    windows = windows + 1;
end
if exist('sh1','var') || exist('sh2','var')
    windows = windows + 1;
end
if exist('P','var') || exist('P_dP','var') || exist('P_hres','var')
    windows = windows + 1;
end
if exist ('Mx','var') || exist ('My','var') || exist ('Mz','var')
    windows = windows + 1;
end
if exist('O2','var')
    windows = windows + 1;
end
if exist ('Incl_X','var') || exist ('Incl_Y','var') || exist('Incl_T','var')
    windows = windows + 1;
end
if (exist ('C1_dC1','var') || exist('C2_dC2','var'))
    windows = windows + 1;
end
if (exist('T1_dT1','var') || exist('T2_dT2','var'))
    windows = windows + 1;
end
window_index = 0;

spectra_legend = {};
spectra_fast   = [];
spectra_slow   = [];
F_slow         = [];

%
if exist('Ax','var') || exist('Ay','var') || exist('Az','var')
    window_index = window_index + 1;
    h(window_index)=subplot(windows,1,window_index);
    plot_vector = [];
    legend_string = {};
    if exist('Ax','var')
        plot_vector = [plot_vector Ax];
        legend_string{end+1} =  'A_x';
        [junk, F]  = csd_odas(Ax, Ax, n_fft, fs_fast, [], n_fft/2,'linear');
        spectra_fast = [spectra_fast junk];
        spectra_legend{end+1} = 'Ax';
    end
    if exist('Ay','var')
        plot_vector = [plot_vector Ay];
        legend_string{end+1} =  'A_y';
        [junk, F]  = csd_odas(Ay, Ay, n_fft, fs_fast, [], n_fft/2,'linear');
        spectra_fast = [spectra_fast junk];
        spectra_legend{end+1} = 'Ay';
    end
    if exist('Az','var')
        plot_vector = [plot_vector Az];
        legend_string{end+1} =  'A_z';
        [junk, F]  = csd_odas(Az, Az, n_fft, fs_fast, [], n_fft/2,'linear');
        spectra_fast = [spectra_fast junk];
        spectra_legend{end+1} = 'Az';
    end
    plot(t_fast, plot_vector); 
    legend(legend_string,'location','eastoutside')
    ylabel ('[counts]')
    set(gca,'xticklabel',[])
    set(gca, 'fontsize', 16)
end
%
%____________
if (exist('C1_dC1','var') || exist('C2_dC2','var'))
    window_index = window_index + 1;
    h(window_index)=subplot(windows,1,window_index);
    plot_vector = [];
    legend_string = {};
    if exist('C1_dC1','var')
        mean_C1_dC1 = mean(C1_dC1);
        plot_vector = C1_dC1 - mean_C1_dC1;
        value = num2str(round(abs(mean_C1_dC1)));
        if mean_C1_dC1<0;
            C1_string = ['C1\_dC1+' value];
        else
            C1_string = ['C1\_dC1-' value];
        end
        legend_string{end+1} = C1_string;
        [junk, F]  = csd_odas(C1_dC1, C1_dC1, n_fft, fs_fast, [], n_fft/2,'linear');
        spectra_fast = [spectra_fast junk];
        spectra_legend{end+1} = 'C1';
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
        legend_string{end+1} =  C2_string;
        [junk, F]  = csd_odas(C2_dC2, C2_dC2, n_fft, fs_fast, [], n_fft/2,'linear');
        spectra_fast = [spectra_fast junk];
        spectra_legend{end+1} = 'C2';
    end

    plot(t_fast, plot_vector); 
    legend(legend_string, 'location','eastoutside')
    ylabel ('[counts]')
    set(gca,'xticklabel',[])
    set(gca, 'fontsize', 16)
end

%
if exist('Mx','var') || exist('My','var') || exist('Mz','var')
    window_index = window_index + 1;
    h(window_index)=subplot(windows,1,window_index);
    plot_vector = [];
    legend_string = {};
    if exist('Mx','var')
        plot_vector = [plot_vector Mx];
        legend_string{end+1} = 'M_x';
    end
    if exist('My','var')
        plot_vector = [plot_vector My];
        legend_string{end+1} = 'M_y';
    end
    if exist('Mz','var')
        plot_vector = [plot_vector Mz];
        legend_string{end+1} = 'M_z';
    end
    plot(t_slow, plot_vector); 
    legend(legend_string, 'location','eastoutside')
    ylabel ('[counts]')
    set(gca, 'fontsize', 16)
    set(gca,'xticklabel',[])
end

%
if exist('O2','var')
    window_index = window_index + 1;
    h(window_index)=subplot(windows,1,window_index);
    plot(t_slow, O2); 
    legend('O_2 [Hz]', 'location','eastoutside')
    set(gca,'xticklabel',[])
    set(gca, 'fontsize', 16)
end

%
if exist('Incl_X','var') || exist('Incl_Y','var') || exist('Incl_T','var')
    window_index = window_index + 1;
    h(window_index)=subplot(windows,1,window_index);
    plot_vector = [];
    legend_string = {};
    if exist('Incl_X','var')
        plot_vector = [plot_vector Incl_X];
        legend_string{end+1} = 'Incl\_X';
    end
    if exist('Incl_Y','var')
        plot_vector = [plot_vector Incl_Y];
        legend_string{end+1} = 'Incl\_Y';
    end
    if exist('Incl_T','var')
        plot_vector = [plot_vector Incl_T];
        legend_string{end+1} = 'Incl\_T';
    end
    plot(t_slow, plot_vector); 
    legend(legend_string, 'location','eastoutside')
    ylabel ('[degrees]')
    set(gca,'xticklabel',[])
    set(gca, 'fontsize', 16)
end

%
if (exist('T1_dT1','var') || exist('T2_dT2','var'))
    window_index = window_index + 1;
    h(window_index)=subplot(windows,1,window_index);
    plot_vector = [];
    legend_string = {};
    if exist('T1_dT1','var')
        mean_T1_dT1 = mean(T1_dT1);
        plot_vector = T1_dT1 - mean_T1_dT1;
        value = num2str(round(abs(mean_T1_dT1)));
        if mean_T1_dT1<0;
            T1_string = ['T1\_dT1+' value];
        else
            T1_string = ['T1\_dT1-' value];
        end
        legend_string{end+1} = T1_string;
        [junk, F]  = csd_odas(T1_dT1, T1_dT1, n_fft, fs_fast, [], n_fft/2,'linear');
        spectra_fast = [spectra_fast junk];
        spectra_legend{end+1} = 'T1';
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
        legend_string{end+1} = T2_string;
        [junk, F]  = csd_odas(T2_dT2, T2_dT2, n_fft, fs_fast, [], n_fft/2,'linear');
        spectra_fast = [spectra_fast junk];
        spectra_legend{end+1} = 'T2';
    end

    plot(t_fast, plot_vector); 
    legend(legend_string, 'location','eastoutside')
    ylabel ('[counts]')
    set(gca,'xticklabel',[])
    set(gca, 'fontsize', 16)
end

%
if (exist('sh1','var') || exist('sh2','var'))
    window_index = window_index + 1;
    h(window_index)=subplot(windows,1,window_index);
    plot_vector = [];
    legend_string = {};
    if exist('sh1','var')
        plot_vector = sh1;
        legend_string{end+1} = 'sh1';
        [junk, F]  = csd_odas(sh1, sh1, n_fft, fs_fast, [], n_fft/2,'linear');
        spectra_fast = [spectra_fast junk];
        spectra_legend{end+1} = 'sh1';
    end
    if exist('sh2','var')
        plot_vector = [plot_vector sh2];
        legend_string{end+1} = 'sh2';
        [junk, F]  = csd_odas(sh2, sh2, n_fft, fs_fast, [], n_fft/2,'linear');
        spectra_fast = [spectra_fast junk];
        spectra_legend{end+1} = 'Sh2';
    end

    plot(t_fast, plot_vector); 
    legend(legend_string, 'location','eastoutside')
    ylabel ('[counts]')
    set(gca,'xticklabel',[])
    set(gca, 'fontsize', 16)
end

%
if (exist('P','var') || exist('P_dP','var') || exist('P_hres','var'))
    window_index = window_index + 1;
    h(window_index)=subplot(windows,1,window_index);
    plot_vector = [];
    legend_string = {};
    if exist('P','var')
        plot_vector = P;
        legend_string{end+1} = 'P';
    end
    if exist('P_dP','var')
        plot_vector = [plot_vector P_dP];
        legend_string{end+1} = 'P\_dP';
        [junk, F_slow]  = ...
            csd_odas(P_dP, P_dP, n_fft_slow, fs_slow, [], n_fft_slow/2,'linear');
        spectra_slow = junk;
        spectra_legend{end+1} = 'P\_dP';
    end
    if exist('P_hres','var')
        plot_vector = [plot_vector P_hres];
        legend_string{end+1} = 'P\_hres';
    end

    plot(t_slow, plot_vector); 
    legend(legend_string, 'location','eastoutside')
    ylabel ('[counts]')
    set(gca,'xticklabel',[])
    set(gca, 'fontsize', 16)
end

subplot(windows,1,1)
title(title_string)

subplot(windows,1,windows)
xlabel('t [s]')
set(gca,'xticklabelmode','auto')

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
if exist('export_fig.m', 'file') == 2
    export_fig(gcf, fig_name, '-pdf', '-transparent') 
end
saveas(gcf, fig_name, 'fig')

%
% Plot if JAC C, T Chlorophyll or Turbidity sensors exist
windows = 0;
if exist('Turbidity','var')
    windows = windows + 1;
end
if exist('Chlorophyll','var')
    windows = windows + 1;
end
if exist('JAC_T','var')
    windows = windows + 1;
end
if exist('JAC_C','var')
    windows = windows + 1;
end
window_index = 0;
if windows > 0
    figure_position = [20 150  landscape_size];
    figure_num = figure_num + 1;
    figure_position(1) = figure_position(1)+50;
    fig_handle = figure(figure_num);
    clf(fig_handle)
    set(fig_handle,'Position',figure_position);
    
    if exist('Turbidity','var')
        window_index = window_index + 1;
        h(window_index)=subplot(windows,1,window_index);
        plot(t_fast, Turbidity); 
        legend('Turbidity', 'location','eastoutside')
        set(gca,'xticklabel',[])
        ylabel('[counts]')
        set(gca, 'fontsize', 16)
    end
    if exist('Chlorophyll','var')
        window_index = window_index + 1;
        h(window_index)=subplot(windows,1,window_index);
        plot(t_fast, Chlorophyll); 
        legend('Chlorophyll', 'location','eastoutside')
        set(gca,'xticklabel',[])
        ylabel('[counts]')
        set(gca, 'fontsize', 16)
    end
    if exist('JAC_T','var')
        window_index = window_index + 1;
        h(window_index)=subplot(windows,1,window_index);
        plot(t_slow, JAC_T); 
        legend('JAC\_T', 'location','eastoutside')
        set(gca,'xticklabel',[])
        ylabel('[counts]')
        set(gca, 'fontsize', 16)
    end
    if exist('JAC_C','var')
        window_index = window_index + 1;
        h(window_index)=subplot(windows,1,window_index);
        plot(t_slow, JAC_C); 
        legend('JAC\_C', 'location','eastoutside')
        set(gca,'xticklabel',[])
        ylabel('[counts]')
        set(gca, 'fontsize', 16)
    end
    
    subplot(windows,1,1)
    title(title_string)

    subplot(windows,1,windows)
    xlabel('t [s]')
    set(gca,'xticklabelmode','auto')

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
    %fig2pdf( fig_handle, fig_name );
    saveas(gcf, fig_name, 'fig')
end
%----------------------------------------------------------------
%spectra

title_string = { texstr(sprintf('%s; %d_%02d_%02d, %02d:%02d UTC', ...
                                filename, Year, Month, Day, Hour, Minute)),
                 texstr(sprintf('SN_%s, Spectra', SN))};

figure_num = figure_num + 1;
figure_position(1) = figure_position(1)+50;
fig_handle = figure(figure_num);
%set(fig_handle,'paperorientation','landscape','Position',figure_position);clf
set(fig_handle,'Position',figure_position);clf

h = loglog(F,spectra_fast, F_slow, spectra_slow);
legend(spectra_legend, 'location','eastoutside');
title(title_string)
for index = 1:length(h)
    set(h(index),'linewidth', 1.5)
end

xlim_lower = 0.9/fft_length_in_seconds;
xlim_upper = 1.1*fs_fast/2;
set(gca, 'ylim', [1e-4 1e2], 'xlim',[xlim_lower xlim_upper])
ylabel('[counts Hz^{-1}]')
xlabel('\it f \rm [Hz]')
set(gca, 'fontsize', 16)

fig_name = [figure_file_name '_Fig_' num2str(figure_num)];
if exist('export_fig.m', 'file') == 2
    export_fig(gcf, fig_name, '-pdf', '-transparent') 
end
saveas(gcf, fig_name, 'fig')
