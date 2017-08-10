function mooring_conditions(whichmooring)
% function mooring_conditions(whichmooring)
%
% where "whichmooring" is wether 'riki', 'n080' or 'both'
%    
% This function was mooring_quickdisplay.m at the beginning, but I
% change it since I quite upgraded it to display conditions of the
% mooring (wind, water temp., etc.). See svn revision70 for mooring_quickdisplay

% Which Mooring?
riki = strcmp(whichmooring, 'riki');  
n080 = strcmp(whichmooring, 'n080');  
both = strcmp(whichmooring, 'both');

if  riki+n080+both > 1
    disp(['Error with input! The only parameter must be ''riki'' or ' ...
          '''n080'' or ''both'''])
end
    
    
% ------------------------------ Wind --------------------------- %
disp('winds...')
data = load(['~/WINDEX/data_processing/sept_2011_mission/' ...
             'Mouillages/wind_tides/wind_mooring_6hF.dat']);

time_wind = data(:,1);
wind = data(:,2);
winddir = data(:,3);
clear data
% --------------------------------------------------------------- %



% ----------------------------- Tides ---------------------------- %
disp('tides...')
data = load(['~/WINDEX/data_processing/sept_2011_mission/' ...
             'Mouillages/wind_tides/tide_moorings.dat']);

time_tide = datenum(data(:,1), data(:,2), data(:,3), data(:,4), data(:,5), 0);
tide = data(:,6);
clear data
% --------------------------------------------------------------- %


% -------------------------- Temperature ------------------------ %
disp('temperature...')
if riki == 1
    load('~/WINDEX/data_processing/sept_2011_mission/Mouillages/M_RIKI/Tfield_M_RIKI.mat')
    Tgrid_riki = Tgrid;
elseif n080 == 1
    load('~/WINDEX/data_processing/sept_2011_mission/Mouillages/M_N080/Tfield_M_N080.mat')
    Tgrid_n080 = Tgrid;
else
    riki=1;
    load('~/WINDEX/data_processing/sept_2011_mission/Mouillages/M_RIKI/Tfield_M_RIKI.mat')
    Tgrid_riki = Tgrid;
    load('~/WINDEX/data_processing/sept_2011_mission/Mouillages/M_N080/Tfield_M_N080.mat')
    Tgrid_n080 = Tgrid;
end
% ---------------------------------------------------------------- %

% plot
%
timelim = [datenum(2011, 09, 20 ,0, 0 ,0) datenum(2011, 10, 12, 12, 0, 0)];
%timelim = [datenum(2011, 09, 27 ,0, 0 ,0) datenum(2011, 09, 30, 0, 0, 0)];
plot_mooring % for MN080
%plot_mooring2 for MRIKI



% $$$ % save figure
% $$$ print('-dpng', '-r300','mooring_conditions_both.png')
% $$$ set(gcf, 'renderer', 'painters')
% $$$ print('-depsc', 'mooring_conditions_both.eps')
% $$$ 
% $$$ % save figure
% $$$ print('-dpng', '-r300','mooring_conditions_N080_27-30.png')
% $$$ set(gcf, 'renderer', 'painters')
% $$$ print('-depsc', 'mooring_conditions_N080_27-30.eps')
% $$$ 

keyboard