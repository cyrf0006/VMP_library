clear


paperwidth = 14;%cm
paperheight = 5;%cm
month = 4:11;
TIK = [datenum(999, month,1, 0,0,0)]; %ticklabel for plot
TIK15 = [datenum(999, month,15, 0,0,0)]; %ticklabel for plot
T_LIM = [datenum(999,4,15, 0,0,0) datenum(999,11,15, 0,0,0)]; %XLIM
                                                              %for
                                                              %plot(apr15                                                   
                                                              %-15nov)

% Parameters for figure costumization
offset1 = 0.01; % offset between figure and colorbar
offset2 = 0.01; % offset between heigth of colorbar and heigth of figure
offset3 = 0.0005; % offset beteween figure and custom Xlabel 
offset4 = 0.02; % offset between Xlabel and OuterPosition at bottom

T_LIM2 = [datenum(999,4,1, 0,0,0) datenum(999,11,30, 0,0,0)]; %XLIM for plot
no_days = T_LIM(2)-T_LIM(1);
no_days2 = T_LIM2(2)-T_LIM2(1);
sub = 0; % 1 for subplot, 0 for independant plots


% Set time vector
for i = 1:length(month) 

    n(i) = datenum(999, month(i), 15, 0,0,0); % climatology, 15th of the month, year 999
    n = n';
       
end

fromKobs = load('T_core_mod.dat');
fromKcst = load('T_core_5e-5.dat');
obs = load('T_core_obs.dat');

%regression of observations
[PP,S]=polyfit(obs(:,1), obs(:,2), 1);
Tmin_o = PP(2)+obs(:,1).*PP(1);





figure(1)
clf
set(gcf,'PaperUnits','centimeters','PaperPosition',[10 10 paperwidth paperheight])


plot(obs(:,1), Tmin_o, 'k');
hold on

plot(fromKobs(:,1), fromKobs(:,2), '--k*');
plot(fromKcst(:,1), fromKcst(:,2), '--k+');

hold off

ylabel({'^{\circ}C'})
set(gca, 'xtick', obs(:,1))
datetick('x', 3, 'keepticks')

%set(gca, 'xticklabel', [], 'fontsize', 10)
%set(gca, 'xtick', TIK)
set(gca, 'fontsize', 10)
set(gca, 'xgrid', 'on','ygrid', 'on' )
set(gca, 'ylim', [-1 1.5])

% $$$ 
% $$$ %replace datetick for custumization
% $$$ for i=2:length(month)-1
% $$$     mm=datestr(n(i), 3);
% $$$     day_norm = (n(i)-T_LIM(1))/no_days;
% $$$     text(day_norm,-offset3,mm,'units','normalized', 'HorizontalAlignment', ...
% $$$          'center', 'verticalalignment', 'top', 'fontsize', 10)
% $$$ end


%save fig
set(gcf, 'renderer', 'painters')
print('-deps2', 'slope_compa.eps')
