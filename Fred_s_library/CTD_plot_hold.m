function CTD_plot_hold(fnames)

% CTD_plot_hold(fnames)
% usage ex: CTD_plot_hold(['profile001_07-21'; 'profile001_07-23'; 'profile001_07-27'; 'profile001_09-09'; 'profile001_09-10'; 'profile001_09-16'])


colorlegend = ['b'; 'r'; 'k'; 'y'; 'm'; 'g'];
file =1;
%keyboard

S = size(fnames);

while file<=S(1)

load(fnames(file, :))

dat = date(1,:);
h = str2num(time(1, 1:2));
m = str2num(time(1, 4:5)); 
s = str2num(time(1, 7:8)); 
dd = str2num(dat(1:2));
mm = str2num(dat(4:5));
yyyy = str2num(dat(7:10));

n = datenum(yyyy,mm,dd,h,m,s);

L(file, :) = datestr(n,1);

subplot(1,3,1)
plot(SBT, P, colorlegend(file))
set(gca, 'ydir', 'reverse')
axis([-1 10 0 180])
%title('Temperature')
%title({'profil de temperature'; datestr(n,1)}, 'FontSize', 15)
ylabel('Profondeur (m)', 'FontSize', 20)
xlabel('Temperature (^{\circ}C)', 'FontSize', 20)
set(gca, 'FontSize', 20)
hold on

subplot(1,3,2)
plot(SBS, P, colorlegend(file))
set(gca, 'ydir', 'reverse')
axis([24 35 0 180])
set(gca, 'yticklabel', [])
%title('Salinite')
%title({'profil de salinite'; datestr(n,1)}, 'FontSize', 20)
%ylabel('Profondeur (m)', 'FontSize', 20)
xlabel('Salinite (psu)', 'FontSize', 20)
set(gca, 'FontSize', 20)
hold on
DENS = sw_dens(SBS, SBT, P); %from CSIRO
SIG_T = DENS-1000;

subplot(1,3,3)
plot(SIG_T, P, colorlegend(file))
set(gca, 'ydir', 'reverse')
axis([18 28 0 180])
%xlim([18 28])
set(gca, 'yticklabel', [])
%title('profil de densite')
%title({'profil de densite'; datestr(n,1)}, 'FontSize', 20)
%ylabel('Profondeur (m)', 'FontSize', 20)
xlabel('\sigma_T', 'FontSize', 20)
set(gca, 'FontSize', 20)
hold on

file = file+1;

end

% subplot(1,3,1)
% hold off
% legend(L,'Location','SouthWest')

subplot(1,3,2)
hold off
legend(L,'Location','SouthWest')

% subplot(1,3,3)
% hold off
% legend(L,'Location','SouthWest')
