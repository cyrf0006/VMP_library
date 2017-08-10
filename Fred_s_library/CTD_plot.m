function CTD_plot(fname)
% CTD_plot(fname)

load(fname)

dat = date(1,:);
h = str2num(time(1, 1:2));
m = str2num(time(1, 4:5)); 
s = str2num(time(1, 7:8)); 
dd = str2num(dat(1:2));
mm = str2num(dat(4:5));
yyyy = str2num(dat(7:10));

n = datenum(yyyy,mm,dd,h,m,s);


subplot(1,3,1)
plot(SBT, P)
set(gca, 'ydir', 'reverse')
%title('profil de temperature')
title({'profil de temperature'; datestr(n,0)})
ylabel('Profondeur (m)')
xlabel('T (degC)')

subplot(1,3,2)
plot(SBS, P)
set(gca, 'ydir', 'reverse')
title('profil de salinite')
title({'profil de salinite'; datestr(n,0)})
ylabel('Profondeur (m)')
xlabel('S (psu)')

DENS = sw_dens(SBS, SBT, P); %from CSIRO
SIG_T = DENS-1000;

subplot(1,3,3)
plot(SIG_T, P)
set(gca, 'ydir', 'reverse')
title('profil de densite')
title({'profil de densite'; datestr(n,0)})
ylabel('Profondeur (m)')
xlabel('SIGMA-T')

