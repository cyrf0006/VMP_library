clear

%%%%%%%%%%%%%%%%%%
% ---- ADCP ----- %
%%%%%%%%%%%%%%%%%%

%adcp_file = 'adcp_2009-09-29.mat';
adcp_file = 'IML4_ADCP09.mat';

% load ADCP infos
load(adcp_file);
v = ADCP.north_vel;
u = ADCP.east_vel;
z = ADCP.config.ranges;
mtime = ADCP.mtime;

dz = z(2)-z(1);
dt = mtime(2)-mtime(1);
clear ADCP

% correct ADCP year in mtime
mm = str2num(datestr(mtime, 5));
yyyy = mm;
yyyy(:) = 2009;
dd = str2num(datestr(mtime, 7));
hour = datestr(mtime, 13);
hh = str2num(hour(:,1:2));
mi = str2num(hour(:,4:5));

time_adcp = datenum(yyyy, mm, dd, hh, mi, 0);


%%%%%%%%%%%%%%%%%%%%%%%%%
% ---- T-S-density ---- %
%%%%%%%%%%%%%%%%%%%%%%%%%

% load T-S infos
% $$$ load tclim_vector.dat;
% $$$ load Tclim_matrix.dat;
% $$$ load Sclim_matrix.dat;
% $$$ Tclim_matrix = Tclim_matrix';
% $$$ Sclim_matrix = Sclim_matrix';

% load T-S infos for 2009
load('datprofiles.dat')
load('tprofiles.dat')
load('sprofiles.dat')

yyyy = str2num(datestr(datprofiles, 10));
I=find(yyyy==2009);

% time vector
tvec = datprofiles(I);

% T-S matrix
Tmat = tprofiles(:,I);
Smat = sprofiles(:,I);
[X,Pmat] = meshgrid(1:21, 1:300);

% density matrix
Dmat = sw_dens(Smat, Tmat, Pmat);


%%%%%%%%%%%%%%%%%%%%%%%%%
% ---- GRIDDING ---- %
%%%%%%%%%%%%%%%%%%%%%%%%%

%make sure no double date in CTD
J = find(diff(tvec)==0);
if ~isempty(J)==1
    for i=1:length(J)
        tvec(J(i)+1) = [];
        Dmat(:,J(i)+1) = [];
        Tmat(:,J(i)+1) = [];
     end
end

%time vector
t1 = max([tvec(1); time_adcp(1)]);
%t2 = min([tvec(end); time_adcp(end)]);
t2 = datenum(2009, 10 ,28);
tv = t1:dt:t2;

%pressure vector
pv = z;

% depth interpolation of density to ADCP resolution
for i = 1:size(Dmat,2)
    
    J = find(~isnan(Dmat(:,i))==1); %ignore NaNs
    DM_0(:,i) = interp1(Pmat(J,1), Dmat(J,i), pv);
    TM_0(:,i) = interp1(Pmat(J,1), Tmat(J,i), pv);
    
end

% time interpolation of density to ADCP resolution
for i = 1:size(DM_0,1)
    
    J = find(~isnan(DM_0(i,:))==1); %ignore NaNs
3    DM(i,:) = interp1(tvec(J), DM_0(i,J), tv);
    TM(i,:) = interp1(tvec(J), TM_0(i,J), tv);
    
end

% restrict ADCP values to the same grid (press is OK, but time...)
I1 = nearestpoint(tv(1), time_adcp);
I2 = nearestpoint(tv(end), time_adcp);

u = u(:,I1:I2);
v = v(:,I1:I2);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ---- Richardson Number ---- %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%param
g = 9.81; %m2/s
rho_0 = 1015; %kg/m3

% shear prod.
dudz = diff(u,1)./dz;
dvdz = diff(v,1)./dz;
S2 = dudz.^2 + dvdz.^2;


%try to remove NaNs
for i = 1:size(S2,1)
    
    J = find(~isnan(S2(i,:))==1); %ignore NaNs
    S2(i,:) = interp1(tv(J), S2(i,J), tv);
    
end

%remove zeros in S2
I = find(S2==0);
S2(I)=NaN;


% buoy. prod.
N2 = g./rho_0./dz.*diff(DM,1);

% Richardson Number
Ri = N2./S2;
I = find(isnan(Ri)==1);
Ri(I)=1;
I = find(Ri>1);
Ri(I)=1;

figure(1)
clf
set(gcf,'PaperUnits','centimeters','PaperPosition',[10 10 30 20])
set(gcf, 'renderer', 'painters')
load Ri_colormap;
imagesc(tv, pv(1:end-1), Ri)
colormap(flipud(mycolormap));
caxis([0 1])
datetick('x')
colorbar
title('Richardson number at IML4')
xlabel('year 2009')
ylabel('depth (m)')

hold on
contour(tv, pv, DM, 'k')
hold off

print('-depsc2', 'Ri_IML4.eps') % save without Xlabel


figure(2)
clf
set(gcf,'PaperUnits','centimeters','PaperPosition',[10 10 30 20])
set(gcf, 'renderer', 'painters')
imagesc(tv, pv(1:end-1), S2)
colormap(jet);
caxis([0 0.01])
datetick('x')
colorbar
title('Shear production at IML4')
xlabel('year 2009')
ylabel('depth (m)')

hold on
[c, h] = contour(tv, pv, TM, [0 1 5 10], 'k');
clabel(c,h)
hold off

print('-depsc2', 'S2_IML4.eps')


figure(3)
clf
set(gcf,'PaperUnits','centimeters','PaperPosition',[10 10 30 20])
set(gcf, 'renderer', 'painters')
imagesc(tv(3900:4500), pv(1:end-1), S2(:,3900:4500))
colormap(jet);
caxis([0 0.01])
datetick('x', 7)
colorbar
title('Shear production at IML4')
xlabel('year 2009')
ylabel('depth (m)')

hold on
[c, h] = contour(tv(3900:4500), pv, TM(:,3900:4500), [0 1 5 10], 'k');
clabel(c,h)
hold off

figure(4)
clf
set(gcf,'PaperUnits','centimeters','PaperPosition',[10 10 30 20])
set(gcf, 'renderer', 'painters')
load Ri_colormap;
imagesc(tv(3000:5000), pv(1:end-1), Ri(:,3000:5000))
colormap(flipud(mycolormap));
caxis([0 1])
datetick('x')
colorbar
title('Richardson number at IML4')
xlabel('year 2009')
ylabel('depth (m)')

hold on
contour(tv(3000:5000), pv, TM(:,3000:5000), 'k')
hold off

