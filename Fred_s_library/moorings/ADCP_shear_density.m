% Script that plot the shear between layer in function of the
% density profile (isopycnals instead of the depth). The beginning
% of the script comes from rich_IML4.m
%
% for the moment it must be run in rich_IML4 folder (because of
% tprofiles.dat and sprofiles.dat)
%
% author: F. Cyr - 29 sept. 2010
%
% ----------------------------------------------------- %

clear


%%%%%%%%%%%%%%%%%%
% ---- ADCP ----- %
%%%%%%%%%%%%%%%%%%00

%adcp_file = 'adcp_2009-09-29.mat';
adcp_file = 'IML4_ADCP09.mat';
%adcp_file = 'IML4_ADCP10.mat';
year = 2009;
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
yyyy(:) = year;
dd = str2num(datestr(mtime, 7));
hour = datestr(mtime, 13);
hh = str2num(hour(:,1:2));
mi = str2num(hour(:,4:5));

time_adcp = datenum(yyyy, mm, dd, hh, mi, 0);


%%%%%%%%%%%%%%%%%%%%%%%%%
% ---- T-S-density ---- %
%%%%%%%%%%%%%%%%%%%%%%%%%

% load T-S infos for 2009-2010
load('datprofiles.dat')
load('tprofiles.dat')
load('sprofiles.dat')

yyyy = str2num(datestr(datprofiles, 10));
I=find(yyyy==year);

% time vector
tvec = datprofiles(I);

% T-S matrix
Tmat = tprofiles(:,I);
Smat = sprofiles(:,I);
[X,Pmat] = meshgrid(1:size(Tmat,2), 1:300);

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
t2 = min([tvec(end); time_adcp(end)]);
%t2 = datenum(year, 06 ,28);
tv = t1:dt:t2;

%pressure vector
pv = z;

% regular density vector
%dmin = min(min(Dmat));
%dmax = max(max(Dmat));
dmin = 1021;
dmax = 1027;
dd = 0.1;
dv = [dmin:dd:dmax]';


% depth interpolation of density to ADCP resolution
for i = 1:size(Dmat,2)
    
    J = find(~isnan(Dmat(:,i))==1); %ignore NaNs
    DM_0(:,i) = interp1(Pmat(J,1), Dmat(J,i), pv);
    TM_0(:,i) = interp1(Pmat(J,1), Tmat(J,i), pv);
    
end

% time interpolation of density to ADCP resolution
for i = 1:size(DM_0,1)
    
    J = find(~isnan(DM_0(i,:))==1); %ignore NaNs
    DM(i,:) = interp1(tvec(J), DM_0(i,J), tv);
    TM(i,:) = interp1(tvec(J), TM_0(i,J), tv);
    
end

% restrict ADCP values to the same grid (pressure is OK, but time...)
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
    
    J = find(~isnan(u(i,:))==1); %ignore NaNs
    u(i,:) = interp1(tv(J), u(i,J), tv);

    J = find(~isnan(v(i,:))==1); %ignore NaNs
    v(i,:) = interp1(tv(J), v(i,J), tv);
    
    J = find(~isnan(dudz(i,:))==1); %ignore NaNs
    dudz(i,:) = interp1(tv(J), dudz(i,J), tv);
    
    J = find(~isnan(dvdz(i,:))==1); %ignore NaNs
    dvdz(i,:) = interp1(tv(J), dvdz(i,J), tv);
end

%remove zeros in S2, dudz, dvdz
I = find(S2==0);
S2(I)=NaN;

I = find(u==0);
u(I)=NaN;

I = find(v==0);
v(I)=NaN;

I = find(dudz==0);
dudz(I)=NaN;

I = find(dvdz==0);
dvdz(I)=NaN;


% here S2 is the shear matrix 
% (25 depth-cells of 4m bin and 8616 time step)
% DM matrix has almost the same size (26 cdepth-cells)
% The idea is to fill an empty matrix with values of S2
% corresponding to closest values of DM vs dv

S2d = nan(length(dv),size(S2,2));
Ud = nan(length(dv),size(u,2));
Vd = nan(length(dv),size(v,2));
DUd = nan(length(dv),size(dudz,2));
DVd = nan(length(dv),size(dvdz,2));
TMd = nan(length(dv),size(TM,2));


for i = 1:size(S2,2)
        
    I = nearestpoint(dv, DM(1:end-1,i));
    
    for j = 1:length(dv)
        if DM(I(j),i)-dv(j)<0.5
            S2d(j,i) = S2(I(j),i);
            Ud(j,i) = u(I(j),i);
            Vd(j,i) = v(I(j),i);
            DUd(j,i) = dudz(I(j),i);
            DVd(j,i) = dvdz(I(j),i);
        end
    end
    
    I = nearestpoint(dv, DM(:,i));

    for j = 1:length(dv)
        if DM(I(j),i)-dv(j)<0.5
            TMd(j,i) = TM(I(j),i);
        end        
    end
    
end

% while S2d is the matrix, S2d_mean is the average in time
S2d_mean=nanmean(S2d,2);
Ud_mean=nanmean(Ud,2);
Vd_mean=nanmean(Vd,2);
DUd_mean=nanmean(DUd,2);
DVd_mean=nanmean(DVd,2);
TMd_mean=nanmean(TMd,2);

% $$$ figure(1)
% $$$ clf
% $$$ set(gcf,'PaperUnits','centimeters','PaperPosition',[10 10 30 20])
% $$$ set(gcf, 'renderer', 'painters')
% $$$ imagesc(tv, dv, S2d)
% $$$ colormap(jet);
% $$$ caxis([0 0.01])
% $$$ datetick('x', 3)
% $$$ colorbar
% $$$ title('Shear production at IML4')
% $$$ xlabel('year 2009')
% $$$ ylabel('density (kg m^{-3})')
% $$$ 
% $$$ hold on
% $$$ [c, h] = contour(tv(3900:4500), pv, TM(:,3900:4500), [0 1 5 10], 'k');
% $$$ clabel(c,h)
% $$$ hold off

figure(2)
clf
set(gcf,'PaperUnits','centimeters','PaperPosition',[10 10 30 20])
set(gcf, 'renderer', 'painters')
%plot(DUd_mean, dv)
plot(Ud_mean, dv)
hold on
plot(Vd_mean, dv, 'r')
%plot(DVd_mean, dv, 'r')
%plot(S2d_mean, dv, '--k')
hold off
set(gca, 'ydir', 'reverse')
xlabel('mean velocity (m ^{s-1})')
ylabel('density (kg m^{-3})')
%legend('dudz', 'dvdz', 'S2')
legend('eastward', 'northward')

figure(3)
clf
set(gcf,'PaperUnits','centimeters','PaperPosition',[10 10 30 20])
set(gcf, 'renderer', 'painters')
xlabels{1} = '\rho (kg m^{-3})';
xlabels{2} = 'T (^{\circ}C)';
ylabels{1} = 'Depth(m)';
ylabels{2} = 'Depth(m)';
[ax,h1,h2] = plotxx(nanmean(DM,2), pv,nanmean(TM,2), pv, xlabels, ylabels);
set(ax, 'ydir', 'reverse')

figure(4)
clf
set(gcf,'PaperUnits','centimeters','PaperPosition',[10 10 30 20])
set(gcf, 'renderer', 'painters')
xlabels{1} = 'S (s^{-2})';
xlabels{2} = 'T (^{\circ}C)';
ylabels{1} = '\rho (kg m^{-3})';
ylabels{2} = '\rho (kg m^{-3})';
[ax,h1,h2] = plotxx(S2d_mean, dv, TMd_mean, dv, xlabels, ylabels);
set(ax, 'ydir', 'reverse')


figure(5)
clf
set(gcf,'PaperUnits','centimeters','PaperPosition',[10 10 30 20])
set(gcf, 'renderer', 'painters')
%plot(DUd_mean, dv)
plot(nanmean(u,2), pv)
hold on
plot(nanmean(v,2), pv, 'r')
hold off
set(gca, 'ydir', 'reverse')
xlabel('mean velocity (m ^{s-1})')
ylabel('z (m)')
%legend('dudz', 'dvdz', 'S2')
legend('eastward', 'northward')

figure(6)
clf
set(gcf,'PaperUnits','centimeters','PaperPosition',[10 10 30 20])
set(gcf, 'renderer', 'painters')
xlabels{1} = 'S (s^{-2})';
xlabels{2} = 'T (^{\circ}C)';
ylabels{1} = 'z (m)';
ylabels{2} = 'z (m)';
[ax,h1,h2] = plotxx(nanmean(S2,2), pv(1:end-1),nanmean(TM,2), pv , xlabels, ylabels);
set(ax, 'ydir', 'reverse')


figure(7)
clf
set(gcf,'PaperUnits','centimeters','PaperPosition',[10 10 10 10])

[valong vcross] = rotate_vec(nanmean(u,2), nanmean(v,2), 46);

xlabels{2} = 'U (m s^{-1})';
xlabels{1} = 'T (^{\circ}C)';
ylabels{1} = 'z (m)';
ylabels{2} = 'z (m)';
[ax,h1,h2] = plotxx(nanmean(TM,2), pv,valong, pv , xlabels, ylabels);
set(ax, 'ydir', 'reverse')

hold on
plot(vcross, pv, '--r')
hold off
legend('along-shore', 'cross-shore', 'location', 'southeast')

keyboard

print('-dpng', '-r300', 'along_cross_vel-and-T-vs-z.png')

set(gcf, 'renderer', 'painters')
print('-depsc2', 'along_cross_vel-and-T-vs-z.eps')



% $$$ figure(8)
% $$$ clf
% $$$ set(gcf,'PaperUnits','centimeters','PaperPosition',[10 10 10 10])
% $$$ 
% $$$ xlabels{1} = 'cross-shore velocity (m s^{-1})';
% $$$ xlabels{2} = 'T (^{\circ}C)';
% $$$ ylabels{1} = 'z (m)';
% $$$ ylabels{2} = 'z (m)';
% $$$ [ax,h1,h2] = plotxx(vcross, pv,nanmean(TM,2), pv , xlabels, ylabels);
% $$$ set(ax, 'ydir', 'reverse')
% $$$ 
% $$$ print('-dpng', '-r300', 'cross_vel-and-T-vs-z.png')
% $$$ 
% $$$ set(gcf, 'renderer', 'painters')
% $$$ print('-depsc2', 'cross_vel-and-T-vs-z.eps')


figure(9)
clf
set(gcf,'PaperUnits','centimeters','PaperPosition',[10 10 10 10])

[valong vcross] = rotate_vec(Ud_mean, Vd_mean, 46);

xlabels{2} = 'U (m s^{-1})';
xlabels{1} = 'T (^{\circ}C)';
ylabels{1} = '\rho (m^3 s^{-1})';
ylabels{2} = '\rho (m^3 s^{-1})';
[ax,h1,h2] = plotxx(TMd_mean, dv, valong, dv , xlabels, ylabels);
set(ax, 'ydir', 'reverse')

hold on
plot(vcross, dv, '--r')
hold off
legend('along-shore', 'cross-shore', 'location', 'southeast')


%% -- Test for rotation -- %%

for angle = 1:60
    
    [valong vcross] = rotate_vec(nanmean(u,2), nanmean(v,2), angle);
    
    A(angle, :) = [angle, max(valong)-min(valong)];
    
    angle
    
figure(7)
clf
set(gcf,'PaperUnits','centimeters','PaperPosition',[10 10 10 10])

[valong vcross] = rotate_vec(nanmean(u,2), nanmean(v,2), angle);

xlabels{2} = 'U (m s^{-1})';
xlabels{1} = 'T (^{\circ}C)';
ylabels{1} = 'z (m)';
ylabels{2} = 'z (m)';
[ax,h1,h2] = plotxx(nanmean(TM,2), pv,valong, pv , xlabels, ylabels);
set(ax, 'ydir', 'reverse')

hold on
plot(vcross, pv, '--r')
hold off
legend('along-shore', 'cross-shore', 'location', 'southeast')
pause
end
