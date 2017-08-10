function sl_yoyotimeseries(ctdFile)

% The name says everything
%usage ex:
% > sl_yoyotimeseries('SBE19plus_01906786_2013_06_10_0070.cnv')
% Could be easily modified to include multiple input files.    

% F. Cyr - june 2013

load(ctdFile);
    
dp = .5;
dt = 5/86400;   
PVec = 0.25:dp:50;
timeVec = min(mtime):dt:max(mtime);


P = data(:,2);
T = data(:,3);
S = data(:,6);


dP = diff(P);
I = find(dP<0);
P(I) = [];
T(I) = [];
S(I) = [];
mtime(I) = [];
I = find(P<1.5);
P(I) = [];
T(I) = [];
S(I) = [];
mtime(I) = [];
zmax = max(P);


Tgrid = griddata((mtime-mean(mtime))/std(mtime),(P-mean(P))/std(P),T,(timeVec-mean(mtime))/std(mtime),(PVec'-mean(P))/std(P));
Sgrid = griddata((mtime-mean(mtime))/std(mtime),(P-mean(P))/std(P),S,(timeVec-mean(mtime))/std(mtime),(PVec'-mean(P))/std(P));

figure(1)
clf
pcolor(timeVec, PVec,Tgrid)
shading interp
hold on
contour(timeVec, PVec,Tgrid, 'k')
set(gca, 'ydir', 'reverse')
datetick
hold on
plot(mtime, P, '.k')
ylim([0 zmax])
colorbar 

figure(2)
clf
pcolor(timeVec, PVec,Sgrid)
shading interp
hold on
contour(timeVec, PVec,Sgrid, 'k')
set(gca, 'ydir', 'reverse')
datetick
hold on
plot(mtime, P, '.k')
ylim([0 zmax])
colorbar
