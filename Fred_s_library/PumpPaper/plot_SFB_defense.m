clear

ncast = [1:12];

Hmax = 115;

Method = 'linear';

%root_data = '/Users/db/Documents/Cours/ISMER/Oceano_Experimentale/2012/MS_2012/Data/Leg_3/CTD/';
root_data = '/home/cyrf0006/PhD/Nitrates/CTD_station_fixe_MS_2012/';

fname_prefix = 'MS2012_SFB';

BM1 = datenum(2012,09,23,07,23,00);
HM1 = datenum(2012,09,23,13,05,00);
BM2 = datenum(2012,09,23,19,43,00);

year  = 2012;
month = 09;
day   = 29;

dt = 60;
dz = 0.5;

TIME   = [];
ZZ     = [];
TT     = [];
SS     = [];
O2O2   = [];
TRTR   = [];
PARPAR = [];
FLFL   = [];

for cast = ncast

  %fname = [root_data,fname_prefix,num2str(cast,'%2.2d'),'_moy.cnv']
  fname = [root_data,fname_prefix,num2str(cast),'_moy.cnv']
  fid1 = fopen(fname);
  
  tline(1:5) = 'aaaaa';
  nheader = 0;
  while 1
      
    tline = fgetl(fid1);
    
    if length(tline) >= 17 & tline(1:17) == '* NMEA UTC (Time)';
        year  = tline(28:31);
        month = tline(21:23);
        day   = tline(25:26);
        time  = tline(33:40);
        datestring  = ([day,'-',month,'-',year,' ',time])
        datestring0 = ([day,'-',month,'-',year,' ','00:00:00'])
        mtime  = datenum(datestring);
        mtime0 = datenum(datestring0);
        time = (mtime - mtime0)*86400;
    elseif length(tline) >= 17 & tline(1:15) == '* NMEA Latitude';
        lat  = str2num(tline(19:20)) + str2num(tline(22:26))/60;
    elseif length(tline) >= 17 & tline(1:16) == '* NMEA Longitude';
        lon  = str2num(tline(20:22)) + str2num(tline(24:29))/60;
    end
    
    if tline(1:5) == '*END*'
      fclose(fid1);
      nheader = nheader + 1;
      break;
    else
      nheader = nheader + 1;
    end
  end
  %data = importdata([root_data,fname_prefix,num2str(cast,'%2.2d'),'_moy.cnv'],' ',nheader);
  data = importdata([root_data,fname_prefix,num2str(cast),'_moy.cnv'],' ',nheader);

  z   = data.data(:,1);
  kmax = find(z == max(z));
  z = z(1:kmax);
  zmax = z(kmax);
  
  T   = data.data(1:kmax,2);
  S   = data.data(1:kmax,4);
  P   = data.data(1:kmax,13);
  % By-pass salinity by the Absolute Salinity - FC
  [SA, in_ocean] = gsw_SA_from_SP(S,P,lon,lat);
  S = SA;
  
  
  O2  = data.data(1:kmax,6);
  
  Tr    = data.data(1:kmax,9);
  PAR   = data.data(1:kmax,8);
  Fl    = data.data(1:kmax,7);
  time  = repmat(time,length(z),1);
  
  TIME   = [TIME;time];
  ZZ     = [ZZ;z];
  TT     = [TT;T];
  SS     = [SS;S];
  O2O2   = [O2O2;O2];
  TRTR   = [TRTR;Tr];
  PARPAR = [PARPAR;PAR];
  FLFL   = [FLFL;Fl];
  
end

% ---- Recover Nitrates data from bottles ---- %
disp('Recover Nutrients data from bottles')
data = load('NO3_SFB.dat');

% Average all bottles at same station (uncommoent below to check)
I = find(diff(data(:,1))==0 & diff(data(:,2))==0); 
while ~isempty(I)
    data(I(1),:) = nanmean(data(I:I+1,:), 1);
    data(I(1)+1, :) = [];
    %    disp(sprintf('Averaged lines %d and %d', I(1), I(1)+1))
    I = find(diff(data(:,1))==0 & diff(data(:,2))==0);     
end

% loop on data to build vectors
timeVec = data(:,1); % not time yet!
zVec = data(:,2);
PO4Vec = data(:,3);
NO3Vec = data(:,4);
SiVec = data(:,5);
[castTime, Icast] = unique(TIME);

% put time instead station no.
if length(castTime) ~= timeVec(end)
    disp('check, maybe error... [K]')
    keyboard
end

for i = 1:length(castTime)
    I = find(timeVec==i);
    if ~isempty(I) == 1
        timeVec(I) = castTime(i);
    end
end
% --------------------------------------------- %

% Griddatat Interpolation
timei = [min(TIME):dt:max(TIME)];
zi = [0:dz:max(ZZ)];
  
Ti  = griddata(TIME,ZZ,TT,timei',zi,Method); % from cast
Si  = griddata(TIME,ZZ,SS,timei',zi,Method);
O2i = griddata(TIME,ZZ,O2O2,timei',zi,Method);
Tri = griddata(TIME,ZZ,TRTR,timei',zi,Method);
Fli = griddata(TIME,ZZ,FLFL,timei',zi,Method);

Ni  = griddata(timeVec,zVec,NO3Vec,timei',zi,Method); % from bottles
Pi  = griddata(timeVec,zVec,NO3Vec,timei',zi,Method);
Sii  = griddata(timeVec,zVec,NO3Vec,timei',zi,Method);


% Interpolation by columns and lines
profTime = unique(timeVec); 
Nitp = nan(length(zi), length(profTime));
for i = 1:length(profTime)
    I = find(timeVec==profTime(i));
    Nitp(:,i) = interp1(zVec(I), NO3Vec(I), zi, 'cubic'); 
end

Nitp2 = nan(length(zi), length(timei));
for i = 1:length(zi)
    Nitp2(i,:) = interp1(profTime, Nitp(i,:), timei, 'cubic');
end

% by-pass griddata interpolation
disp('By-pass griddata interpolation for nitrates')
Ni=Nitp2;

ms = 1.0;
fs = 18;
nint = 20;
zbottom = Hmax;
timei = timei/3600;

HM1 = (HM1 - mtime0)*86400/3600;
BM1 = (BM1 - mtime0)*86400/3600;
%HM2 = (HM2 - mtime0)*86400/3600;


orient tall

%x0 = 0.1;
%y0 = 0.8;
%dy = 0.02
%length_fig = 0.9;
%height_fig = 0.15;

x0 = 0.1;
y0 = 0.7;
dy = 0.02;
length_fig = 0.875;
height_fig = 0.27;


xtext = 11.5;
ytext = 110;

lw = 0.5;



h = figure(1);
clf
set(h,'PaperUnits','centimeters','PaperPosition',[1 1 18 18])
% *********************** Adjust_space.m ************************ %
% Fields required by the function adjust_space.m. Please fill every
% of the following and call "adjust_space" in the script whenever
% you want. Do not touch four last fields
ncol = 1; % no. subplot column
nrow = 3; % no. subplot row
dx = 0.03 ; % horiz. space between subplots
dy = 0.03; % vert. space between subplots
lefs = 0.12; % very left of figure
rigs = 0.15; % very right of figure
tops = 0.01; % top of figure
bots = 0.12; % bottom of figure
figw = (1-(lefs+rigs+(ncol-1)*dx))/ncol;
figh = (1-(tops+bots+(nrow-1)*dy))/nrow;
count_col = 1;
count_row = 1;
% *************************************************************** %
cbVertOffset = .05;

s1 = subplot(311);
%contourf(timei,zi,Ti,[1:0.25:10],'linestyle','none');
contourf(timei,zi,Ti,50,'linestyle','none');
hold on
contour(timei,zi,Ti,[3 4 4.5 5],'k', 'linewidth', 2);
%contour(timei,zi,Ti,[1:0.25:4,5:2:10],'k','linewidth',lw);
plot(TIME/3600,ZZ,'k.','markersize',ms);
set(gca,'ydir','reverse','fontname','Times','fontsize',fs,'tickdir','out');
set(gca,'XMinorTick','on','YMinorTick','on','XTicklabel',[]);
axis([timei(1)-0.25 timei(end)+0.25 0 zbottom]);
%text(xtext,ytext,'T (^oC)','Color',[1 1 1],'fontname','Times','fontsize',20);
hold off
cb1 = colorbar('fontsize', fs);
ti = ylabel(cb1,'T (^{\circ}C)', 'FontSize', 18);
adjust_space
% $$$ ti_pos = get(ti, 'pos');
% $$$ ti_pos(1) = 10;
% $$$ set(ti, 'pos', ti_pos)
drawnow

cbPos = get(cb1, 'pos');
Pos = get(s1, 'pos');
cbPos(1) = Pos(1)+Pos(3)*1.02;
%cb1VertOffset = cb1Pos(4)-cb1Pos(4)*.9;
cbPos(2) = cbPos(2)+cbVertOffset/2;
cbPos(3) = cbPos(3)*.5;
cbPos(4) = cbPos(4)-cbVertOffset;
set(cb1, 'pos', cbPos)


s2 = subplot(312);
contourf(timei,zi,Si,50,'linestyle','none');
hold on
contour(timei,zi,Si,[29.5 30 31 32 32.5],'k','linewidth',2);
plot(TIME/3600,ZZ,'k.','markersize',ms);
set(gca,'ydir','reverse','fontname','Times','fontsize',fs,'tickdir','out');
set(gca,'XMinorTick','on','YMinorTick','on','XTicklabel',[]);
axis([timei(1)-0.25 timei(end)+0.25 0 zbottom]);
ylabel('Depth (m)','fontname','Times','fontsize',fs)
%text(xtext,ytext,'S (psu)','Color',[1 1 1],'fontname','Times','fontsize',20);
hold off
cb2 = colorbar('fontsize', fs);
ti = ylabel(cb2,'S_A (g kg^{-1})', 'FontSize', 18);
adjust_space
ti_pos = get(ti, 'pos');
ti_pos(1) = 10;
set(ti, 'pos', ti_pos)
drawnow

cbPos = get(cb2, 'pos');
Pos = get(s2, 'pos');
cbPos(1) = Pos(1)+Pos(3)*1.02;
%cbVertOffset = cbPos(4)-cbPos(4)*.9;
cbPos(2) = cbPos(2)+cbVertOffset/2;
cbPos(3) = cbPos(3)*.5;
cbPos(4) = cbPos(4)-cbVertOffset;
set(cb2, 'pos', cbPos)



s3 = subplot(313);
%contourf(timei,zi,Ni,50,'linestyle','none');
contourf(timei,zi,Ni,[nanmin(nanmin(Ni)):.5:nanmax(nanmax(Ni))],'linestyle','none');
hold on
contour(timei,zi,Ni,[5 10 15],'k','linewidth',2);
%plot(TIME/3600,ZZ,'k.','markersize',ms);
plot(timeVec/3600,zVec,'*','markersize',8, 'color', 'k', 'linewidth', 1);
set(gca,'ydir','reverse','fontname','Times','fontsize',fs,'tickdir','out');
set(gca,'XMinorTick','on','YMinorTick','on');
axis([timei(1)-0.25 timei(end)+0.25 0 zbottom]);
xlabel('Time (UTC) - 29 Sept. 2012','fontname','Times','fontsize',fs);
%text(xtext,ytext,'NO_3 (mmol m^{-3})','Color','w','fontname','Times','fontsize',20);
hold off
cb3 = colorbar('fontsize', fs);
ti = ylabel(cb3,'[NO_3] (mmol m^{-3})', 'FontSize', 18);
adjust_space
ti_pos = get(ti, 'pos');
ti_pos(1) = 10;
set(ti, 'pos', ti_pos)
drawnow

cbPos = get(cb3, 'pos');
Pos = get(s3, 'pos');
cbPos(1) = Pos(1)+Pos(3)*1.02;
%cbVertOffset = cbPos(4)-cbPos(4)*.9;
cbPos(2) = cbPos(2)+cbVertOffset/2;
cbPos(3) = cbPos(3)*.5;
cbPos(4) = cbPos(4)-cbVertOffset;
set(cb3, 'pos', cbPos)

print(h, '-dpng','-r300',['SFB_NO3_defense.png']);


% $$$ axes('Position',[x0,y0,length_fig,height_fig]);
% $$$ %contourf(timei,zi,Ti,[1:0.25:10],'linestyle','none');
% $$$ contourf(timei,zi,Ti,50,'linestyle','none');
% $$$ hold
% $$$ %contour(timei,zi,Ti,[1:0.25:4,5:2:10],'w','linewidth',lw);
% $$$ plot(TIME/3600,ZZ,'k.','markersize',ms);
% $$$ set(gca,'ydir','reverse','fontname','Times','fontsize',fs,'tickdir','out');
% $$$ set(gca,'XMinorTick','on','YMinorTick','on','XTicklabel',[]);
% $$$ axis([timei(1)-0.25 timei(end)+0.25 0 zbottom]);
% $$$ ylabel('Profondeur (m)','fontname','Times','fontsize',fs)
% $$$ text(xtext,ytext,'T (^oC)','Color','k','fontname','Times','fontsize',20);
% $$$ colorbar;
% $$$ text(HM1,-15,'HM','fontname','Times','fontsize',14,'HorizontalAlignment','center');
% $$$ text(BM1,-15,'BM','fontname','Times','fontsize',14);
% $$$ %text(HM2,-15,'HM','fontname','Times','fontsize',14);
% $$$ 
% $$$ axes('Position',[x0,y0-(height_fig+dy),length_fig,height_fig]);
% $$$ contourf(timei,zi,Si,50,'linestyle','none');
% $$$ hold;
% $$$ %contour(timei,zi,Si,[18:1:27,27.25:0.25:32],'w');
% $$$ plot(TIME/3600,ZZ,'k.','markersize',ms);
% $$$ set(gca,'ydir','reverse','fontname','Times','fontsize',fs,'tickdir','out');
% $$$ set(gca,'XMinorTick','on','YMinorTick','on','XTicklabel',[]);
% $$$ axis([timei(1)-0.25 timei(end)+0.25 0 zbottom]);
% $$$ ylabel('Profondeur (m)','fontname','Times','fontsize',fs)
% $$$ text(xtext,ytext,'S (psu)','Color','k','fontname','Times','fontsize',20);
% $$$ colorbar;
% $$$ 
% $$$ axes('Position',[x0,y0-2*(height_fig+dy),length_fig,height_fig]);
% $$$ %contourf(timei,zi,O2i,[30:10:120],'linestyle','none');
% $$$ contourf(timei,zi,O2i,[nanmin(nanmin(O2i)):10:nanmax(nanmax(O2i))],'linestyle','none');
% $$$ %contourf(timei,zi,O2i,[50:25:300],'linestyle','none');
% $$$ %contourf(timei,zi,O2i,[8:0.1:9]);
% $$$ hold;
% $$$ plot(TIME/3600,ZZ,'k.','markersize',ms);
% $$$ set(gca,'ydir','reverse','fontname','Times','fontsize',fs,'tickdir','out');
% $$$ set(gca,'XMinorTick','on','YMinorTick','on');
% $$$ axis([timei(1)-0.25 timei(end)+0.25 0 zbottom]);
% $$$ ylabel('Profondeur (m)','fontname','Times','fontsize',fs);
% $$$ xlabel('Heure (UTC) le 29 septembre 2012','fontname','Times','fontsize',fs);
% $$$ text(xtext,ytext,'O_2 (\mu mol/kg)','Color','k','fontname','Times','fontsize',20);
% $$$ colorbar;
% $$$ 
% $$$ 
% $$$ print('-dpng','-r300',['SFB.png']);

