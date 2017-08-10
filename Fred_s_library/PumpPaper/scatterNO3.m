function scatterNO3(ctd_files, outfile_prefix, zbin, zmax, yearmin)
    
% scatterNO3('dat_profiles','scatterNO3_stat25', 10, 500, 2000)
   
figure(1)
clf
set(gcf,'PaperUnits','centimeters','PaperPosition',[1 1 12 15])
% *********************** Adjust_space.m ************************ %
% Fields required by the function adjust_space.m. Please fill every
% of the following and call "adjust_space" in the script whenever
% you want. Do not touch four last fields
ncol = 1; % no. subplot column
nrow = 1; % no. subplot row
dx = 0.03 ; % horiz. space between subplots
dy = 0.04; % vert. space between subplots
lefs = 0.1; % very left of figure
rigs = 0.05; % very right of figure
tops = 0.02; % top of figure
bots = 0.12; % bottom of figure
figw = (1-(lefs+rigs+(ncol-1)*dx))/ncol;
figh = (1-(tops+bots+(nrow-1)*dy))/nrow;
count_col = 1;
count_row = 1;
% *************************************************************** %


% -- Some parameters -- %
%zbin = 10; % would need to be adjusted if not already binned
zmin=zbin/2;
%zmax=450;
P = [zmin:zbin:zmax]';  
nboot = 500;

% -- Getting infos on profiles -- %
% load file names
fid = fopen(ctd_files);
C = textscan(fid, '%s', 'delimiter', '\n');
ctd = char(C{1});


n = datenum(ctd(:,1:11)); % in files from ODF2ASCII_bottle.sh

% year restriction
I = find(str2num(datestr(n, 10))<yearmin);
n(I) = [];
ctd(I,:) = [];

% month restriction (remove dec. to march)
I = find(str2num(datestr(n, 5)) < 4 | str2num(datestr(n, 5)) == 12);
count = length(I);
n(I) = [];
ctd(I,:) = [];


no_files = size(ctd,1); 
NO2_mat = nan(length(P), no_files);
NO3_mat = nan(length(P), no_files);
PO4_mat = nan(length(P), no_files);
zRawVec = [];
NO2RawVec = [];
NO3RawVec = [];
PO4RawVec = [];
for j = 1:length(n)       

    fname = ctd(j, :);
    I = find(fname==' ');   
    fname(I) = [];
    file = load(fname);
    
    % Quality control + save into matrix T and S
    if ~isempty(file)==1 % check if file empty
        if size(file, 2)~=5
            disp(['Are you sure you have the right file fortmat? in ' ctd(j,:)])
        end
            
        pres = file(:,1);
        NO2 = file(:,2);
        NO3 = file(:,3);
        PO4 = file(:,4);

        %remove -99
        I = find(pres==-99);
        pres(I) = []; %remove when no pressure
        NO2(I) = [];
        NO3(I) = [];
        PO4(I) = [];
        
        % store raw variables
        zRawVec = [zRawVec; pres];
        NO2RawVec = [NO2RawVec; NO2];
        NO3RawVec = [NO3RawVec; NO3];
        PO4RawVec = [PO4RawVec; PO4];
        
        % sort just in case (problem in Stat25, 24-APR-2001_14:44:27.dat)
        [pres, I] = sort(pres);
        NO2 = NO2(I);
        NO3 = NO3(I);
        PO4 = PO4(I);
        
        I = find(diff(pres)==0); %if 2 bottles at same depth
        pres(I) = []; % average them
        NO2(I) = nanmean(NO2([I I+1]));
        NO3(I) = nanmean(NO3([I I+1]));
        PO4(I) = nanmean(PO4([I I+1]));
        NO2(I+1) = []; % remove xtra entry
        NO3(I+1) = [];
        PO4(I+1) = [];
        
        for i = 1:length(pres)
            [Y, I] = min(abs(pres(i)-P));
            NO2_mat(I,j) = NO2(i);
            NO3_mat(I,j) = NO3(i);
            PO4_mat(I,j) = PO4(i);
        end
    end
end  % for j

% Bootles from MS2012
if strcmp(outfile_prefix, 'scatter_stat25') == 1
    disp('Add bottles from MS2012')
    data1 = load('~/PhD/Nitrates/CTD_station_fixe_MS_2012/NO3_SFA.dat');
    data2 = load('~/PhD/Nitrates/CTD_station_fixe_MS_2012/NO3_SFB.dat');
    % Average all bottles at same station (uncommoent below to check)
    I = find(diff(data1(:,1))==0 & diff(data1(:,2))==0); 
    while ~isempty(I)
        data1(I(1),:) = nanmean(data1(I:I+1,:), 1);
        data1(I(1)+1, :) = [];
        %    disp(sprintf('Averaged lines %d and %d', I(1), I(1)+1))
        I = find(diff(data1(:,1))==0 & diff(data1(:,2))==0);     
    end
    I = find(diff(data2(:,1))==0 & diff(data2(:,2))==0); 
    while ~isempty(I) 
        data2(I(1),:) = nanmean(data2(I:I+1,:), 1);
        data2(I(1)+1, :) = [];
        %    disp(sprintf('Averaged lines %d and %d', I(1), I(1)+1))
        I = find(diff(data2(:,1))==0 & diff(data2(:,2))==0);     
    end
    data = [data1; data2];

    % vectors
    MS.zVec = data(:,2);
    MS.PO4Vec = data(:,3);
    MS.NO3Vec = data(:,4);
    MS.SiVec = data(:,5);   

    NO3RawVec = [NO3RawVec; MS.NO3Vec];
    zRawVec = [zRawVec; MS.zVec];    
end


% -- Quality control -- %
[X, Z] = meshgrid(1:length(n), P);
zVecNO2 = Z(:);
zVecNO3 = Z(:);
zVecPO4 = Z(:);
NO2Vec = NO2_mat(:);
NO3Vec = NO3_mat(:);
PO4Vec = PO4_mat(:);

I = find(isnan(NO2Vec)==1);
NO2Vec(I) = [];
zVecNO2(I) = [];
I = find(NO2Vec==-99);
NO2Vec(I) = [];
zVecNO2(I) = [];

I = find(isnan(NO3Vec)==1);
NO3Vec(I) = [];
zVecNO3(I) = [];
I = find(NO3Vec==-99);
NO3Vec(I) = [];
zVecNO3(I) = [];

I = find(isnan(PO4Vec)==1);
PO4Vec(I) = [];
zVecPO4(I) = [];
I = find(PO4Vec==-99);
PO4Vec(I) = [];
zVecPO4(I) = [];

% For raw vectors
I = find(isnan(NO3RawVec)==1);
zRawVec(I) = [];
NO3RawVec(I) = [];
I = find(NO3RawVec==-99);
zRawVec(I) = [];
NO3RawVec(I) = [];



% -- bootstrap -- %
disp('bootstrap...')
NO2_boot_b = nan(length(P), nboot);
NO3_boot_b = nan(length(P), nboot);
PO4_boot_b = nan(length(P), nboot);
for i = 1:length(P)
    
    % previous version
    I = find(zVecNO2>=P(i)-zbin/2 & zVecNO2<P(i)+zbin/2);
    NO2 = NO2Vec(I);
% $$$     I = find(zVecNO3>=P(i)-zbin/2 & zVecNO3<P(i)+zbin/2);
% $$$     NO3 = NO3Vec(I);
    I = find(zVecPO4>=P(i)-zbin/2 & zVecPO4<P(i)+zbin/2);
    PO4 = PO4Vec(I);
        
    I = find(zRawVec>=P(i)-zbin/2 & zRawVec<P(i)+zbin/2);
    NO3 = NO3RawVec(I);
    
    NNO2 = length(NO2); % NO2, NO3, PO4 may not be equal...
    NNO3 = length(NO3);
    NPO4 = length(PO4);
    
    for b = 1:nboot        
        if NNO2>1
            r = rand(NNO2,1);
            r = ceil(r*NNO2/1);
            NO2_boot_b(i,b) = nanmean(NO2(r));
        end
        
        if NNO3>1
            r = rand(NNO3,1);
            r = ceil(r*NNO3/1);
            NO3_boot_b(i,b) = nanmean(NO3(r));
        end
        
        if NPO4>1
            r = rand(NPO4,1);
            r = ceil(r*NPO4/1);
            PO4_boot_b(i,b) = nanmean(PO4(r));
        end
    end
end



% --- mean values --- %
NO2_boot_dot = nanmean(NO2_boot_b,2);
NO3_boot_dot = nanmean(NO3_boot_b,2);
PO4_boot_dot = nanmean(PO4_boot_b,2);

% --- Compute 95% confidence interval --- %
NO2_sort = sort(NO2_boot_b, 2);
NO3_sort = sort(NO3_boot_b, 2);
PO4_sort = sort(PO4_boot_b, 2);

CI_2p5 = round(2.5/100*nboot);
CI_97p5 = round(97.5/100*nboot);

NO2_2p5 = NO2_sort(:,CI_2p5);
NO2_97p5 = NO2_sort(:,CI_97p5);
NO2_ave = NO2_boot_dot;

NO3_2p5 = NO3_sort(:,CI_2p5);
NO3_97p5 = NO3_sort(:,CI_97p5);
NO3_ave = NO3_boot_dot;

PO4_2p5 = PO4_sort(:,CI_2p5);
PO4_97p5 = PO4_sort(:,CI_97p5);
PO4_ave = PO4_boot_dot;
% -------------------------------------- %


% ----- Interpolate between average pts ----- %
I = find(~isnan(NO2_boot_dot)==1);
NO2_ave = interp1(P(I), NO2_boot_dot(I), P, 'linear');
I = find(~isnan(NO3_boot_dot)==1);
NO3_ave = interp1(P(I), NO3_boot_dot(I), P, 'linear');
I = find(~isnan(NO3_boot_dot)==1);
PO4_ave = interp1(P(I), PO4_boot_dot(I), P, 'linear');

I = find(~isnan(NO2_2p5)==1);
NO2_2p5 = interp1(P(I), NO2_2p5(I), P, 'linear');
I = find(~isnan(NO3_2p5)==1);
NO3_2p5 = interp1(P(I), NO3_2p5(I), P, 'linear');
I = find(~isnan(NO3_2p5)==1);
PO4_2p5 = interp1(P(I), PO4_2p5(I), P, 'linear');

I = find(~isnan(NO2_97p5)==1);
NO2_97p5 = interp1(P(I), NO2_97p5(I), P, 'linear');
I = find(~isnan(NO3_97p5)==1);
NO3_97p5 = interp1(P(I), NO3_97p5(I), P, 'linear');
I = find(~isnan(NO3_97p5)==1);
PO4_97p5 = interp1(P(I), PO4_97p5(I), P, 'linear');
% -------------------------------------- %


% ----- tanh fits ----- %
I = find(~isnan(NO3_boot_dot)==1);
S = [NO3_ave(I); NO3_2p5(I); NO3_97p5(I)]; Z=[P(I); P(I); P(I)];
Z3 = [P(I); P(I); P(I)];
Save = NO3_ave(I);
S2p5 = NO3_2p5(I); 
S97p5 = NO3_97p5(I); 
Z = P(I);
save testfit S Save S2p5 S97p5 Z Z3

clear S Z
global S Z
a0 = [5 10 .1 100];
load testfit
[Z, I] = sort(Z);  
Pfit = [round(min(Z)):1:round(max(Z))]';

S = Save(I);
options = optimset('MaxFunEvals', 1000, 'display', 'off');
[asol, funVal, exitFlag] = fminsearch ( 'fitFun10', a0, options);
Sfit_ave = asol(1) + asol(2).*tanh((Pfit-asol(3))./asol(4));
%plot(S3, Pfit, 'r--')
asol
S = S2p5(I);
options = optimset('MaxFunEvals', 1000, 'display', 'off');
[asol, funVal, exitFlag] = fminsearch ( 'fitFun10', a0, options);
Sfit2p5 = asol(1) + asol(2).*tanh((Pfit-asol(3))./asol(4));
%plot(S3, Pfit, 'r--')
asol
S = S97p5(I);
options = optimset('MaxFunEvals', 1000, 'display', 'off');
[asol, funVal, exitFlag] = fminsearch ( 'fitFun10', a0, options);
Sfit97p5 = asol(1) + asol(2).*tanh((Pfit-asol(3))./asol(4));
%plot(S3, Pfit, 'r--')
% ------------------------ %
asol


% -- Figure -- %
plot(NO3RawVec, zRawVec, '.k') 
hold on

% bin bootstrap
% $$$ x1 = NO3_2p5;
% $$$ x2 = NO3_97p5;
% $$$ I = find(~isnan(x1)==1 & ~isnan(x2)==1);
% $$$ patch([x1(I); flipud(x2(I)); x1(I(1))], [P(I); flipud(P(I)); P(I(1))], [1 1 1]*.6, 'edgecolor', 'none');

% tanh fit on bootstrapped bins
x1 = Sfit2p5;
x2 = Sfit97p5;
I = find(~isnan(x1)==1 & ~isnan(x2)==1);
patch([x1(I); flipud(x2(I)); x1(I(1))], [Pfit(I); flipud(Pfit(I)); Pfit(I(1))], [1 1 1]*.6, 'edgecolor', 'none');

plot(NO3RawVec, zRawVec, '.k') 
I = find(~isnan(NO3_boot_dot)==1);
herrorbar(NO3_ave(I), P(I), NO3_ave(I)-NO3_2p5(I), NO3_97p5(I)-NO3_ave(I), '.r');

% COR-1302 Nitrates profiles
if strcmp(outfile_prefix, 'scatter_stat22') == 1
    disp('Add profiles from COR-1302')
    fid = fopen('~/PhD/Nitrates/COR-1302/cor-1302.list');
    C = textscan(fid, '%s', 'delimiter', '\n');
    ctd = char(C{1});

    for i = 1:size(ctd,1)
        data = load(ctd(i,:));
        z = data(:,1);
        x = data(:,3);
        % 1-m bins
        dz = 1;
        zbin = dz/2:dz:max(z); 
        xbin = nan(size(zbin));
        for j = 1:length(zbin)
            I = find(z>=zbin(j)-dz/2 & z<zbin(j)+dz/2);
            xbin(j) = nanmean(x(I));
        end

        I = find(isnan(xbin)==1);
        xbin(I) = [];
        zbin(I) = [];
        if i==1 | i==2 | i==4 %| i==3 | i==5 | i==6
        %if  i==3 | i==5 | i==6
            % 1-m bins            
            plot(runmean(xbin, 10),zbin, 'm', 'linewidth', 2)
        end
    end        
    
    
    disp('Add bottles from COR-1203')
    data = load('~/PhD/Nitrates/COR-1302/BLT.dat');
    
    % Average all bottles at same station (uncommoent below to check)
    I = find(diff(data(:,1))==0 & diff(data(:,2))==0); 
    while ~isempty(I)
        data(I(1),:) = nanmean(data(I:I+1,:), 1);
        data(I(1)+1, :) = [];
        %    disp(sprintf('Averaged lines %d and %d', I(1), I(1)+1))
        I = find(diff(data(:,1))==0 & diff(data(:,2))==0);     
    end

    % loop on data to build vectors
    statVec = data(:,1)-1;
    zVec = data(:,2);
    PO4Vec = data(:,3);
    NO3Vec = data(:,4);
    SiVec = data(:,5);   

    for i = 1:length(zVec)
        if statVec(i) == 1 | statVec(i)==2 | statVec(i)==4
            %if statVec(i) == 3 | statVec(i)==5 | statVec(i)==6
            plot(NO3Vec(i), zVec(i), '*m', 'linewidth', 3)
        end
    end
        
end


% MS-2012 Nitrates bottles
if strcmp(outfile_prefix, 'scatter_stat25') == 1
    disp('Add bottles from MS2012')
    data = load('~/PhD/Nitrates/CTD_station_fixe_MS_2012/NO3_SFA.dat');
    
    % Average all bottles at same station (uncommoent below to check)
    I = find(diff(data(:,1))==0 & diff(data(:,2))==0); 
    while ~isempty(I)
        data(I(1),:) = nanmean(data(I:I+1,:), 1);
        data(I(1)+1, :) = [];
        %    disp(sprintf('Averaged lines %d and %d', I(1), I(1)+1))
        I = find(diff(data(:,1))==0 & diff(data(:,2))==0);     
    end

    % loop on data to build vectors
    zVec = data(:,2);
    PO4Vec = data(:,3);
    NO3Vec = data(:,4);
    SiVec = data(:,5);   

    for i = 1:length(zVec)
        plot(NO3Vec, zVec, '.m')
    end
end



set(gca, 'ydir', 'reverse')
hold off
xlabel('[NO_3] (mmol m^{-3})')
ylabel('Depth(m)')
ylim([0 350])
hold off




adjust_space


outfile = [outfile_prefix '.png'];
print('-dpng', outfile)
set(gcf, 'renderer', 'painters')
outfile = [outfile_prefix '.eps'];
print('-depsc2', outfile)

% save output for transect file
outfile = [outfile_prefix '.dat'];
dlmwrite(outfile, [P NO3_ave NO3_2p5 NO3_97p5],'delimiter',' ','precision',6)

% save smoothed profile
outfile = ['smooth_' outfile_prefix '.mat'];
save(outfile, 'Pfit', 'Sfit_ave', 'Sfit2p5', 'Sfit97p5', 'NO3RawVec', ...
     'zRawVec', 'P', 'NO3_ave', 'NO3_2p5', 'NO3_97p5', 'NO3_boot_dot')



% Fit from tanh on all gulf profile
% $$$ hold on
% $$$ asol = [15.5846 8.1804  112.9972   62.2792];
% $$$ asol = [15.8793    7.7332  114.9883   53.3929];
% $$$ S3 = asol(1) + asol(2).*tanh((Pfit-asol(3))./asol(4));
% $$$ plot(S3, Pfit, 'm--')
% $$$ hold off


