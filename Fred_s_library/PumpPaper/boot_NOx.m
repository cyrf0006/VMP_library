function boot_NOx(ctd_files, outfile_suffix, zmax, yearmin)

% function boot_NO3_all_stat(ctd_files, outfile_suffix, zmax, yearmin)
%
% This function performs bootstrap analysis and monthly average of
% profiles in a list. It also computes the confidence intervals on
% the monthly average. There is 3 main features:
% ex in /home/cyrf0006/PhD/CTD_MPO/Stat16:
%
% >  boot_NOx('dat_profiles','stat16.dat', 500, 1996)
%
% 
%
% "ls -1 *.dat > dat_profiles"
%
% Author: Frederic Cyr, Oct. 2011
% 
% ------------------------------------------------------------- %
global Z S

% -- Some parameters -- %
zbin = 10; % would need to be adjusted if not already binned
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
bvctd(I,:) = [];

no_files = size(ctd,1); 
NO2_mat = nan(length(P), no_files);
NO3_mat = nan(length(P), no_files);
PO4_mat = nan(length(P), no_files);
zRawVec = [];
NO2RawVec = [];
NO3RawVec = [];
PO4RawVec = [];
for j = 1:length(n)       

    file = load(ctd(j,:));
    
    % Quality control + save into matrix T and S
    if ~isempty(file)==1 % check if file empty
        
        if size(file, 2)~=5
            disp('Are you sure you have the right file fortmat?')
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


%plot(NO3Vec, zVec, '.k')


% -- bootstrap -- %
disp('bootstrap...')
NO2_boot_b = nan(length(P), nboot);
NO3_boot_b = nan(length(P), nboot);
PO4_boot_b = nan(length(P), nboot);
for i = 1:length(P)
    I = find(zVecNO2>=P(i)-zbin/2 & zVecNO2<P(i)+zbin/2);
    NO2 = NO2Vec(I);
    I = find(zVecNO3>=P(i)-zbin/2 & zVecNO3<P(i)+zbin/2);
    NO3 = NO3Vec(I);
    I = find(zVecPO4>=P(i)-zbin/2 & zVecPO4<P(i)+zbin/2);
    PO4 = PO4Vec(I);
    NNO2 = length(NO2); % NO2, NO3, PO4 may not be equal...
    NNO3 = length(NO3);
    NPO4 = length(PO4);
    
    for b = 1:nboot        
        r = rand(NNO2,1);
        r = ceil(r*NNO2/1);
        NO2_boot_b(i,b) = nanmean(NO2(r));
        
        r = rand(NNO3,1);
        r = ceil(r*NNO3/1);
        NO3_boot_b(i,b) = nanmean(NO3(r));
        
        r = rand(NPO4,1);
        r = ceil(r*NPO4/1);
        PO4_boot_b(i,b) = nanmean(PO4(r));
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


% -- Figure -- %
plot(NO3_ave, P, 'k')
hold on
x1 = NO3_2p5;
x2 = NO3_97p5;
I = find(~isnan(x1)==1 & ~isnan(x2)==1);
patch([x1(I); flipud(x2(I)); x1(I(1))], [P(I); flipud(P(I)); P(I(1))], [1 1 1]*.6, 'edgecolor', 'none');
plot(NO3RawVec, zRawVec, '.k') 
% $$$ plot(NO3_ave, P, 'k', 'linewidth', 2)
% $$$ plot(NO3_2p5, P, '--k', 'linewidth', 2)
% $$$ plot(NO3_97p5, P, '--k', 'linewidth', 2)
I = find(~isnan(NO2_boot_dot)==1);

% just a test...
S=[NO3_ave(I); NO3_2p5(I); NO3_97p5(I)]; Z=[P(I); P(I); P(I)]; 
save testfit S Z

herrorbar(NO3_ave(I), P(I), NO3_ave(I)-NO3_2p5(I), NO3_97p5(I)-NO3_ave(I), '.r');
set(gca, 'ydir', 'reverse')
hold off
xlabel('[NO_3] (mmol m^{-3})')
ylabel('Depth(m)')
ylim([0 350])


% Fit from tanh on this station  profile
hold on
clear S Z
global S Z
load testfit
[Z, I] = sort(Z);  
S = S(I);
a0 = [5 10 .1 100];
asol = fminsearch ( 'fitFun10', a0);
S3 = asol(1) + asol(2).*tanh((Z-asol(3))./asol(4));
plot(S3, Z, 'r--')
hold off

adjust_space

print('-dpng', 'scatter_NO3.png')
set(gcf, 'renderer', 'painters')
print('-depsc2', 'scatter_NO3.eps')

hold on


% $$$ % Fit from NOx-S relation
% $$$ load /home/cyrf0006/PhD/CTD_MPO/Stat25_10km/boot_S.mat
% $$$ fit = (S_ave-30.74)/.1403;
% $$$ plot(fit, P, 'k')

% Fit from tanh on all gulf profile
asol = [15.5846 8.1804  112.9972   62.2792];
asol = [15.8793    7.7332  114.9883   53.3929];
Z = zmin:zmax;
S3 = asol(1) + asol(2).*tanh((Z-asol(3))./asol(4));
plot(S3, Z, 'm--')

hold off
























% -- Some parameters -- %
zbin = 1; % would need to be adjusted if not already binned
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
bvctd(I,:) = [];

no_files = size(ctd,1); 
NO2_mat = nan(length(P), no_files);
NO3_mat = nan(length(P), no_files);
PO4_mat = nan(length(P), no_files);
zRawVec = [];
NO2RawVec = [];
NO3RawVec = [];
PO4RawVec = [];
for j = 1:length(n)       

    file = load(ctd(j,:));
    
    % Quality control + save into matrix T and S
    if ~isempty(file)==1 % check if file empty
        
        if size(file, 2)~=5
            disp('Are you sure you have the right file format?')
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
        
    end
end  % for j


% For raw vectors
I = find(isnan(NO2RawVec)==1);
zNO2 = zRawVec;
zNO2(I) = [];
NO2RawVec(I) = [];
I = find(NO2RawVec==-99);
zNO2(I) = [];
NO2RawVec(I) = [];

I = find(isnan(NO3RawVec)==1);
zNO3 = zRawVec;
zNO3(I) = [];
NO3RawVec(I) = [];
I = find(NO3RawVec==-99);
zNO3(I) = [];
NO3RawVec(I) = [];

I = find(isnan(PO4RawVec)==1);
zPO4 = zRawVec;
zPO4(I) = [];
PO4RawVec(I) = [];
I = find(PO4RawVec==-99);
zPO4(I) = [];
PO4RawVec(I) = [];



% -- bootstrap on tanh fit -- %
disp('bootstrap...')
NO2_boot_b = nan(nboot, 4); % 4 for four coefficients
NO3_boot_b = nan(nboot, 4);
PO4_boot_b = nan(nboot, 4);
NNO2 = length(NO2RawVec); % NO2, NO3, PO4 may not be equal...
NNO3 = length(NO3RawVec);
NPO4 = length(PO4RawVec);
count = 0;
count2 = 0;
for iboot = 1:nboot        
    % NO3 and PO4 only since NO2 has not the same pattern
    if mod(iboot, 100) == 0
        disp(sprintf('iteration no. %d / %d', iboot, nboot))
    end
    r = rand(NNO3,1);
    r = ceil(r*NNO3/1);
    NO3Vec = NO3RawVec(r);       
    [Z, I] = sort(zNO3(r));  
    S = NO3Vec(I);
    a0 = [5 10 .1 100];
    options = optimset('MaxFunEvals', 1500, 'display', 'off');
    [asol, funVal, exitFlag] = fminsearch ( 'fitFun10', a0, options);
    if exitFlag == 1
        NO3_boot_b(iboot,:) = asol;
    else
        count = count+1;
    end
    
        
    r = rand(NPO4,1);
    r = ceil(r*NPO4/1);
    PO4Vec = PO4RawVec(r);       
    [Z, I] = sort(zPO4(r));  
    S = PO4Vec(I);
    a0 = [1 1 .1 10];
    options = optimset('MaxFunEvals', 1500, 'display', 'off');
    [asol, funVal, exitFlag] = fminsearch ( 'fitFun10', a0, options);
    if exitFlag == 1
        PO4_boot_b(iboot,:) = asol;
    else
        count2 = count2+1;
    end
end

% Result idea
colorscale = 0:.1:1;
figure(1)
clf
count = 1;
plot(NO3RawVec, zNO3, '.k')
set(gca, 'ydir', 'reverse')
hold on
NO3Mat = nan(length(P), nboot);
for i = 1:nboot
    asol = NO3_boot_b(i,:); 
    S3 = asol(1) + asol(2).*tanh((P-asol(3))./asol(4));
    NO3Mat(:,i) = S3';
    
    plot(S3, P, 'color', [1 1 1].*colorscale(count))
    if count == length(colorscale)
        count = 1;
    else
        count = count+1;
    end

end
plot(NO3RawVec, zNO3, '.k')
hold off


disp(['The function is not finished since preliminary results are not satisfying.'])
keyboard


