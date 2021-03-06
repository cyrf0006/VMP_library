function boot_NO3_all_stat(ctd_files, outfile_suffix, zmax, yearmin)

% function boot_NO3_all_stat(ctd_files, outfile_suffix, zmax, yearmin)
%
% This function performs bootstrap analysis and monthly average of
% profiles in a list. It also computes the confidence intervals on
% the monthly average. There is 3 main features:
% ex in /home/cyrf0006/PhD/CTD_MPO/Stat16:
%
% >  boot_NO3_all_stat('dat_profiles','stat16.dat', 500, 2000)
%
% 
%
% "ls -1 *.dat > dat_profiles"
%
% Author: Frederic Cyr, Oct. 2011
% 
% ------------------------------------------------------------- %

figure(1)
clf
set(gcf,'PaperUnits','centimeters','PaperPosition',[1 1 15 15])

% -- Some parameters -- %
zbin = 1; % would need to be adjusted if not already binned
zmin=1;
%zmax=450;
P = [zmin:zbin:zmax]';   
nboot = 500;

% -- Getting infos on profiles -- %
% load file names
fid = fopen(ctd_files);
C = textscan(fid, '%s', 'delimiter', '\n');
ctd = char(C{1});

n = datenum(ctd(:,1:11)); % in files from ODF2ASCII_ctd.sh

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

for j = 1:length(n)       

    file = load(ctd(j,:));
    
    % Quality control + save into matrix T and S
    if ~isempty(file)==1 % check if file empty
        
        if size(file, 2)~=5
            continue
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

        if length(pres)<=1 % need at least 2 pts for itp
            disp(sprintf('only one bottle, ignore cast no.%d', j))
        end
        
                
        % interp profiles and store in matrix        
        I = find(NO2==-99);
        p = pres;
        p(I) = [];
        NO2(I) = [];
        if length(p)>1
            NO2_mat(:,j) = interp1(p, NO2, P);
        end

        I = find(NO3==-99);
        p=pres;
        p(I) = [];
        NO3(I) = [];
        if length(p)>1        
            NO3_mat(:,j) = interp1(p, NO3, P);
        end
        
        I = find(PO4==-99);
        p=pres;
        p(I) = [];
        PO4(I) = [];
        if length(p)>1    
            PO4_mat(:,j) = interp1(p, PO4, P);
        end

    end
    
end  % for j

% -- Quality control -- %
% Remove empty columns
I = find(sum(~isnan(NO2_mat),1)>0);
NO2 = NO2_mat(:,I);
I = find(sum(~isnan(NO3_mat),1)>0);
NO3 = NO3_mat(:,I);
I = find(sum(~isnan(PO4_mat),1)>0);
PO4 = PO4_mat(:,I);

% $$$ % ml/l to mg/l from SeaBird note
% $$$ O2 = O2.*1.42903;
% $$$ % mg/l to umol/l
% $$$ O2 = O2.*1000/32;


% -- bootstrap -- %
disp('bootstrap...')
NNO2 = size(NO2,2); % NO2, NO3, PO4 may not be equal...
NNO3 = size(NO3,2); 
NPO4 = size(PO4,2); 
for b = 1:nboot        
    r = rand(NNO2,1);
    r = ceil(r*NNO2/1);
    NO2_boot_b(:,b) = nanmean(NO2(:,r),2);

    r = rand(NNO3,1);
    r = ceil(r*NNO3/1);
    NO3_boot_b(:,b) = nanmean(NO3(:,r),2);
    
    r = rand(NPO4,1);
    r = ceil(r*NPO4/1);
    PO4_boot_b(:,b) = nanmean(PO4(:,r),2);
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

% --- Max and Min of the distribution (not bootstrapped) --- %
NO2_sort = sort(NO2, 2);
for i = 1:size(NO2, 1)
    I = find(~isnan(NO2_sort(i, :))==1);
    if length(I) > 3
        good_bin = NO2_sort(i, I);
        CI_2p5 = round(2.5/100*length(I));
        CI_97p5 = round(97.5/100*length(I));
        if CI_2p5==0
            CI_2p5=1;
        end
        NO2_2p5_raw(i) = good_bin(CI_2p5);
        NO2_97p5_raw(i) = good_bin(CI_97p5);
    else
        NO2_2p5_raw(i) = NaN;
        NO2_97p5_raw(i) = NaN;
    end
end
NO2_2p5_raw = NO2_2p5_raw';
NO2_97p5_raw = NO2_97p5_raw';

NO3_sort = sort(NO3, 2);
for i = 1:size(NO3, 1)
    I = find(~isnan(NO3_sort(i, :))==1);
    if length(I) > 3
        good_bin = NO3_sort(i, I);
        CI_2p5 = round(2.5/100*length(I));
        CI_97p5 = round(97.5/100*length(I));
        if CI_2p5==0
            CI_2p5=1;
        end
        NO3_2p5_raw(i) = good_bin(CI_2p5);
        NO3_97p5_raw(i) = good_bin(CI_97p5);
    else
        NO3_2p5_raw(i) = NaN;
        NO3_97p5_raw(i) = NaN;
    end
end
NO3_2p5_raw = NO3_2p5_raw';
NO3_97p5_raw = NO3_97p5_raw';

PO4_sort = sort(PO4, 2);
for i = 1:size(PO4, 1)
    I = find(~isnan(PO4_sort(i, :))==1);
    if length(I) > 3
        good_bin = PO4_sort(i, I);
        CI_2p5 = round(2.5/100*length(I));
        CI_97p5 = round(97.5/100*length(I));
        if CI_2p5==0
            CI_2p5=1;
        end
        PO4_2p5_raw(i) = good_bin(CI_2p5);
        PO4_97p5_raw(i) = good_bin(CI_97p5);
    else
        PO4_2p5_raw(i) = NaN;
        PO4_97p5_raw(i) = NaN;
    end
end
PO4_2p5_raw = PO4_2p5_raw';
PO4_97p5_raw = PO4_97p5_raw';
% --------------------------------------------------------- %

figure(1)

% NO2 plot
subplot(1,3,1)
plot(NO2_ave, P, 'k', 'linewidth', 1)
hold on

x1 = NO2_97p5;
x2 = NO2_2p5;
x3 = NO2_97p5_raw;
x4 = NO2_2p5_raw;

nann = find(isnan(x3)==1);
if ~isempty(nann)==1 % remove NaNs
    
    x3(nann)=[];
    x4(nann)=[];
    PP = P;
    PP(nann)=[];
    patch([x3; flipud(x4); x3(1)], [PP; flipud(PP); PP(1)], [.9 .9 .9], 'edgecolor', 'none');
else
    patch([x3; flipud(x4); x3(1)], [P; flipud(P); P(1)], [.9 .9 .9], 'edgecolor', 'none');
end

nann = find(isnan(x1)==1);
if ~isempty(nann)==1 % remove NaNs
    
    x1(nann)=[];
    x2(nann)=[];
    PP = P;
    PP(nann)=[];
    patch([x1; flipud(x2); x1(1)], [PP; flipud(PP); PP(1)], [.8 .8 .8], 'edgecolor', 'none');
    plot(NO2_ave, P, 'k', 'linewidth', 1)
    hold off
else
    patch([x1; flipud(x2); x1(1)], [P; flipud(P); P(1)], [.8 .8 .8], 'edgecolor', 'none');
    plot(NO2_ave, P, 'k', 'linewidth', 1)
    hold off
end


set(gca, 'ydir', 'reverse')
set(gca, 'ygrid', 'on')
set(gca, 'xgrid', 'on')
set(gca, 'yminortick', 'on')
set(gca, 'xminortick', 'on')
set(gca, 'box', 'on')
%axis([0 15 0 zmax])
ylabel('P (dbar)')    
xlabel('NO2 (mmol m^{-3})')



% NO3 plot
subplot(1,3,2)
plot(NO3_ave, P, 'k', 'linewidth', 1)
hold on

x1 = NO3_97p5;
x2 = NO3_2p5;
x3 = NO3_97p5_raw;
x4 = NO3_2p5_raw;

nann = find(isnan(x3)==1);
if ~isempty(nann)==1 % remove NaNs
    
    x3(nann)=[];
    x4(nann)=[];
    PP = P;
    PP(nann)=[];
    patch([x3; flipud(x4); x3(1)], [PP; flipud(PP); PP(1)], [.9 .9 .9], 'edgecolor', 'none');
else
    patch([x3; flipud(x4); x3(1)], [P; flipud(P); P(1)], [.9 .9 .9], 'edgecolor', 'none');
end

nann = find(isnan(x1)==1);
if ~isempty(nann)==1 % remove NaNs
    
    x1(nann)=[];
    x2(nann)=[];
    PP = P;
    PP(nann)=[];
    patch([x1; flipud(x2); x1(1)], [PP; flipud(PP); PP(1)], [.7 .7 .7], 'edgecolor', 'none');
    plot(NO3_ave, P, 'k', 'linewidth', 1)
    hold off
else
    patch([x1; flipud(x2); x1(1)], [P; flipud(P); P(1)], [.7 .7 .7], 'edgecolor', 'none');
    plot(NO3_ave, P, 'k', 'linewidth', 1)
    hold off
end

set(gca, 'ydir', 'reverse')
set(gca, 'ygrid', 'on')
set(gca, 'xgrid', 'on')
set(gca, 'yminortick', 'on')
set(gca, 'xminortick', 'on')
set(gca, 'box', 'on')
%axis([50 450 0 zmax])
ylabel('P (dbar)')    
xlabel('NO3(mmol m^{-3})')


% PO4 plot
subplot(1,3,3)
plot(PO4_ave, P, 'k', 'linewidth', 1)
hold on

x1 = PO4_97p5;
x2 = PO4_2p5;
x3 = PO4_97p5_raw;
x4 = PO4_2p5_raw;

nann = find(isnan(x3)==1);
if ~isempty(nann)==1 % remove NaNs
    
    x3(nann)=[];
    x4(nann)=[];
    PP = P;
    PP(nann)=[];
    patch([x3; flipud(x4); x3(1)], [PP; flipud(PP); PP(1)], [.9 .9 .9], 'edgecolor', 'none');
else
    patch([x3; flipud(x4); x3(1)], [P; flipud(P); P(1)], [.9 .9 .9], 'edgecolor', 'none');
end

nann = find(isnan(x1)==1);
if ~isempty(nann)==1 % remove NaNs
    
    x1(nann)=[];
    x2(nann)=[];
    PP = P;
    PP(nann)=[];
    patch([x1; flipud(x2); x1(1)], [PP; flipud(PP); PP(1)], [.7 .7 .7], 'edgecolor', 'none');
    plot(PO4_ave, P, 'k', 'linewidth', 1)
    hold off
else
    patch([x1; flipud(x2); x1(1)], [P; flipud(P); P(1)], [.7 .7 .7], 'edgecolor', 'none');
    plot(PO4_ave, P, 'k', 'linewidth', 1)
    hold off
end

set(gca, 'ydir', 'reverse')
set(gca, 'ygrid', 'on')
set(gca, 'xgrid', 'on')
set(gca, 'yminortick', 'on')
set(gca, 'xminortick', 'on')
set(gca, 'box', 'on')
%axis([50 450 0 zmax])
ylabel('P (dbar)')    
xlabel('PO4(mmol m^{-3})')


outfileNO2 = ['boot_NO2_' outfile_suffix];
outfileNO3 = ['boot_NO3_' outfile_suffix];
outfilePO4 = ['boot_PO4_' outfile_suffix];
% save bootstrap mean profile
dlmwrite(outfileNO2, [P NO2_ave NO2_2p5 NO2_97p5 NO2_2p5_raw NO2_97p5_raw],'delimiter',' ','precision',6)
dlmwrite(outfileNO3, [P NO3_ave NO3_2p5 NO3_97p5 NO3_2p5_raw NO3_97p5_raw],'delimiter',' ','precision',6)
dlmwrite(outfilePO4, [P PO4_ave PO4_2p5 PO4_97p5 PO4_2p5_raw NO2_97p5_raw],'delimiter',' ','precision',6)

disp(sprintf('%d NO3 profiles used', size(NO3,2)))


keyboard
% save figure
figname = [outfile(1:end-4) '.png'];
print('-dpng', '-r300', figname)


% $$$ % save data
% $$$ outname = [outfile(1:end-4) '.mat'];
% $$$ save(outname, '')