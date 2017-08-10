function boot_o2_all_stat(ctd_files, outfile, zmax)

% function 
%
% This function performs bootstrap analysis and monthly average of
% profiles in a list. It also computes the confidence intervals on
% the monthly average. There is 3 main features:
% ex in /home/cyrf0006/PhD/CTD_MPO/Stat16:
%
% >  boot_o2_all_stat('dat_profiles','boot_T-O2_stat16_2005-11.dat', 300)
%
% 
%
% "ls -1 *.bin > ctd_files"
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

no_files = size(ctd,1); 

n = datenum(ctd(:,1:11)); % in files from ODF2ASCII_ctd.sh



for j = 1:length(n)       
    
    file = load(ctd(j,:));
    
    % Quality control + save into matrix T and S
    if ~isempty(file)==1 % check if file empty
        
        if size(file, 2)~=5
            continue
        end
        
        pres = file(:,1);
        oxyg = file(:,5);
        temp = file(:,2);
        sal = file(:,3);
        
        %remove -99
        I = find(oxyg==-99);
        oxyg(I) = NaN;
        I = find(temp==-99);
        temp(I) = NaN;
        I = find(sal==-99);
        sal(I) = NaN;
        
        % bin profiles
        for k = 1:length(P)
            I = find(pres>(zmin+(k-1)*zbin) & pres<(zmin+(k)*zbin));
            oxyg_bin(k, :) = nanmean(oxyg(I));
            temp_bin(k, :) = nanmean(temp(I));
            sal_bin(k, :) = nanmean(sal(I));

        end
        oxyg = oxyg_bin;
        temp = temp_bin;
        sal = sal_bin;
        pres = [round(min(pres)):zbin:round(max(pres))]';
        
        if max(pres) >= max(P) & min(pres)<=min(P)
            I = find(pres>=min(P) & pres<=max(P));
            
            % fill matrix with empty value
            O2(P,j)=NaN;
            T(P,j)=NaN;
            S(P,j)=NaN;
            
            
            O2(pres(I),j) = oxyg(pres(I));
            T(pres(I),j) = temp(pres(I));
            S(pres(I),j) = sal(pres(I));

            
        else % CTD profile doesnt start at p=1m or pmax smaller than 300m
            
            I = find(pres>=max(min(pres), min(P)) & pres<=min(max(pres), max(P)));
            
            % fill matrix with empty value
            O2(P,j)=NaN;
            T(P,j)=NaN;
            S(P,j)=NaN;

            % substituate empty value where it is possible
            O2(pres(I),j) = oxyg(pres(I));
            T(pres(I),j) = temp(pres(I));
            S(pres(I),j) = sal(pres(I));

        end
        
    else %file empty
        O2(P,j)=NaN;
        T(P,j)=NaN;
        S(P,j)=NaN;
    end
    
end  % for j

% -- Quality control -- %
% Remove empty columns
I = find(sum(~isnan(O2),1)>0);
O2 = O2(:,I);
I = find(sum(~isnan(T),1)>0);
T = T(:,I);
S = S(:,I);
% ml/l to mg/l from SeaBird note
O2 = O2.*1.42903;
% mg/l to umol/l
O2 = O2.*1000/32;

% Eliminate P > 450m
I = find(P>458);
O2(I,:) = nan(zmax-450, size(O2,2));
T(I,:) = nan(zmax-450, size(T,2));
S(I,:) = nan(zmax-450, size(S,2));


% vertical interpolation to remove NaNs
for j = 1:size(O2,2)
    I = find(isnan(O2(:,j))==1);
    if ~isempty(I)==1 
        p = P;
        o = O2(:,j);
        o(I) = [];
        p(I) = [];
        O2(:,j) = interp1(p, o, P);       
    end
    I = find(isnan(T(:,j))==1);
    if ~isempty(I)==1 
        p = P;
        t = T(:,j);
        t(I) = [];
        p(I) = [];
        T(:,j) = interp1(p, t, P);       
    end
     I = find(isnan(S(:,j))==1);
    if ~isempty(I)==1 
        p = P;
        s = S(:,j);
        s(I) = [];
        p(I) = [];
        S(:,j) = interp1(p, s, P);       
    end
end

% Remove column with only zeros
I = find(nansum(O2, 1)==0);
O2(:,I) = [];
I = find(nansum(T, 1)==0);
T(:,I) = [];
I = find(nansum(S, 1)==0);
S(:,I) = [];

% -- bootstrap -- %
disp('bootstrap...')
NO = size(O2,2); % NT and NS may not be equal...
NT = size(T,2); 
NS = size(S,2); 
for b = 1:nboot        
    r = rand(NO,1);
    r = ceil(r*NO/1);
    O2_boot_b(:,b) = nanmean(O2(:,r),2);
    
    r = rand(NT,1);
    r = ceil(r*NT/1);
    T_boot_b(:,b) = nanmean(T(:,r),2);
    
    r = rand(NS,1);
    r = ceil(r*NS/1);
    S_boot_b(:,b) = nanmean(S(:,r),2);
end




% --- mean values --- %
O2_boot_dot = nanmean(O2_boot_b,2);
T_boot_dot = nanmean(T_boot_b,2);
S_boot_dot = nanmean(S_boot_b,2);


% --- Compute 95% confidence interval --- %
O2_sort = sort(O2_boot_b, 2);
T_sort = sort(T_boot_b, 2);
S_sort = sort(S_boot_b, 2);

CI_2p5 = round(2.5/100*nboot);
CI_97p5 = round(97.5/100*nboot);

O2_2p5 = O2_sort(:,CI_2p5);
O2_97p5 = O2_sort(:,CI_97p5);
O2_ave = O2_boot_dot;
T_2p5 = T_sort(:,CI_2p5);
T_97p5 = T_sort(:,CI_97p5);
T_ave = T_boot_dot;
S_2p5 = S_sort(:,CI_2p5);
S_97p5 = S_sort(:,CI_97p5);
S_ave = S_boot_dot;

% --- Max and Min of the distribution (not bootstrapped) --- %
T_sort = sort(T, 2);
for i = 1:size(T, 1)
    I = find(~isnan(T_sort(i, :))==1);
    if length(I) > 3
        good_bin = T_sort(i, I);
        CI_2p5 = round(2.5/100*length(I));
        CI_97p5 = round(97.5/100*length(I));
        if CI_2p5==0
            CI_2p5=1;
        end
        T_2p5_raw(i) = good_bin(CI_2p5);
        T_97p5_raw(i) = good_bin(CI_97p5);
    else
        T_2p5_raw(i) = NaN;
        T_97p5_raw(i) = NaN;
    end
end
T_2p5_raw = T_2p5_raw';
T_97p5_raw = T_97p5_raw';

O2_sort = sort(O2, 2);
for i = 1:size(O2, 1)
    I = find(~isnan(O2_sort(i, :))==1);
    if length(I) > 3
        good_bin = O2_sort(i, I);
        CI_2p5 = round(2.5/100*length(I));
        CI_97p5 = round(97.5/100*length(I));
        if CI_2p5==0
            CI_2p5=1;
        end
        O2_2p5_raw(i) = good_bin(CI_2p5);
        O2_97p5_raw(i) = good_bin(CI_97p5);
    else
        O2_2p5_raw(i) = NaN;
        O2_97p5_raw(i) = NaN;
    end
end
O2_2p5_raw = O2_2p5_raw';
O2_97p5_raw = O2_97p5_raw';

S_sort = sort(S, 2);
for i = 1:size(S, 1)
    I = find(~isnan(S_sort(i, :))==1);
    if length(I) > 3
        good_bin = S_sort(i, I);
        CI_2p5 = round(2.5/100*length(I));
        CI_97p5 = round(97.5/100*length(I));
        if CI_2p5==0
            CI_2p5=1;
        end
        S_2p5_raw(i) = good_bin(CI_2p5);
        S_97p5_raw(i) = good_bin(CI_97p5);
    else
        S_2p5_raw(i) = NaN;
        S_97p5_raw(i) = NaN;
    end
end
S_2p5_raw = S_2p5_raw';
S_97p5_raw = S_97p5_raw';
% --------------------------------------------------------- %

figure(1)

%Temperature plot
subplot(1,2,1)
plot(T_ave, P, 'k', 'linewidth', 1)
hold on

x1 = T_97p5;
x2 = T_2p5;
x3 = T_97p5_raw;
x4 = T_2p5_raw;

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
    plot(T_ave, P, 'k', 'linewidth', 1)
    hold off
else
    patch([x1; flipud(x2); x1(1)], [P; flipud(P); P(1)], [.8 .8 .8], 'edgecolor', 'none');
    plot(T_ave, P, 'k', 'linewidth', 1)
    hold off
end


set(gca, 'ydir', 'reverse')
set(gca, 'ygrid', 'on')
set(gca, 'xgrid', 'on')
set(gca, 'yminortick', 'on')
set(gca, 'xminortick', 'on')
set(gca, 'box', 'on')
axis([0 15 0 zmax])
ylabel('P (dbar)')    
xlabel('T (^{\circ}C)')



% Oxygen plot
subplot(1,2,2)
plot(O2_ave, P, 'k', 'linewidth', 1)
hold on

x1 = O2_97p5;
x2 = O2_2p5;
x3 = O2_97p5_raw;
x4 = O2_2p5_raw;

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
    plot(O2_ave, P, 'k', 'linewidth', 1)
    hold off
else
    patch([x1; flipud(x2); x1(1)], [P; flipud(P); P(1)], [.7 .7 .7], 'edgecolor', 'none');
    plot(O2_ave, P, 'k', 'linewidth', 1)
    hold off
end

% hypoxia level
hold on
%plot([2 2], [150 350], '--k')
plot([62.5 62.5], [150 350], '--k')
%plot([150 150], [150 500], '--k')
hold off

set(gca, 'ydir', 'reverse')
set(gca, 'ygrid', 'on')
set(gca, 'xgrid', 'on')
set(gca, 'yminortick', 'on')
set(gca, 'xminortick', 'on')
set(gca, 'box', 'on')
axis([50 450 0 zmax])
ylabel('P (dbar)')    
xlabel('diss. O_2(\mu mol/l)')


% save bootstrap mean profile
dlmwrite(outfile, [P T_ave T_2p5 T_97p5 T_2p5_raw T_97p5_raw O2_ave O2_2p5 O2_97p5 O2_2p5_raw O2_97p5_raw],'delimiter',' ','precision',6)

% Save salinity (for Nitrate paper)
save boot_S.mat P S_ave S_2p5 S_97p5

% save figure
figname = [outfile(1:end-4) '.png'];
print('-dpng', '-r300', figname)


