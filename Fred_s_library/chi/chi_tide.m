function chi_tide(chifiles, time2file, hab)

% I need to identify if a probe is missing....
% usage ex: chi_tide('chi_files', 'time_to_hightide.mat', 20)

        
% ------------------ compute Chi --------------------- % 
fid = fopen(chifiles);
C = textscan(fid, '%s', 'delimiter', '\n');
chi_files = char(C{1});

no_profiles = size(chi_files, 1); %number of profile* files 

load(time2file);
time2 = A;

count = 1;
clear Chi
for i = 1:no_profiles

    fname_in = chi_files(i, :);
    I = find(fname_in==' ');
    fname_in(I)=[];
    load(fname_in);

    I = find(isnan(chi1) == 1);
    chi1(I) = [];
    p_chi1(I) = [];
    
    %    I = find(p_chi1 > max(p_chi1-hab));
    I = find(p_chi1 > max(p_chi1)-hab);

    %    mean_chi = nanmean(chi1(I);
    Chi(count:count+length(I)-1) = chi1(I);
    Time2(count:count+length(I)-1) = time2(i);
    count = count+length(I);
end
% ------------------------------------------------------ % 

% ------------------ compute time relative to high tide --------------------- % 

% regular tide vector (-6, 6)
reg_tide = unique(round(time2));

reg_Chi = nan(round(length(Chi)/2), length(reg_tide)); % /2 its just approx
for i = 1: length(Time2)
    [Y, I] = min(abs(Time2(i)-reg_tide));
    J = find(isnan(reg_Chi(:, I))==1);
    reg_Chi(J(1), I)=Chi(i);
end
% ---------------------------------------------------------------------------- % 
%disp('Loess is too long (too much points! (2012-10-06))')
%keyboard

% --------------------- LOESS Statistical test ------------------------ %
II = find(isnan(Chi)==1);
Chi(II)=[]; Time2(II)=[];

% Wrap the timeserie
[Time2, I] = sort(Time2);
Chi = Chi(I);
TTime2 = [Time2(end:1) Time2 Time2(end:1)];
CChi = [Chi(end:1) Chi Chi(end:1)];
[Time2, I] = sort(TTime2);
CChi = CChi(I);
I = find(TTime2>-8 & TTime2<8);
Time2 = TTime2(I);
Chi = CChi(I);

% Loess runmean
smooth_chi = loess(Time2, Chi, sort(Time2), 0.3, 1);
time_chi = sort(Time2);
% $$$ 
% $$$ 
% $$$ II = find(isnan(mean_S2)==1);
% $$$ mean_S2(II)=[]; time2(II)=[];
% $$$ %smooth_shear2 = loess(time2, mean_S2, sort(time2), 0.3, 2);
% $$$ smooth_shear1 = loess(time2, mean_S2, sort(time2), 0.3, 1);
% $$$ time_shear = sort(time2);
% --------------------------------------------------------------------- %

% ------------------ Bootstrap --------------------- %
nboot = 500;
for i = 1:length(reg_tide)
    I = find(~isnan(reg_Chi(:,i)));
    samp = reg_Chi(I,i);
    NN = length(samp);
    clear chi_boot_b
    % create random sampling
    for b = 1:nboot
        r = rand(NN,1);
        r = ceil(r*NN/1);
        chi_boot_b(b) = nanmean(samp(r));
    end

    % mean of random sampling
    chi_boot_dot(i) = nanmean(chi_boot_b);

    % compute error (EFRON & GONG 83)
    chi_error(i) = sqrt(sum((diff([chi_boot_b chi_boot_dot(i)],1, 2)).^2, 2)./(nboot-1));
end


figure(1)
clf
set(gcf,'PaperUnits','centimeters','PaperPosition',[1 1 12 8])   

errorbar(reg_tide, chi_boot_dot, chi_error, '.k')
hold on
plot(time_chi, smooth_chi, 'k', 'linewidth', 2)
hold off
xlim([-6.5 6.5])
ylabel('\chi (^{\circ}C^2 s^{-1})')
xlabel('time versus high tide (hour)')
%set(gca, 'xticklabel', [])
set(gca, 'xgrid', 'on')
set(gca, 'fontsize', 10)
set(gca, 'yscale', 'log')
%ylim([.082 .09])

keyboard

print('-dpng', '-r300', 'chi_M2.png')
set(gcf, 'renderer', 'painters')
print('-depsc', 'chi_M2.eps')

% ---------------------------------------------------- %

