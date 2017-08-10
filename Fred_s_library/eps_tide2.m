function eps_tide2(epsfiles, time2file, hab)

% I need to identify if a probe is missing....
% usage ex: eps_tide('eps_files', 'time_to_hightide.mat', 20)

%2012-11-16 copy-past of eps_tide, adapted.
        
% ------------------ compute Eps --------------------- % 
fid = fopen(epsfiles);
C = textscan(fid, '%s', 'delimiter', '\n');
eps_files = char(C{1});

no_profiles = size(eps_files, 1); %number of profile* files 

load(time2file);
time2 = A;

count = 1;
clear Eps
for i = 1:no_profiles

    fname_in = eps_files(i, :);
    I = find(fname_in==' ');
    fname_in(I)=[];
    load(fname_in);

    I = find(isnan(eps1) == 1);
    eps1(I) = [];
    p_eps1(I) = [];
    
    %    I = find(p_eps1 > max(p_eps1-hab));
    I = find(p_eps1 > max(p_eps1)-hab);

    %    mean_eps = nanmean(eps1(I);
    Eps(count:count+length(I)-1) = eps1(I);
    Time2(count:count+length(I)-1) = time2(i);
    count = count+length(I);
end
% ------------------------------------------------------ % 

% ------------------ compute time relative to high tide --------------------- % 

% regular tide vector (-6, 6)
reg_tide = unique(round(time2));

reg_Eps = nan(round(length(Eps)/2), length(reg_tide)); % /2 its just approx
for i = 1: length(Time2)
    [Y, I] = min(abs(Time2(i)-reg_tide));
    J = find(isnan(reg_Eps(:, I))==1);
    reg_Eps(J(1), I)=Eps(i);
end
% ---------------------------------------------------------------------------- % 
%disp('Loess is too long (too much points! (2012-10-06))')
%keyboard

% --------------------- LOESS Statistical test ------------------------ %
II = find(isnan(Eps)==1);
Eps(II)=[]; Time2(II)=[];

% $$$ % Wrap the timeserie
% $$$ [Time2, I] = sort(Time2);
% $$$ Eps = Eps(I);
% $$$ TTime2 = [Time2(end:1) Time2 Time2(end:1)];
% $$$ CEps = [Eps(end:1) Eps Eps(end:1)];
% $$$ [Time2, I] = sort(TTime2);
% $$$ CEps = CEps(I);
% $$$ I = find(TTime2>-8 & TTime2<8);
% $$$ Time2 = TTime2(I);
% $$$ Eps = CEps(I);

smooth_eps = loess(Time2, Eps, sort(Time2), 0.3, 1);
time_eps = sort(Time2);
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
    I = find(~isnan(reg_Eps(:,i)));
    samp = reg_Eps(I,i);
    NN = length(samp);
    clear eps_boot_b
    % create random sampling
    for b = 1:nboot
        r = rand(NN,1);
        r = ceil(r*NN/1);
        eps_boot_b(b) = nanmean(samp(r));
    end

    % mean of random sampling
    eps_boot_dot(i) = nanmean(eps_boot_b);

    % compute error (EFRON & GONG 83)
    eps_error(i) = sqrt(sum((diff([eps_boot_b eps_boot_dot(i)],1, 2)).^2, 2)./(nboot-1));
end


figure(1)
clf
set(gcf,'PaperUnits','centimeters','PaperPosition',[1 1 12 8])   

errorbar(reg_tide, eps_boot_dot, eps_error, '.k')
hold on
plot(time_eps, smooth_eps, 'k', 'linewidth', 2)
hold off
xlim([-6.5 6.5])
ylabel('\eps (^{\circ}C^2 s^{-1})')
xlabel('time versus high tide (hour)')
%set(gca, 'xticklabel', [])
set(gca, 'xgrid', 'on')
set(gca, 'fontsize', 10)
set(gca, 'yscale', 'log')
%ylim([.082 .09])

keyboard

print('-dpng', '-r300', 'eps_M2.png')
set(gcf, 'renderer', 'painters')
print('-depsc', 'eps_M2.eps')

% ---------------------------------------------------- %

