function boot_VMP_O2(file_names, o2_bootfile, zbin, varargin)

% function boot_VMP_O2.m
%
% ex: boot_VMP_O2('file_names_riki', 'boot_T-O2_stat23_2005-11.dat', 5) 
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% *********************** Adjust_space.m ************************ %
% Fields required by the function adjust_space.m. Please fill every
% of the following and call "adjust_space" in the script whenever
% you want. Do not touch four last fields
ncol = 3; % no. subplot column
nrow = 1; % no. subplot row
dx = 0.03 ; % horiz. space between subplots
dy = 0.05; % vert. space between subplots
lefs = 0.1; % very left of figure
rigs = 0.05; % very right of figure
tops = 0.02; % top of figure
bots = 0.18; % bottom of figure
figw = (1-(lefs+rigs+(ncol-1)*dx))/ncol;%(do not touch!)
figh = (1-(tops+bots+(nrow-1)*dy))/nrow;%(do not touch!)
count_col = 1; %(do not touch!)
count_row = 1; %(do not touch!)
% *************************************************************** %

% Deal with varargin
if isempty(varargin)==1
    sav = 0;
    xtra_contour = 0;
elseif size(varargin,2)==1
    sav = 1;
    xtra_contour = 0;
    outfile = varargin{1};
elseif size(varargin,2)==2
    sav = 1;
    outfile = varargin{1};
    xtra_contour = 1;
    file_names2 = varargin{2};
else
    disp('Wrong input... try "help mean_K"')
    return
end

% some constants
GAMMA = 0.2; %mixing efficiency

% depth range for extracting the average
zmin = 0; % top of the fisrt bin
zmax = 340;
nboot = 500;    
P_bin = [zmin+zbin/2:zbin:zmax-zbin/2]';

% load eps_files names (file in which are recorded .P files)
fid = fopen(file_names);
C = textscan(fid, '%s', 'delimiter', '\n');
epsfiles = char(C{1});

siz = size(epsfiles);
no_profile = siz(1); %number of eps_files 


% raw matrix to fill
mat_eps_bin = nan(length(P_bin), no_profile);
mat_K_bin = mat_eps_bin;
mat_N2_bin = mat_eps_bin;
    
N2_count = 0;
N2_min = -6;

% ----------------------------------------- %
% -------- loop on epsilon  profiles ------ %
% ----------------------------------------- %
for profile = 1:no_profile

    fname = epsfiles(profile, :);
    I = find(fname==' ');   
    fname(I) = [];
    
    load(fname)

    %%%%%%%%%%%%%%%%%%%%%%
    % - EPSILON, N2, K - %
    %%%%%%%%%%%%%%%%%%%%%%
    
    % ---- average of 2 epsilon profiles ---- %
    if length(p_eps1) == length(p_eps2) % no problem, easy
        p_k = p_eps1;
        p_N = p_eps1;
        EPS2=nan(length(p_eps1),2);
        EPS2(:,1) = eps1;
        EPS2(:,2) = eps2;
    elseif min(p_eps1)<min(p_eps2) % p1 drive, nomatter size of p2!
        p_k = p_eps1;        
        p_N = p_eps1;
        EPS2=nan(length(p_eps1),2);
        EPS2(:,1) = eps1;
        ind = find(p_eps2 == p_eps1(1));
        EPS2(ind:length(p_eps2),2) = eps2;
    else % p2 drive, nomatter size of p1!
        p_k = p_eps2;
        p_N = p_eps1;
        EPS2=nan(length(p_eps2),2);
        EPS2(:,2) = eps2;
        ind = find(p_eps1 == p_eps2(1));
        EPS2(ind:length(p_eps1),1) = eps1;
    end
    
    % "selection" average (uncomment plotting portion to inspect)
    I = find(EPS2(:,1)>10*EPS2(:,2));
    EPS2(I,1)=NaN;
    I = find(EPS2(:,2)>10*EPS2(:,1));
    EPS2(I,2)=NaN;
   
    EPS = nanmean(EPS2,2); %MEAN EPSILON
                        
    
    % if only nans
    if sum(~isnan(EPS))==0;
        continue
    end
    
    % Homestyle despike
    [Y, No] = Dspike(EPS, 5, 8);

    % uses home despike
    EPS = Y;
    
    % ---- mean N2 ---- %
    N2=nan(length(p_k),1);
    if length(p_N)==length(p_k);
        N2(:) = N.^2;
    else
        N2(ind:length(p_N))=N.^2;
    end
    
   
    % ---- compute K_rho ---- %
    K_rho = GAMMA.*EPS./(N2);

   
    % ---- Remove unrealistic diffusivity ---- %
    I = find(log10(N2)<N2_min);
    if ~isempty(I)==1
        K_rho_ignored(N2_count+1:N2_count+length(I))=K_rho(I)';
        K_rho(I)=NaN;
        N2_count = N2_count+length(I);
    end
    
  
    % ---------- bin profile -------------- %
    for i = 1:length(P_bin)
        I = find(p_k >= P_bin(i)-zbin/2 & p_k <= P_bin(i)+zbin/2);
        mat_eps_bin(i, profile) = nanmean(EPS(I));
        mat_K_bin(i, profile) = nanmean(K_rho(I));
        mat_N2_bin(i, profile) = nanmean(N2(I));
    end
    % ------------------------------------ %
    
    % ------ remove last bins ------------- %
   
    if profile<814
    I = find(~isnan(mat_eps_bin(:,profile))==1);
    
    %maxP(profile) = P_bin(I(end));
    
    % keyboard
    no_remove = 4;
    if length(I)<no_remove
        mat_eps_bin(:, profile)=NaN;
        mat_K_bin(:, profile)=NaN;
        mat_N2_bin(:, profile)=NaN;
    elseif no_remove ~= 0
        mat_eps_bin(I(end-no_remove+1:end), profile)=NaN;
        mat_K_bin(I(end-no_remove+1:end), profile)=NaN;
        mat_N2_bin(I(end-no_remove+1:end), profile)=NaN;
    end
    end
    % ------------------------------------ %
end 


%%%%%%%%%%%%%%%%%%%%
% --- bootstrap --- %
%%%%%%%%%%%%%%%%%%%%
% just rename
mat_eps  = mat_eps_bin;
mat_K = mat_K_bin;
mat_N2 = mat_N2_bin;

% -- bootstrap -- %
disp('bootstrap...')
N = size(mat_eps,2);

% create random sampling
for b = 1:nboot
    r = rand(N,1);
    r = ceil(r*N/1);
    
    eps_boot_b(:,b) = nanmean(mat_eps(:,r),2);
    K_boot_b(:,b) = nanmean(mat_K(:,r),2);
    N2_boot_b(:,b) = nanmean(mat_N2(:,r),2);
end

% mean of random sampling
eps_boot_dot = nanmean(eps_boot_b,2);
K_boot_dot = nanmean(K_boot_b,2);
N2_boot_dot = nanmean(N2_boot_b,2);

% Compute 95% confidence interval
eps_sort = sort(eps_boot_b, 2);
K_sort  = sort(K_boot_b, 2);
N2_sort  = sort(N2_boot_b, 2);

CI_2p5 = round(2.5/100*nboot);
CI_97p5 = round(97.5/100*nboot);

eps_2p5 = eps_sort(:,CI_2p5);
K_2p5 = K_sort(:,CI_2p5);
N2_2p5 = N2_sort(:,CI_2p5);
eps_97p5 = eps_sort(:,CI_97p5);
K_97p5 = K_sort(:,CI_97p5);
N2_97p5 = N2_sort(:,CI_97p5);

% -------- Mean Profiles (3 methods) -------- %
% $$$ %  1) sort of trimmean
% $$$ eps_ave = nan(size(mat_eps,1),1);
% $$$ K_ave = nan(size(mat_eps,1),1);
% $$$ N2_ave = nan(size(mat_eps,1),1);
% $$$ 
% $$$ for i = 1:size(mat_eps,1)
% $$$     
% $$$     ep = mat_eps(i,:);
% $$$     ka = mat_K(i,:);
% $$$     en = mat_N2(i,:);
% $$$     
% $$$     I = find(isnan(ep));
% $$$ 
% $$$     if length(I)<length(ep)-4
% $$$         ep(I) = [];
% $$$         ka(I) = [];
% $$$         en(I) = [];
% $$$         ep = sort(ep);
% $$$         ka = sort(ka);
% $$$         en = sort(en);
% $$$         
% $$$         x1 = round(5/100*length(ep));
% $$$         x2 = round(95/100*length(ep));
% $$$ 
% $$$         if x1==0
% $$$             x1=1;
% $$$         end
% $$$         %[i x1 x2]
% $$$         eps_ave(i) = mean(ep(x1:x2)); 
% $$$         K_ave(i) = mean(ka(x1:x2)); 
% $$$         N2_ave(i) = mean(en(x1:x2));  
% $$$     end
% $$$ end


% 2) average of all profiles
% $$$ eps_ave = nanmean(mat_eps,2);
% $$$ K_ave = nanmean(mat_K,2);
% $$$ N2_ave = nanmean(mat_N2,2);

% 3) Mean bootstrapped profile
eps_ave = eps_boot_dot;
K_ave = K_boot_dot;
N2_ave = N2_boot_dot;
% ---------------------------------------------- %



% ----------------- Load O2 data ------------------- %
% format: [P T_ave T_2p5 T_97p5 T_2p5_raw T_97p5_raw O2_ave O2_2p5 O2_97p5 O2_2p5_raw O2_97p5_raw]
O2 = load(o2_bootfile); % already in umol/l
p_o2 = O2(:,1);
o2_mean = O2(:,7);
o2_2p5 = O2(:,8);
o2_97p5 = O2(:,9);

% bin to VMP resolution
for i = 1:length(P_bin)
    I = find(p_o2 >= P_bin(i)-zbin/2 & p_o2 <= P_bin(i)+zbin/2);
    O2_mean_bin(i, :) = nanmean(o2_mean(I));
    O2_2p5_bin(i, :) = nanmean(o2_2p5(I));
    O2_97p5_bin(i, :) = nanmean(o2_97p5(I));

end
% ----------------- end O2 ------------------------- %

% -----------------  O2 flux ------------------- %
DO = O2_mean_bin;
DO_2p5 = O2_2p5_bin;
DO_97p5 = O2_97p5_bin;
dodz = diff(DO)./zbin;
P_itp =  P_bin(1:end-1)+zbin/2;
Kitp = interp1(P_bin, K_ave, P_itp);
Fo2 = -Kitp.*dodz; %(in umol*m/l/s)
Fo2 = Fo2*3153600; % conversion to umol/cm2/yr
% ----------------- end O2 ------------------------- %



%%%%%%%%%%%%%%%%%%%%
% --- plotting --- %
%%%%%%%%%%%%%%%%%%%%


figure(2)    
clf
set(gcf,'PaperUnits','centimeters','PaperPosition',[10 10 14 8])

subplot(1, 3, 1);
plot(DO, P_bin, 'k', 'linewidth', 0.25);

% shade area
hold on
x2 = DO_2p5;
x1 = DO_97p5;

nann = find(isnan(x1)==1);
if ~isempty(nann)==1 % remove NaNs    
    x1(nann)=[];
    x2(nann)=[];
    PP = P_bin;
    PP(nann)=[];
    patch([x1; flipud(x2); x1(1)], [PP; flipud(PP); PP(1)], [.8 .8 .8], 'edgecolor', 'none');
else
    patch([x1; flipud(x2); x1(1)], [P_bin; flipud(P_bin); P_bin(1)], [.8 .8 .8], 'edgecolor', 'none');
end

plot(DO, P_bin, 'k', 'linewidth', 0.25);
plot([62.5 62.5], [250 340], '--k')
hold off
set(gca, 'ydir', 'reverse')
axis([50 350 0 340])
%set(gca, 'xtick', [1e-9 1e-8 1e-7 1e-6])
set(gca, 'xminortick', 'off')
ylabel('P (dbar)', 'FontSize', 9)
xlabel('[O_2] (\mu mol l^{-1})', 'FontSize', 9)
set(gca, 'FontSize', 9)
set(gca, 'xminortick', 'on')
set(gca, 'yminortick', 'on')
set(gca, 'tickdir', 'out')
set(gca, 'ygrid', 'on')

% Write letter identification
text(340, 310, 'a', ...
         'verticalalignment', 'bottom', 'horizontalalignment', 'right', 'fontsize', 9, 'fontweight','bold');


adjust_space

subplot(1, 3, 2)

semilogx(K_ave, P_bin, 'k', 'linewidth', 0.25)
hold on
%semilogx(K_ave*0+5.5e-5, P_bin, '--k', 'linewidth', 0.25)
x2 = K_2p5;
x1 = K_97p5;
nann = find(isnan(x1)==1);
if ~isempty(nann)==1 % remove NaNs    
    x1(nann)=[];
    x2(nann)=[];
    PP = P_bin;
    PP(nann)=[];
    patch([x1; flipud(x2); x1(1)], [PP; flipud(PP); PP(1)], [.8 .8 .8], 'edgecolor', 'none');
else
    patch([x1; flipud(x2); x1(1)], [P_bin; flipud(P_bin); P_bin(1)], [.8 .8 .8], 'edgecolor', 'none');
end
semilogx(K_ave, P_bin, 'k', 'linewidth', 0.25)
hold off

set(gca, 'ydir', 'reverse')
axis([5e-6 1e-3 0 340])
set(gca, 'xtick', [1e-5 1e-4 1e-3])
set(gca, 'xminortick', 'off')
xlabel('K (m^2 s^{-1})', 'FontSize', 9)
set(gca, 'FontSize', 9)
set(gca, 'yticklabel', [])
set(gca, 'xminortick', 'on')
set(gca, 'yminortick', 'on')
set(gca, 'tickdir', 'out')
set(gca, 'ygrid', 'on')

% Write letter identification
Ylim = [0 340];
text(5e-4, 310, 'b', ...
         'verticalalignment', 'bottom', 'horizontalalignment', 'right', 'fontsize', 9, 'fontweight','bold');

adjust_space


subplot(1, 3, 3)

semilogx(Fo2, P_itp, 'k', 'linewidth', 0.25)
hold on
plot([0 0], [0 340], '--k','linewidth', 0.25)
%semilogx(K_ave*0+5.5e-5, P_bin, '--k', 'linewidth', 0.25)
%x2 = K_2p5;
%x1 = K_97p5;

%patch([x1; flipud(x2); x1(1)], [P_bin; flipud(P_bin); P_bin(1)], [.8 .8 .8], 'edgecolor', 'none');
hold off
set(gca, 'ydir', 'reverse')
axis([1e-1 1e3 0 340])
%set(gca, 'xtick', [1e-5 1e-4 1e-3])
set(gca, 'xminortick', 'off')
xlabel('F_{DO} (\mu mol cm^{-2} yr^{-1})', 'FontSize', 9)
set(gca, 'FontSize', 9)
set(gca, 'yticklabel', [])
set(gca, 'xminortick', 'on')
set(gca, 'yminortick', 'on')
set(gca, 'tickdir', 'out')
set(gca, 'ygrid', 'on')

% Write letter identification
Ylim = [0 340];
text(5e2, 310, 'c', ...
         'verticalalignment', 'bottom', 'horizontalalignment', 'right', 'fontsize', 9, 'fontweight','bold');

adjust_space




% $$$ % ------ bin after this ------- %
% $$$ zbin=5; P_bin2 = [zmin+zbin/2:zbin:zmax-zbin/2]';
% $$$ for i = 1:length(P_bin2)
% $$$ I = find(P_bin >= P_bin2(i)-zbin/2 & P_bin <= P_bin2(i)+zbin/2); 
% $$$ KK(i) = nanmean(K_ave(I));
% $$$ end
% $$$ figure
% $$$ semilogx(KK, P_bin2)       
% $$$ set(gca, 'ydir', 'reverse')
% $$$ 
% $$$ disp('end')
% $$$ keyboard
% $$$ % ------------------------------- %



% save figure
disp('Saving figure...')
print('-dpng', '-r300','OKF_boot_5mbinCI.png')
set(gcf, 'renderer', 'painters')
print('-depsc2', 'OKF_boot_5mbinCI.eps')


% $$$ 
% $$$ % -------------- In the bottom layer -------------- %
% $$$ 
% $$$ zmin_cil = 100;
% $$$ zmax_cil = 320;
% $$$ I = find(P_bin>zmin_cil & P_bin<zmax_cil);
% $$$ 
% $$$ disp('\epsilon_{200-300m}...')
% $$$ % $$$ x1 = (eps_ave+eps_error);
% $$$ % $$$ x2 = (eps_ave-eps_error);
% $$$ x2 = eps_2p5;
% $$$ x1 = eps_97p5;
% $$$ [nanmean(eps_ave(I)) nanmean(x2(I)) nanmean(x1(I))]
% $$$ 
% $$$ disp('K_{200-300m}...')
% $$$ % $$$ x1 = (K_ave+K_error);
% $$$ % $$$ x2 = (K_ave-K_error);
% $$$ x2 = K_2p5;
% $$$ x1 = K_97p5;
% $$$ [nanmean(K_ave(I)) nanmean(x2(I)) nanmean(x1(I))]
% $$$ % -------------------------------------------------- %
% $$$ 
% $$$ % ---------  O2 flux (constant Kz) ----------------- %
% $$$ DO = O2_mean_bin*44660/1.42903/1028; % conversion to umol/Kg
% $$$ dodz = diff(DO)./dz;
% $$$ Fo2 = -nanmean(K_ave(I)).*dodz; %(in umol*m/kg/s)
% $$$ Fo2 = Fo2*3.2419e5; % conversion to umol/cm2/yr
% $$$ I = find(P_itp>=100);
% $$$ % -------------------------------------------------- %
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ subplot(1, 3, 3)
% $$$ 
% $$$ semilogx(Fo2(I), P_itp(I), 'k', 'linewidth', 0.25)
% $$$ hold on
% $$$ plot([0 0], [0 320], '--k','linewidth', 0.25)
% $$$ %semilogx(K_ave*0+5.5e-5, P_bin, '--k', 'linewidth', 0.25)
% $$$ %x2 = K_2p5;
% $$$ %x1 = K_97p5;
% $$$ 
% $$$ %patch([x1; flipud(x2); x1(1)], [P_bin; flipud(P_bin); P_bin(1)],
% $$$ %[.8 .8 .8], 'edgecolor', 'none');
% $$$ hold off
% $$$ set(gca, 'ydir', 'reverse')
% $$$ %axis([5e-6 5e-4 0 320])
% $$$ axis([1e0 1e2 0 320])
% $$$ %set(gca, 'xtick', [1e-5 1e-4 1e-3])
% $$$ set(gca, 'xminortick', 'off')
% $$$ xlabel('F_{DO} (\mu mol cm^{-2} yr^{-1})', 'FontSize', 9)
% $$$ set(gca, 'FontSize', 9)
% $$$ set(gca, 'yticklabel', [])
% $$$ set(gca, 'xminortick', 'on')
% $$$ set(gca, 'yminortick', 'on')
% $$$ set(gca, 'tickdir', 'out')
% $$$ set(gca, 'ygrid', 'on')
% $$$ 
% $$$ % Write letter identification
% $$$ Ylim = [0 320
% $$$ text(4e-5, 310, 'c', ...
% $$$          'verticalalignment', 'bottom', 'horizontalalignment', 'right', 'fontsize', 9, 'fontweight','bold');
% $$$ 
% $$$ 
% $$$ % save figure
% $$$ disp('K = cst; press any key to continue')
% $$$ pause
% $$$ disp('Saving figure...')
% $$$ print('-dpng', '-r300','OKF_boot_5mbinCI_Kcst.png')
% $$$ set(gcf, 'renderer', 'painters')
% $$$ print('-depsc2', 'OKF_boot_5mbinCI_Kcst.eps')


%%%%%%%%%%%%%%%%%%%%%%
% - save variables - %
%%%%%%%%%%%%%%%%%%%%%%
keyboard
if sav == 1;
    dlmwrite('boot_OK_5mbin.dat', [P_bin DO DO_2p5 DO_97p5 K_ave K_2p5 K_97p5], 'delimiter',' ','precision',6);
    dlmwrite('boot_F02_5mbin.dat', [P_itp Fo2], 'delimiter',' ','precision',6);                                 
    save boot_OKF.mat P_bin DO DO_2p5 DO_97p5 K_ave K_2p5 K_97p5 P_itp Fo2
end    

disp(sprintf('  %d profiles used', N))


