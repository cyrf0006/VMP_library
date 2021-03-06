function R_transect(zbin, transect, eps_list, gps_file, varargin)


% usage ex: >> R_transect(5, 'IML4transect1500.dat', 'eps_t_j10-b1','gps_b23', 'eps_t_j10-g1','gps_g23')
% usage ex: >> R_transect(5, 'IML4transect1500.dat', 'eps_t_j07-b1','gps_b20', 'eps_t_j07-g1','gps_g20')
% usage ex: >> R_transect(5, 'IML4transect1500.dat', 'eps_t_j08-b1','gps_b21')
% usage ex: >> R_transect(5, 'IML4transect1500.dat', 'eps_t_j08-b2','gps_b21')
% usage ex: >> R_transect(5, 'IML4transect1500.dat', 'eps_t_j08-b3','gps_b21')
% usage ex: >> R_transect(5, 'IML4transect1500.dat', 'eps_t_j05-b1','gps_b16', 'eps_t_j07-g1','gps_g16')
%

%% -- preamble -- %%
zmin = 0;
zmax = 200;
count=1;
gamma = 0.2;
g = 9.8;
rho_0 = 1025;
    
%% -- load files -- %%

%LLDZ = load(transect); %[ln' lt' dv zv];
% By-pass this...
bath = load('/home/cyrf0006/WINDEX/data_processing/MStage_2011/Transect_Simrad/bathym_latlon.mat');
ln = bath.lon';
lt = bath.lat';
dv = bath.x';
zv = bath.H';
LLDZ = [ln lt dv zv];


% load eps_profiles* names 
fid = fopen(eps_list);
C = textscan(fid, '%s', 'delimiter', '\n');
eps_files = char(C{1});

% xtract GPS infos
[lat_gps lon_gps n_gps] = xtract_track(gps_file);
lon_gps = lon_gps*-1;

%% -- Compute 1st profile -- %%
% load first file to compute pressure vector
load(eps_files(1,:));

dp = p_eps1(3)-p_eps1(1); % skip one value out of 2
z1 = zmin+dp/2;                          

p = p_eps1;
eps = eps1;
p_bin = [z1:zbin:zmax]';

% epsilon binned
for i = 1:length(p_bin)
    I = find(p>(zmin+(i-1)*zbin) & p<(zmin+(i)*zbin));
    eps_bin(i) = nanmean(eps(I));
end


% empty epsilon matrix
eps_bin_mat = nan(length(p_bin),size(eps_files,1));
T_bin_mat = eps_bin_mat;
S_bin_mat = eps_bin_mat;
R_bin_mat = eps_bin_mat;

% ranging 1st profile
eps_bin_mat(1:length(p_bin),1) = eps_bin;

% profile coord.
[mi n_ind] = min(abs(n_gps-mtime_eps(1)));
lat_prof = lat_gps(n_ind);
lon_prof = lon_gps(n_ind);
    
% faut minimiser lat lon avec x    
[mi x_ind] = min(abs(lat_prof-LLDZ(:,2))+abs(lon_prof-LLDZ(:,1)));    
x_prof(count) = LLDZ(x_ind,3);    
z_prof(count) = LLDZ(x_ind,4);


count = count+1;

% loop on epsilon profiles
for prof = 2:size(eps_files,1)
    
    load(eps_files(prof,:));

    p = p_eps1;
    eps = eps1;

    clear eps_bin
    % epsilon binned
    for i = 1:length(p_bin)
        I = find(p>(zmin+(i-1)*zbin) & p<(zmin+(i)*zbin));
        eps_bin(i) = nanmean(eps(I));
    end

    % ranging profile
    eps_bin_mat(1:length(p_bin),count) = eps_bin; 
   
    % profile coord.
    [mi n_ind] = min(abs(n_gps-mtime_eps(1)));
    lat_prof = lat_gps(n_ind);
    lon_prof = lon_gps(n_ind);

    % faut minimiser lat lon avec x
    [mi x_ind] = min(abs(lat_prof-LLDZ(:,2))+abs(lon_prof-LLDZ(:,1)));
    x_prof(count) = LLDZ(x_ind,3);
    z_prof(count) = LLDZ(x_ind,4);
    
    count = count+1;
    
end

           
%% -- Call function to have T,S,R -- %%
base = eps_list;
base(1:3)=[];
p_list = ['p', base]; 
[T, S, R] = teps_profile(p_bin, transect, p_list, gps_file);

T_bin_mat = T;
S_bin_mat = S;
R_bin_mat = R;


%% -- if varargin is not empty (other epsilon files) -- %%

% -> same shit as before! (except no need to compute press. vector)

if ~isempty(varargin)==1 % varargin not empty
    
    
    if mod(size(varargin,2),2)~=0
        disp('Problem with input files, please provide epsilon_fileslist AND gps_tracks')
        return
    else

        for j = 1:2:size(varargin,2)-1
            
            
            % new file names
            eps_list = varargin{j};
            gps_file = varargin{j+1};
          
         
            % load eps_profiles* names 
            fid = fopen(eps_list);
            C = textscan(fid, '%s', 'delimiter', '\n');
            eps_files = char(C{1});

            % xtract GPS infos
            [lat_gps lon_gps n_gps] = xtract_track(gps_file);
            lon_gps = lon_gps*-1;


            % loop on epsilon profiles
            for prof = 1:size(eps_files,1)
    
                load(eps_files(prof,:));

                p = p_eps1;
                eps = eps1;

                clear eps_bin
                % epsilon binned
                for i = 1:length(p_bin)
                    I = find(p>(zmin+(i-1)*zbin) & p<(zmin+(i)*zbin));
                    eps_bin(i) = nanmean(eps(I));
                end

                % ranging profile
                eps_bin_mat(1:length(p_bin),count) = eps_bin; 

                % profile coord.
                [mi n_ind] = min(abs(n_gps-mtime_eps(1)));
                lat_prof = lat_gps(n_ind);
                lon_prof = lon_gps(n_ind);

                % faut minimiser lat lon avec x
                [mi x_ind] = min(abs(lat_prof-LLDZ(:,2))+abs(lon_prof-LLDZ(:,1)));
                
                x_prof(count) = LLDZ(x_ind,3);
                z_prof(count) = LLDZ(x_ind,4);
                
                count = count+1;
            end %for prof   
            
            %% -- Call function to have T,S,R __ %%
            base = eps_list;
            base(1:3)=[];
            p_list = ['p', base]; 

            [T, S, R] = teps_profile(p_bin, transect, p_list, gps_file);
            
            T_bin_mat = [T_bin_mat T];
            S_bin_mat = [S_bin_mat S];
            R_bin_mat = [R_bin_mat R];
            
            
        end %for j
    end %modulo
end %isempty

% sorting profiles (2nd version, average when repetition)
[x_prof_unique, I] = unique(x_prof);
z_prof = z_prof(I);

for i=1:length(I)
    II = find(x_prof==x_prof_unique(i));
    eps_bin_matu(:,i) = nanmean(eps_bin_mat(:,II), 2);
    T_bin_matu(:,i) = nanmean(T_bin_mat(:,II), 2);
    S_bin_matu(:,i) = nanmean(S_bin_mat(:,II), 2);
    R_bin_matu(:,i) = nanmean(R_bin_mat(:,II), 2);
end

eps_bin_mat = eps_bin_matu;
T_bin_mat = T_bin_matu;
S_bin_mat = S_bin_matu;
R_bin_mat = sort(R_bin_matu, 1); % sorting rho
x_prof = x_prof_unique;

x_reg = floor(min(x_prof)):500/1000:ceil(max(x_prof));
z_reg = interp1(x_prof, z_prof, x_reg);

%% -- Compute N2 and K -- %%
N2_bin_mat = g*diff(R_bin_mat)/zbin/rho_0;
K_bin_mat = gamma*eps_bin_mat(1:end-1, :)./N2_bin_mat;


%% -- Interpolation of T -- %%
disp('interpolation...')
T_mat_itp = nan(size(T_bin_mat,1),length(x_reg));
for i = 1:size(T_bin_mat,1)
    I = find(~isnan(T_bin_mat(i,:))==1);
    if ~isempty(I)==1 & length(I)>1
        T_mat_itp(i,:) = interp1(x_prof(I), T_bin_mat(i,I), x_reg);
    end
end



%% -- Interpolation of Rho -- %%
R_mat_itp = nan(size(R_bin_mat,1),length(x_reg));
for i = 1:size(R_bin_mat,1)
    I = find(~isnan(R_bin_mat(i,:))==1);
    if ~isempty(I)==1 & length(I)>1
        R_mat_itp(i,:) = interp1(x_prof(I), R_bin_mat(i,I), x_reg);
    end
end

outfile = ['rho_' datestr(mtime_eps(1),1)];
save(outfile, 'R_mat_itp', 'x_reg', 'p_bin');


%% -- PLotting section -- %%
figure(1) % for T
clf
set(gcf,'PaperUnits','centimeters','PaperPosition',[1 1 15 15])

% - temperature plot - %
%load BR_colormap;

contourf(x_reg, p_bin, T_mat_itp, [-1:0.1:5], 'linestyle', 'none')
hold on
contour(x_reg, p_bin, T_mat_itp, [1 1 1], 'edgecolor', [0 0 0], 'LineStyle', '-','linewidth', 0.25 )
hold off 

set(gca, 'ydir', 'reverse')
%set(gca, 'xticklabel', [], 'fontsize', 10)
set(gca, 'fontsize', 10)
set(gca, 'XGrid', 'on')
%set(gca, 'xticklabel', [])
set(gca, 'YGrid', 'on')
xlabel('x (km)','FontSize', 12)
ylabel('depth(m)', 'FontSize', 12)
tit = p_list;
tit(1:2) = [];
tit(2) = [];
title(sprintf('Temperature for %s', tit))   

%axis([28 40 0 200]) %zoom
axis([0 45 0 200]) 

% colorbar and its title
caxis([-2 5])
%colormap(mycolormap)
%colormap(flipud(gray));

set(gca, 'xdir', 'normal')
colorbar

set(gcf, 'renderer', 'painters')
print('-depsc2', 'T_t_d.eps');      

% Find extras contour for rho_CIL
for i = 1:size(T_mat_itp,2)
    I = find(T_mat_itp(:,i)<=1);
      
    if isempty(I)==1
       cil_top(i) = NaN;
       cil_bot(i) = NaN;
    else        
        cil_top(i) = p_bin(I(1));
        cil_bot(i) = p_bin(I(end));
    end
end


figure(2) % for D
clf
set(gcf,'PaperUnits','centimeters','PaperPosition',[1 1 15 15])

% - temperature plot - %
%load BR_colormap;

contourf(x_reg, p_bin, R_mat_itp, [1018:0.1:1027], 'linestyle', 'none')
hold on
plot(x_reg, cil_top, 'k');
plot(x_reg, cil_bot, 'k');
%ontour(x_reg, p_bin, R_mat_itp, [rho_cil rho_cil rho_cil], 'edgecolor', [0 0 0], 'LineStyle', '-','linewidth', 0.25 )
hold off 

set(gca, 'ydir', 'reverse')
%set(gca, 'xticklabel', [], 'fontsize', 10)
set(gca, 'fontsize', 10)
set(gca, 'XGrid', 'on')
%set(gca, 'xticklabel', [])
set(gca, 'YGrid', 'on')
xlabel('x (km)','FontSize', 12)
ylabel('depth(m)', 'FontSize', 12)
title(sprintf('Density for %s', tit))   

%axis([28 40 0 200]) %zoom
axis([0 45 0 200]) 




% colorbar and its title
%colormap(mycolormap)
%colormap(flipud(gray));
set(gca, 'xdir', 'normal')
colorbar
caxis([1022 1027])


set(gcf, 'renderer', 'painters')
print('-depsc2', 'R_t_d.eps');                                    



% $$$ 
% $$$ freezeColors
% $$$ 
% $$$ 
% $$$ % $$$ %set(gcf, 'renderer', 'painters')
% $$$ % $$$ print('-dpng', '-r300', 'TK_transect.png');  
% $$$ % $$$ print('-depsc2', 'TK_transect.eps');                                    
% $$$ 
% $$$ % eps_bin_mat, x_prof, z_prof, p_bin
% $$$ epsMin = 1e-7;
% $$$ epsMax = 1e-3;
% $$$ % $$$ epsMin = 1e-8;
% $$$ % $$$ epsMax = 1e-5;
% $$$ 
% $$$ % $$$ % dissipation color discretization
% $$$ % $$$ eps_cat = log10(epsMax/epsMin); %number of color categories for epsilon
% $$$ % $$$ 
% $$$ % $$$ %log scale
% $$$ % $$$ for j = 1:eps_cat %
% $$$ % $$$     iexp = -log10(epsMin)-j+1;
% $$$ % $$$     eps_value(j) = 10^(-(iexp));
% $$$ % $$$ end
% $$$ % $$$ 
% $$$ % $$$ % set the colorbar
% $$$ % $$$ mycolormap = colormap(jet(length(eps_value))); %matrix of RGB colors
% $$$ % $$$                                                
% $$$ % $$$ %fstplot=0;
% $$$ 
% $$$ hold on
% $$$ for j = 1:size(K_bin_mat,2) %135
% $$$     j            
% $$$     
% $$$     [mi x_ind] = min(abs(x_prof(j)-x_reg)); % corresponding
% $$$                                             
% $$$     %h=plot(x_reg(x_ind),p_bin(i), 'o'); 
% $$$     %data = eps_bin_mat(:,j);
% $$$     data = K_bin_mat(:,j);
% $$$     x(1:length(data)) = x_reg(x_ind);
% $$$     y = p_bin;
% $$$     
% $$$     
% $$$     %rectbar = [0.9018 0.1100 0.0281 0.8150];
% $$$     rectbar = [0.92 0.15 0.02 0.7];
% $$$     %keyboard
% $$$     %cb = colour_profile(x',y',data,epsMin,epsMax,1,rectbar,x'.*0);
% $$$     
% $$$     if j == size(K_bin_mat,2) % with cbar
% $$$         cb = colour_profile(x',y',log10(data),log10(epsMin),log10(epsMax),1,rectbar, x'.*0);
% $$$     else % nocbar
% $$$         colour_profile(x',y',log10(data),log10(epsMin),log10(epsMax),0,rectbar, x'.*0);
% $$$     end
% $$$  
% $$$ end
% $$$ hold off
% $$$ 
% $$$ keyboard
% $$$ 
% $$$ % colorbar title
% $$$ ti = ylabel(cb,'log(K{m^2 s^{-1}})', 'FontSize', 10);
% $$$ cbpos = get(cb, 'position');
% $$$ set(gca, 'outer', [0 0 0.97 1])
% $$$ set(cb, 'pos', [cbpos(1)-0.03 cbpos(2) cbpos(3) cbpos(4)])
% $$$ 
% $$$ print('-dpng', '-r300', 'K_transect.png');  
% $$$ set(gcf, 'renderer', 'painters')
% $$$ print('-depsc2', 'K_transect.eps');    
% $$$ 
