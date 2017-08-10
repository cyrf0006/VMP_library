
function T_transect(zbin, transect, p_list, gps_file, varargin)


% usage ex: T_transect(5, 'IML4transect.dat', 'p_b19','gps_b19','p_b20','gps_b20', 'p_g19','gps_g19','p_g20','gps_g20') 
%
% Each extra input must have profile_list AND GPS file in order.
% listing example:
% ls -1 eps_profile_b19_* | sed 's/\.mat//' > eps_b19  
% ls -1 eps_profile_b20_* | sed 's/\.mat//' > eps_b20
% ls -1 eps_profile_g20_* | sed 's/\.mat//' > eps_g20
% ls -1 eps_profile_g19_* | sed 's/\.mat//' > eps_g19

% NOTE: only take in account s1 !!!    
    
%% -- preamble -- %%
zmin = 0;
zmax = 200;
count=1;
    
%% -- load files -- %%
LLDZ = load(transect); %[ln' lt' dv zv];


% load eps_profiles* names 
fid = fopen(p_list);
C = textscan(fid, '%s', 'delimiter', '\n');
p_files = char(C{1});

% xtract GPS infos
[lat_gps lon_gps n_gps] = xtract_track(gps_file);
lon_gps = lon_gps*-1;

%% -- Compute 1st profile -- %%
% load first file to compute pressure vector
load(p_files(1,:));
                          
p_bin = [zmin:zbin:zmax]';

% T binned
for i = 1:length(p_bin)
    I = find(P>(zmin+(i-1)*zbin) & P<(zmin+(i)*zbin));
    T_bin(i) = nanmean(SBT(I));
end


% empty epsilon matrix
T_bin_mat = nan(length(p_bin),size(p_files,1));

% ranging 1st profile
T_bin_mat(1:length(p_bin),1) = T_bin;

% profile coord.
[mi n_ind] = min(abs(n_gps-mtime(1)));
lat_prof(count) = lat_gps(n_ind);
lon_prof(count) = lon_gps(n_ind);
    
% faut minimiser lat lon avec x    
[mi x_ind] = min(abs(lat_prof-LLDZ(:,2))+abs(lon_prof-LLDZ(:,1)));    
x_prof(count) = LLDZ(x_ind,3);    
z_prof(count) = LLDZ(x_ind,4);


count = count+1;

% loop on epsilon profiles
for prof = 2:size(p_files,1)
    
    load(p_files(prof,:));

    clear T_bin
    % T binned
    for i = 1:length(p_bin)
        I = find(P>(zmin+(i-1)*zbin) & P<(zmin+(i)*zbin));
        T_bin(i) = nanmean(SBT(I));
    end

    % ranging profiles
    T_bin_mat(1:length(p_bin),count) = T_bin;

    % profile coord.
    [mi n_ind] = min(abs(n_gps-mtime(1)));
    lat_prof = lat_gps(n_ind);
    lon_prof = lon_gps(n_ind);
    
    % faut minimiser lat lon avec x    
    [mi x_ind] = min(abs(lat_prof-LLDZ(:,2))+abs(lon_prof-LLDZ(:,1)));    
    x_prof(count) = LLDZ(x_ind,3);    
    z_prof(count) = LLDZ(x_ind,4);


    count = count+1;

    
end


%% -- if varargin is not empty (other epsilon files) -- %%

% -> same shit as before! (except no need to compute press. vector)

if ~isempty(varargin)==1 % varargin not empty
    
    
    if mod(size(varargin,2),2)~=0
        disp('Problem with input files, please provide epsilon_fileslist AND gps_tracks')
        return
    else

        for j = 1:2:size(varargin,2)-1
            
            % new file names
            p_list = varargin{j};
            gps_file = varargin{j+1};
            
            % load eps_profiles* names 
            fid = fopen(p_list);
            C = textscan(fid, '%s', 'delimiter', '\n');
            p_files = char(C{1});

            % xtract GPS infos
            [lat_gps lon_gps n_gps] = xtract_track(gps_file);
            lon_gps = lon_gps*-1;


            % loop on epsilon profiles
            for prof = 1:size(p_files,1)
                load(p_files(prof,:));

                clear T_bin
                % T binned
                for i = 1:length(p_bin)
                    I = find(P>(zmin+(i-1)*zbin) & P<(zmin+(i)*zbin));
                    T_bin(i) = nanmean(SBT(I));
                end

                % ranging 1st profile
                T_bin_mat(1:length(p_bin),count) = T_bin;

                % profile coord.
                [mi n_ind] = min(abs(n_gps-mtime(1)));
                lat_prof = lat_gps(n_ind);
                lon_prof = lon_gps(n_ind);

                % faut minimiser lat lon avec x    
                [mi x_ind] = min(abs(lat_prof-LLDZ(:,2))+abs(lon_prof-LLDZ(:,1)));    
                x_prof(count) = LLDZ(x_ind,3);    
                z_prof(count) = LLDZ(x_ind,4);


                count = count+1;       
                
            end %for prof   
        end %for j
    end %modulo
end %isempty


% sorting profiles (1st version, ignore repetition)
[x_prof, I] = unique(x_prof);
z_prof = z_prof(I);
T_bin_mat = T_bin_mat(:,I);
x_reg = round(min(x_prof)):500/1000:round(max(x_prof));


%% -- Interpolation -- %%
disp('interpolation...')
T_mat_itp = nan(size(T_bin_mat,1),length(x_reg));
for i = 1:size(T_bin_mat,1)
     i
    I = find(~isnan(T_bin_mat(i,:))==1);
    if ~isempty(I)==1
        T_mat_itp(i,:) = interp1(x_prof(I), T_bin_mat(i,I), x_reg);
    end
end


%% -- PLotting section -- %%
disp('plotting...')
load BR_colormap;

figure(1)
clf
%whitebg('black')
set(gcf,'PaperUnits','centimeters','PaperPosition',[1 1 15 15])

contourf(x_reg, p_bin, T_mat_itp, [-1:0.1:5], 'linestyle', 'none')
hold on
contour(x_reg, p_bin, T_mat_itp, [1 1 1], 'edgecolor', [0 0 0], 'LineStyle', '-','linewidth', 0.25 )
hold off 

set(gca, 'ydir', 'reverse')
set(gca, 'xticklabel', [], 'fontsize', 10)
set(gca, 'xtick', TIK)
set(gca, 'fontsize', 10)
set(gca, 'XGrid', 'on')
set(gca, 'xticklabel', [])
set(gca, 'YGrid', 'on')

% colorbar and its title
caxis([0 4])
%colormap(BR_colormap)
colormap(flipud(gray))
colorbar

keyboard


print('-dpng', '-r300', 'T_transect.png');  
set(gcf, 'renderer', 'painters')
print('-depsc2', 'T_transect.eps');                                    

