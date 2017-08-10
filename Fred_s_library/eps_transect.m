function eps_transect(zbin, transect, eps_list, gps_file, varargin)


% usage ex: eps_transect(5, 'IML4_transect', 'eps_b19','gps_b19','eps_b20','gps_b20', ...)
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
        end %for j
    end %modulo
end %isempty


% sorting profiles
[x_prof, I] = sort(x_prof);
z_prof = z_prof(I);
eps_bin_mat = eps_bin_mat(:,I);


%% -- PLotting section -- %%

figure(1)
clf
set(gcf,'PaperUnits','centimeters','PaperPosition',[1 1 15 15])


% eps_bin_mat, x_prof, z_prof, p_bin
epsMin = 1e-10;
epsMax = 1e-5;
% dissipation color discretization
eps_cat = log10(epsMax/epsMin); %number of color categories for epsilon

%log scale
for j = 1:eps_cat %
    iexp = -log10(epsMin)-j+1;
    eps_value(j) = 10^(-(iexp));
end

% set the colorbar
mycolormap = colormap(jet(length(eps_value))); %matrix of RGB colors
fst_plot=0;


for i = 1:size(eps_bin_mat,1) %40
    i
    for j = 1:size(eps_bin_mat,2) %135
    
        if isnan(eps_bin_mat(i, j))==1
            continue
        else
        
            eps_color_ind = dsearchn(eps_value', eps_bin_mat(i,j)); %indice of the epsilon colorbar
            RGB_color = mycolormap(eps_color_ind, :);
            h=plot(x_prof(j),p_bin(i), 'o');
            set(h, 'MarkerFaceColor', RGB_color,'MarkerEdgeColor', RGB_color)
            axis([0 50 0 200])

            
            if fst_plot==0
                hold on
                fst_plot==1;
            end
        
        end
        
    end
end
hold off
%keyboard
% colorbar label for epsilon must be no_categories+1
eps_label = [eps_value eps_value(eps_cat)*10]; %add the last colorbar label

set(gca, 'ydir', 'reverse')

c = colorbar('YTickLabel', {eps_label}, 'FontSize', 14);
ti = ylabel(c,'{log \epsilon (W/Kg)}', 'FontSize', 14);


ylabel(c,'{log \epsilon (W/Kg)}', 'FontSize', 14)
%title({'Taux de dissipation de TKE'; datestr(n,1)},'FontSize', 14)
xlabel('x (km)','FontSize', 14)
ylabel('depth(m)', 'FontSize', 14)

set(gca, 'FontSize', 14)

%move colorbar and its title
pos_ti = get(ti,'position'); %it gives a position of 0.0500 2.900 1.0001
%pos_c = get(c,'position'); %it gives a position of 0.0500 2.900 1.0001
%set(c,'position', [pos_c(1) pos_c(2) pos_c(3) pos_c(4)]);
set(ti,'position', [0 pos_ti(2) pos_ti(3)]);



hold on
plot(x_prof, z_prof+20, 'k', 'linewidth', 2)
hold off

keyboard

print('-dpng', '-r300', 'eps_transect.png');  
print('-depsc2', 'eps_transect.eps');                                    

