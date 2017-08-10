function [T_bin_mat, S_bin_mat, R_bin_mat]=teps_profile(p_bin, transect, p_list, gps_file)


% function that returns t,s,rho for a certain eps_profile_list
% certainly not efficient in term of calculation!!
% $$$    
% $$$ %% -- preamble -- %%
zmin = min(p_bin);
zmax = max(p_bin);
zbin = p_bin(2)-p_bin(1);
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
                          
% $$$ p_bin = [zmin:zbin:zmax]';

% T binned
for i = 1:length(p_bin)
    I = find(P>(zmin+(i-1)*zbin) & P<(zmin+(i)*zbin));
    T_bin(i) = nanmean(SBT(I));
    S_bin(i) = nanmean(SBS(I));
    R_bin(i) = nanmean(sw_dens(SBS(I), SBT(I), SBT(I)*0));
end


% empty epsilon matrix
T_bin_mat = nan(length(p_bin),size(p_files,1));
S_bin_mat = nan(length(p_bin),size(p_files,1));
R_bin_mat = nan(length(p_bin),size(p_files,1));

% ranging 1st profile
T_bin_mat(1:length(p_bin),1) = T_bin;
S_bin_mat(1:length(p_bin),1) = S_bin;
R_bin_mat(1:length(p_bin),1) = R_bin;

count = count+1;

% loop on epsilon profiles
for prof = 2:size(p_files,1)
    load(p_files(prof,:));
    
    clear T_bin S_bin R_bin
    % T binned
    for i = 1:length(p_bin)
        I = find(P>(zmin+(i-1)*zbin) & P<(zmin+(i)*zbin));
        T_bin(i) = nanmean(SBT(I));
        S_bin(i) = nanmean(SBS(I));
        R_bin(i) = nanmean(sw_dens(SBS(I), SBT(I), SBT(I)*0));        
    end

    % ranging profiles
    T_bin_mat(1:length(p_bin),count) = T_bin;
    S_bin_mat(1:length(p_bin),count) = S_bin;
    R_bin_mat(1:length(p_bin),count) = R_bin;

    count = count+1;

    
end
