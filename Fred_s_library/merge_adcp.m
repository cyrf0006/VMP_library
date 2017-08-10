function merge_adcp(adcp_files, timestep, zbin)
 
% function merge_adcp(adcp_files, timestep, zbin, merged_file)
%                       (list)      (sec.)   (m)  (output file)
%
%ex: >> merge_adcp('adcpfiles', 120, 4)
%     >> merge_adcp('adcpfilesg', 1200, 4)
%
%
% Will compute shear from ADCP files and will merge these file to a
% single regular Z grid file with time vector.
%
% - timestep of the average in sec.
%
% Please have a look at time correction during July 2010 mission

% default parameters
z0 = 1; %m default ADCP depth if no pressure sensor    
mag_var = -18.5;

    
% depth range for extracting the average
zmin = 0; % top of the fisrt bin
zmax = 150;
nboot = 500;    
P_bin = [zmin+zbin/2:zbin:zmax-zbin/2]';
P_S2 = P_bin(1:end-1)+zbin;

% load *.P files names (file in which are recorded .P files)
fid = fopen(adcp_files);
C = textscan(fid, '%s', 'delimiter', '\n');
files = char(C{1});
no_files = size(files, 1); %number of eps_files 




% raw matrix to fill
S2 = nan(length(P_S2), 1); % start with one profile
count = 1;

%%%%%%%%%%%%%%%%%%%%
% loop on profiles %
%%%%%%%%%%%%%%%%%%%%

for i = 1:no_files
    
    % ---- load ADCP file ---- %
    fname = files(i, :);
    I = find(fname==' ');   
    fname(I) = [];
    disp(fname);
    load(fname);
    % ----------------------- %    

    % ---- Extract info ---- %
    I = find(isnan(nansum(ADCP.east_vel))~=1);
    I(end) = []; % just to make sure!
    %I = find(ADCP.depth > 0.25); % ADCP in water
    zadcp = nanmean(ADCP.depth(I));
    if zadcp == 0
        zadcp = z0;
    end
    z = ADCP.config.ranges + zadcp;
    Evel = ADCP.east_vel(:,I);
    Nvel = ADCP.north_vel(:,I);
    timeadcp = ADCP.mtime(I);
    timevec = timeadcp(1):timestep/86400:timeadcp(end);
    if zadcp < 0 | zadcp > 10
        disp('problem with pressure sensor?')
        zadcp
        keyboard
    end 
    % ---------------------- %
 
     % ---- check few things ---- %
    % Correct for magvar
% $$$     if ADCP.magnetic_var < -20 | ADCP.magnetic_var > -18
% $$$         [Evel, Nvel] = rotate_vecd(Evel, Nvel, mag_var)
% $$$     end
    % --------------------------- %
    
    % ---- Bin time and depth ---- %
    % bin Z
    U = nan(length(P_bin), length(timeadcp));
    V = nan(length(P_bin), length(timeadcp));
    for j = 1:length(P_bin)
        I = find(z >= P_bin(j)-zbin/2 & z <= P_bin(j)+zbin/2);
        if ~isempty(I) == 1
            U(j, :) = nanmean(Evel(I, :));
            V(j, :) = nanmean(Nvel(I, :));
        end
    end
    clear Evel Nvel z

   
    % bin T 
    UU = nan(length(P_bin), length(timevec));
    VV = nan(length(P_bin), length(timevec));
    for j = 2:length(timevec)-1
        I = find(timeadcp >= timevec(j)-timestep/2/86400 & timeadcp <= timevec(j)+timestep/2/86400);
        if ~isempty(I) == 1
            UU(:, j) = nanmean(U(:, I), 2);
            VV(:, j) = nanmean(U(:, I), 2);
        end
    end
    clear U V   
    % ------------------------------ %
    
    dudz = diff(UU, 1)./zbin;
    dvdz = diff(VV, 1)./zbin;
    s2 = dudz.^2 + dvdz.^2;
    P_S2 = P_bin(1:end-1)+zbin;
    
    
    S2(:, count:count+length(timevec)-1) = s2;
    t_S2(count:count+length(timevec)-1) = timevec;
    count = count+length(timevec);
    
end


save merged_adcp_shear S2 t_S2 P_S2;