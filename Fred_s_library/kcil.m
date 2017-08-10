function kcil(zbin, hb, CIL, noload, kmax, transect, eps_list, gps_file, varargin)


% usage ex: >> kcil(5, 20, [30 150], 1, 1e-1, 'IML4transect1500.dat','eps_b19','gps_b19','eps_b20','gps_b20','eps_g19','gps_g19','eps_g20','gps_g20')   
% usage full ex: >> kcil(5, 10, [50 150], 1, 1,'IML4transect1500.dat','eps_b12','gps_b12','eps_b13','gps_b13', 'eps_b15','gps_b15','eps_b16','gps_b16', 'eps_b19','gps_b19','eps_b20','gps_b20','eps_b21','gps_b21','eps_b22','gps_b22','eps_b23','gps_b23','eps_g12','gps_g12','eps_g13','gps_g13', 'eps_g14','gps_g14','eps_g16','gps_g16', 'eps_g19','gps_g19','eps_g20','gps_g20','eps_g22','gps_g22','eps_g23','gps_g23' )   
%
%ex2: kcil(5, 10, [30 150], 1, 1e-1,'IML4transect1500.dat','eps_b12','gps_b12')
%    
% Each extra input must have profile_list AND GPS file in order.
% listing example:
% ls -1 eps_profile_b19_* | sed 's/\.mat//' > eps_b19  
% ls -1 eps_profile_b20_* | sed 's/\.mat//' > eps_b20
% ls -1 eps_profile_g20_* | sed 's/\.mat//' > eps_g20
% ls -1 eps_profile_g19_* | sed 's/\.mat//' > eps_g19

% NOTE: only take in account s1 !!!    
% 
% noload == 1 : files already saved, no need to compute everything
% (no need to give varargin)
% noload == 0 : Need to compute the whole matrix
% Ex: use noload==0 for the first run, then you can easily extract
% statistics by doing noload==1;
% kmax is the maximum K value admissible (removed if greater)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- First part, compute matrices --- %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if noload == 0

    
%% -- preamble -- %%
zmin = 0;
zmax = 200;
count=1;
gamma = 0.2;
g = 9.8;
rho_0 = 1025;
    
%% -- load files -- %%
LLDZ = load(transect); %[ln' lt' dv zv];
                       
%topo = load('/home/cyrf0006/data/Topo_Gulf_St_Lawrence.mat');


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
                
                if LLDZ(x_ind,4) <4
                    disp(sprintf('depth = %d in ',LLDZ(x_ind,4))); eps_files(prof,:)
                    
                    %   keyboard
                end
                count = count+1;
            end %for prof   
            
            %% -- Call function to have T,S,R -- %%
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

% sorting profiles
% $$$ [x_prof, I] = unique(x_prof);
% $$$ z_prof = z_prof(I);
% $$$ x_reg = floor(min(x_prof)):500/1000:ceil(max(x_prof));
% $$$ z_reg = interp1(x_prof, z_prof, x_reg);
% $$$ eps_bin_mat = eps_bin_mat(:,I);
% $$$ T_bin_mat = T_bin_mat(:,I);
% $$$ S_bin_mat = S_bin_mat(:,I);
% $$$ R_bin_mat = R_bin_mat(:,I);

%% -- Compute N2 and K -- %%
R_bin_mat = sort(R_bin_mat,1); % sorting density
N2_bin_mat = g*diff(R_bin_mat)/zbin/rho_0;
K_bin_mat = gamma*eps_bin_mat(1:end-1, :)./N2_bin_mat;
% $$$ 
% $$$ 
% $$$ %% -- Interpolation of T -- %%
% $$$ disp('interpolation...')
% $$$ T_mat_itp = nan(size(T_bin_mat,1),length(x_reg));
% $$$ for i = 1:size(T_bin_mat,1)
% $$$     I = find(~isnan(T_bin_mat(i,:))==1);
% $$$     if ~isempty(I)==1 & length(I)>1
% $$$         T_mat_itp(i,:) = interp1(x_prof(I), T_bin_mat(i,I), x_reg);
% $$$     end
% $$$ end
    
dlmwrite('k_matrix_all.dat', K_bin_mat,'delimiter',' ','precision',6);
dlmwrite('eps_matrix_all.dat', eps_bin_mat,'delimiter',' ','precision',6);
dlmwrite('T_matrix_all.dat', T_bin_mat,'delimiter',' ','precision',6);
dlmwrite('xprof_vector_all.dat', x_prof,'delimiter',' ','precision',6); % position relative to transect
dlmwrite('zprof_vector_all.dat', z_prof,'delimiter',' ','precision',6); % depth at transect position
dlmwrite('depth_vector_all.dat', p_bin,'delimiter',' ','precision',6); % depth vector


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ----- 2nd part, only statistics ----- %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif noload == 1
        
   
    % must use a better bathymetry
    Kmat = load('k_matrix_all.dat');
    Emat = load('eps_matrix_all.dat');
    bottom = load('zprof_vector_all.dat');
    z = load('depth_vector_all.dat');
    count = 1;
    count2 = 1;

    
    for i=1:size(Kmat,2)
        
        I = find(~isnan(Kmat(:,i))==1);
        Kprof = Kmat(I,i); % single Kprofile
        Eprof = Emat(I,i); % single Kprofile

        zz = z(I); % with his depth vector
        keyboard
        % Mean K boundary
        if (max(zz) < CIL(1)) | (max(zz) > CIL(2))
            continue % profile not in border
        else           
            J = find(zz > max(zz)-hb); % z in bndry layer  
            Kbndry(count: count+length(J)-1) = Kprof(J);
            Ebndry(count: count+length(J)-1) = Eprof(J);
            count = count+length(J);
        end
        
        % MEan K out of boundary
        if max(zz) < CIL(1)+hb % profile too short
            continue
        elseif max(zz) > CIL(2)+hb % deep profile away from boundary
            J = find(zz > CIL(1) & zz < CIL(2)); % z in bndry layer  
            Kcenter(count: count+length(J)-1) = Kprof(J);
            Ecenter(count: count+length(J)-1) = Eprof(J);
            count2 = count2+length(J);
        else % profile finishes in the CIL
            J = find(zz > CIL(1) & zz<max(zz)-hb); % z in bndry layer  
            Kcenter(count: count+length(J)-1) = Kprof(J);
            Ecenter(count: count+length(J)-1) = Eprof(J);
            count2 = count2+length(J);
        end
        
        
        
    end
    
    
else
    disp('wrong input')
end
    
 
I = find(Kbndry>kmax);
Kbndry(I)=[];

bdry = [mean(Kbndry) mean(Ebndry)]
center = [mean(Kcenter) mean(Ecenter)]

%keyboard

% bootstrap
nboot=1000;
N = length(Kbndry);
% create random sampling
for b = 1:nboot
    r = rand(N,1);
    r = ceil(r*N/1);
    
    K_boot_b(b) = nanmean(Kbndry(r));
end

% mean of random sampling
K_boot_dot = nanmean(K_boot_b);

% Compute 95% confidence interval
K_sort  = sort(K_boot_b);

CI_2p5 = round(2.5/100*nboot);
CI_97p5 = round(97.5/100*nboot);

K_2p5 = K_sort(CI_2p5);
K_97p5 = K_sort(CI_97p5);

keyboard
[K_boot_dot K_2p5 K_97p5]