function isriki(rad, outfile, file_list, gps_file, varargin)
    
% function list = isriki(rad, outfile, file_list, gps_file,varargin)
%
% List all profile that are within a distance rad from Rimouski
% station. Typically, if you give a profile list and GPS file, the
% script will return all profile that have been taken at Riki
% station.
%
% usage ex: >> isriki(5, 'at_riki', 'eps_b19','gps_b19','eps_b20','gps_b20', 'eps_g19','gps_g19','eps_g20','gps_g20')
%
% usage ex2: >> isriki(5, 'at_riki', 'eps_b12','gps_b12','eps_b13','gps_b13', 'eps_b15','gps_b15','eps_b16','gps_b16', 'eps_b19','gps_b19','eps_b20','gps_b20','eps_b21','gps_b21','eps_b22','gps_b22','eps_b23','gps_b23','eps_g12','gps_g12','eps_g13','gps_g13', 'eps_g14','gps_g14','eps_g16','gps_g16', 'eps_g19','gps_g19','eps_g20','gps_g20','eps_g22','gps_g22','eps_g23','gps_g23' )   
%
%  usage ex3: >> isriki(5, 'at_riki', 'prof_20110520', '20110520', 'prof_20110525', '20110525', 'prof_20110622', '20110622', 'prof_20110628', '20110628', 'prof_20110714', '20110714', 'prof_20110721','20110721', 'prof_20110809', '20110809', 'prof_20110920', '20110920', 'prof_20110922', '20110922', 'prof_20110923', '20110923', 'prof_20110926', '20110926', 'prof_20110927', '20110927', 'prof_20110928', '2110928', 'prof_20110929', '20110929', 'prof_20111013', '20111013', 'prof_20111021', '20111021', 'prof_20111109', '20111109')
% where: 
% - rad: the distance in km
% - outfile: the name of the list file
% - file_list: eps_profile* or profile* list
% - gps_file: corresponding gps file
% - varargin: any repetition of the last two input (see K_transect.m)
%
% Each extra input must have profile_list AND GPS file in order.
% listing example:
% ls -1 eps_profile_b19_* | sed 's/\.mat//' > eps_b19  
% ls -1 eps_profile_b20_* | sed 's/\.mat//' > eps_b20
% ls -1 eps_profile_g20_* | sed 's/\.mat//' > eps_g20
% ls -1 eps_profile_g19_* | sed 's/\.mat//' > eps_g19
%
% Frederic Cyr (2011-02-07)
%
% ----------------------------------------------------------- %


    
lat_riki = 48.666666;
lon_riki = -68.583333;
    
outfile_all = 'all_coordinates.dat';
    
% load eps_profiles* names 
fid = fopen(file_list);
C = textscan(fid, '%s', 'delimiter', '\n');
eps_files = char(C{1});

% loop on epsilon profiles
count = 1;
count2 = 1;
for prof = 1:size(eps_files,1)
    
    [lat_prof, lon_prof] = whereis(eps_files(prof,:), gps_file);
    
    if(isthere(lat_riki, lon_riki, lat_prof, lon_prof, rad)==1);
                
        if count == 1
            fid = fopen(outfile,'w+');
        else
            fid = fopen(outfile,'a');
        end
        % write profile to outfile 
        fprintf(fid,'%s\n.mat',eps_files(prof, :));
        fclose(fid);
        
        count = count+1;
    end
    
    % Write all corrdinates anyway
    LAT(count2) = lat_prof;
    LON(count2) = lon_prof;
    count2 = count2+1;
 
end


%% -- if varargin is not empty (other epsilon files) -- %%

% -> same shit as before! 

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

            for prof = 1:size(eps_files,1)
                
                [lat_prof, lon_prof] = whereis(eps_files(prof,:), gps_file);
                
                if(isthere(lat_riki, lon_riki, lat_prof, lon_prof, rad)==1);
                    
                    if count == 1
                        fid = fopen(outfile,'w+');
                    else
                        fid = fopen(outfile,'a');
                    end
                    % write profile to outfile 
                    fprintf(fid,'%s.mat\n',eps_files(prof, :));
                    fclose(fid);
                    
                    count = count+1;
                end
                 
                % Write all corrdinates anyway
                LAT(count2) = lat_prof;
                LON(count2) = lon_prof;
                count2 = count2+1;
                
            end

        end %for j
    end %modulo
end %isempty


keyboard