function xtract_latlon(profile_names, gps_files)

% function boot_VMP(file_names, zbin, varargin)
%
% usage ex: xtract_latlon('profile_names', 'gps_files')
%
% file_names is a file containing a list of profile*.mat files that 
% we want to consider.
% In linux, an easy command to do in folder containing *eps*.mat is:
% 
% "ls -1 *profile*.mat | sed 's/\.mat//' > profile_names"
%
% will save a file containing [lat lon] of all profiles
%
% author: F. Cyr - 2011/04/19
%
% ---------------------------------------------------------- %

    
% load *.P files names (file in which are recorded .P files)
fid = fopen(profile_names);
C = textscan(fid, '%s', 'delimiter', '\n');
epsfiles = char(C{1});

no_profile = size(epsfiles, 1); %number of eps_files 

fid = fopen(gps_files);
C = textscan(fid, '%s', 'delimiter', '\n');
gpsfiles = char(C{1});

no_gpsfiles = size(gpsfiles,1); %number of eps_files 


count = 1;
for profile = 1:no_profile
    %profile
    fname = epsfiles(profile, :);
    I = find(fname==' ');   
    fname(I) = [];
    
    load(fname)
    
    % time of profile
    VMP_mtime = mtime(1);   
    
    for i = 1:no_gpsfiles
        gpsname = gpsfiles(i, :);
        I = find(gpsname==' ');   
        gpsname(I) = [];
        
        track = load(gpsname); % [LAT LATmin LON LONmin DD MM YYYY HH MM SS]
        
        latitude = track(:, 1) + track(:,2)/60;  
        longitude = track(:, 3) + track(:,4)/60; 
        n_GPS = datenum(track(:,7), track(:,6), track(:,5), track(:,8), ...
                        track(:,9), track(:,10));
% $$$         if profile==2
% $$$             keyboard
% $$$         end
        
        [Y, I] = min(abs(n_GPS-VMP_mtime));
        
        if i==1
            best_min = Y;
            best_day = n_GPS(I);
            coord(count, :) = [latitude(I), longitude(I)];
        else
            if Y<best_min
                coord(count, :) = [latitude(I), longitude(I)]; 
                best_min=Y;
                best_day = n_GPS(I);

            end
        end
       
    end % loop on GPSfiles
   
    if best_min<1/24
        count = count+1;
    end
    disp(datestr(best_day))
    

end



    
dlmwrite('coord_veryall.dat', coord, 'delimiter',' ','precision',6);
disp(sprintf('%d coordinates saved', count));
